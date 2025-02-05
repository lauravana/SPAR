
#' Sparse Projected Averaged Regression
#'
#' Apply Sparse Projected Averaged Regression to High-dimensional Data (see Parzer, Vana-Guer and Filzmoser 2023).
#'
#' @param x n x p numeric matrix of predictor variables.
#' @param y quantitative response vector of length n.
#' @param family  a "\code{\link[stats]{family}}" object used for the marginal generalized linear model,
#'        default \code{gaussian("identity")}.
#' @param model function creating a "\code{sparmodel}" object; defaults to \code{spar_glmnet()}.
#' @param rp function creating a "\code{randomprojection}" object.
#' @param screencoef function creating a "\code{screeningcoef}" object
#' @param nfolds number of folds to use for cross-validation; should be at least 2, defaults to 10.
#' @param nnu number of different threshold values \eqn{\nu} to consider for thresholding;
#'        ignored when \code{nus} is provided; defaults to 20.
#' @param nus optional vector of \eqn{\nu}'s to consider for thresholding;
#'         if not provided, \code{nnu} values ranging from 0 to the maximum absolute
#'         marginal coefficient are used.
#' @param nummods vector of numbers of marginal models to consider for
#'        validation; defaults to \code{c(20)}.
#' @param measure loss to use for validation; defaults to \code{"deviance"}
#'        available for all families. Other options are \code{"mse"} or \code{"mae"}
#'         (between responses and predicted means, for all families),
#'         \code{"class"} (misclassification error) and
#'         \code{"1-auc"} (one minus area under the ROC curve) both just for
#'         binomial family.
# #' @param type.rpm  type of random projection matrix to be employed;
# #'        one of \code{"cwdatadriven"},
# #'        \code{"cw"} \insertCite{Clarkson2013LowRankApprox}{spar},
# #'        \code{"gaussian"}, \code{"sparse"} \insertCite{ACHLIOPTAS2003JL}{spar};
# #'        defaults to \code{"cwdatadriven"}.
# #' @param type.screening  type of screening coefficients; one of \code{"ridge"},
# #'        \code{"marglik"}, \code{"corr"}; defaults to \code{"ridge"} which is
# #'        based on the ridge coefficients where the penalty converges to zero.
#' @param ... further arguments mainly to ensure back-compatibility
#' @returns object of class \code{"spar.cv"} with elements
#' \itemize{
#'  \item \code{betas} p x  \code{max(nummods)} sparse matrix of class
#'  \code{"\link[=dgCMatrix-class]{dgCMatrix}"} containing the
#'   standardized coefficients from each marginal model
#'  \item \code{intercepts} used in each marginal model, vector of length \code{max(nummods)}
#'  \item \code{scr_coef} p-vector of coefficients used for screening for standardized predictors
#'  \item \code{inds} list of index-vectors corresponding to variables kept after
#'  screening in each marginal model of length  \code{max(nummods)}
#'  \item \code{RPMs} list of projection matrices used in each marginal model of length \code{max(nummods)}
#'  \item \code{val_sum} \code{data.frame} with CV results (mean and sd validation measure and mean number of active variables) for each element of nus and nummods
#'  \item \code{nus} vector of \eqn{\nu}'s considered for thresholding
#'  \item \code{nummods} vector of numbers of marginal models considered for validation
#'  \item \code{ycenter} empirical mean of initial response vector
#'  \item \code{yscale} empirical standard deviation of initial response vector
#'. \item \code{xcenter} p-vector of empirical means of initial predictor variables
#'  \item \code{xscale} p-vector of empirical standard deviations of initial predictor variables
#'  \item \code{rp} an object of class "\code{randomprojection}"
#'  \item \code{screencoef} an object of class "\code{screeningcoef}"
#' }
#' @examples
#' \dontrun{
#' data("example_data")
#' spar_res <- spar.cv(example_data$x,example_data$y,nummods=c(5,10,15,20,25,30))
#' spar_res
#' coefs <- coef(spar_res)
#' pred <- predict(spar_res,example_data$x)
#' plot(spar_res)
#' plot(spar_res,"Val_Meas","nummod")
#' plot(spar_res,"Val_numAct","nu")
#' plot(spar_res,"coefs",prange=c(1,400))}
#' @seealso [spar],[coef.spar.cv],[predict.spar.cv],[plot.spar.cv],[print.spar.cv]
#' @export
spar.cv <- function(x, y,
                    family = gaussian("identity"),
                    model = spar_glmnet(),
                    rp = NULL,
                    screencoef = NULL,
                    nfolds = 10,
                    nnu = 20,
                    nus = NULL,
                    nummods = c(20),
                    measure = c("deviance","mse","mae","class","1-auc"),
                    ...
) {
  # Set up and checks ----
  n <- length(y)
  stopifnot("Length of y does not fit nrow(x)." = n == nrow(x))
  stopifnot("Response y must be numeric." = is.numeric(y))
  stopifnot("nfolds must be at least 2." = nfolds >= 2)
  stopifnot("matrix" %in% class(x) |"data.frame" %in% class(x))
  x <- as.matrix(x)
  if (!is.numeric(x[1,1])) {
    stop("There are non-numeric data entries, numerical matrix needed!")
  }

  measure <- match.arg(measure)

  # Ensure back compatibility ----
  args <- list(...)
  arg_list <- check_and_set_args(args, family, model,
                                 screencoef, rp,  measure)
  model <- arg_list$model; rp <- arg_list$rp
  screencoef <- arg_list$screencoef; measure <- arg_list$measure

  # Run initial spar algorithm ----
  SPARres <- spar_algorithm(x = x, y = y,
                        family = family,
                        model = model, rp = rp, screencoef = screencoef,
                        xval = NULL, yval = NULL,
                        nnu = nnu, nus = nus,
                        nummods = nummods,
                        measure = measure,
                        inds = NULL, RPMs = NULL)
  # SPARres <- spar(x, y, family = family, model = model,
  #                 rp = rp,
  #                 screencoef = screencoef,
  #                 nnu = nnu,
  #                 nummods = nummods,
  #                 measure = measure, ...)

  val_res <- SPARres$val_res
  folds <- sample(cut(seq_len(n), breaks = nfolds, labels=FALSE))
  for (k in seq_len(nfolds)) {
    fold_id <- (folds == k)
    foldSPARres <- spar_algorithm(x[!fold_id,SPARres$xscale>0],y[!fold_id],
                        family = family, model = model,
                        xval = x[fold_id,SPARres$xscale>0],
                        yval = y[fold_id],
                        rp = rp, screencoef = screencoef,
                        nus = SPARres$nus,
                        inds = SPARres$inds,
                        RPMs = SPARres$RPMs,
                        nummods = nummods,
                        measure = measure)
    val_res <- rbind(val_res,foldSPARres$val_res)
  }

  val_sum <- dplyr::group_by(val_res, .data$nnu, .data$nu, .data$nummod)
  suppressMessages(
    val_sum <-  dplyr::summarise(val_sum,
                         mMeas = mean(.data$Meas,na.rm=TRUE),
                         sdMeas = sd(.data$Meas,na.rm=TRUE),
                         mNumAct = mean(.data$numAct,na.rm=TRUE))
  )

  res <- list(betas = SPARres$betas, intercepts = SPARres$intercepts,
              scr_coef = SPARres$scr_coef, inds = SPARres$inds,
              RPMs = SPARres$RPMs,
              val_sum = val_sum, nus = SPARres$nus,
              nummods=nummods,
              family = family,
              measure = measure,
              rp = rp, screencoef = screencoef,
              ycenter = SPARres$ycenter, yscale = SPARres$yscale,
              xcenter = SPARres$xcenter, xscale = SPARres$xscale)

  attr(res,"class") <- "spar.cv"
  return(res)
}

#' coef.spar.cv
#'
#' Extract coefficients from spar object
#' @param object result of spar.cv function of class "spar.cv".
#' @param opt_par one of c("1se","best"), chooses whether to select the best
#' pair of nus and nummods according to CV-Meas, or the sparsest solution within
#' one sd of that optimal CV-Meas;
#' ignored when nummod and nu are given
#' @param nummod optional number of models used to form coefficients
#' @param nu optional threshold level used to form coefficients
#' @param ... further arguments passed to or from other methods
#' @return List of coefficients with elements
#' \itemize{
#'  \item intercept
#'  \item beta
#'  \item nummod
#'  \item nu
#' }
#' @export

coef.spar.cv <- function(object,
                         opt_par = c("best","1se"),
                         nummod = NULL,
                         nu = NULL, ...) {
  opt_nunum <- match.arg(opt_par)
  if (is.null(nummod) & is.null(nu)) {
    best_ind <- which.min(object$val_sum$mMeas)
    if (opt_nunum=="1se") {
      allowed_ind <- object$val_sum$mMeas<object$val_sum$mMeas[best_ind]+
        object$val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(object$val_sum$mNumAct[allowed_ind])
      par <- object$val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- object$val_sum[best_ind,]
    }
    nummod <- par$nummod
    nu <- par$nu
  } else if (is.null(nummod)) {
    if (!nu %in% object$val_sum$nu) {
      stop("nu needs to be among the previously fitted values when nummod is not provided!")
    }
    tmp_val_sum <- object$val_sum[object$val_sum$nu==nu,]
    if (opt_nunum=="1se") {
      allowed_ind <- tmp_val_sum$mMeas<tmp_val_sum$mMeas[best_ind]+tmp_val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mMeas),]
    }
    nummod <- par$nummod
  } else if (is.null(nu)) {
    if (!nummod %in% object$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when nu is not provided!")
    }
    tmp_val_sum <- object$val_sum[object$val_sum$nummod==nummod,]
    if (opt_nunum=="1se") {
      allowed_ind <- tmp_val_sum$mMeas<tmp_val_sum$mMeas[best_ind]+tmp_val_sum$sdMeas[best_ind]
      ind_1cv <- which.min(tmp_val_sum$mNumAct[allowed_ind])
      par <- tmp_val_sum[allowed_ind,][ind_1cv,]
    } else {
      par <- tmp_val_sum[which.min(tmp_val_sum$mMeas),]
    }
    nu <- par$nu
  } else {
    if (length(nummod)!=1 | length(nu)!=1) {
      stop("Length of nummod and nu must be 1!")
    }
  }

  if (nummod > ncol(object$betas)) {
    warning("Number of models is too high, maximum of fitted is used instead!")
    nummod <- ncol(object$betas)
  }

  # calc for chosen parameters
  final_coef <- object$betas[object$xscale>0,1:nummod,drop=FALSE]
  final_coef[abs(final_coef)<nu] <- 0
  p <- length(object$xscale)
  beta <- numeric(p)
  beta[object$xscale>0] <- object$yscale*Matrix::rowMeans(final_coef)/(object$xscale[object$xscale>0])
  intercept <- object$ycenter + mean(object$intercepts[1:nummod]) - sum(object$xcenter*beta)
  return(list(intercept=intercept,beta=beta,nummod=nummod,nu=nu))
}

#' predict.spar.cv
#'
#' Predict responses for new predictors from spar object
#' @param object result of spar function of class "spar".
#' @param xnew matrix of new predictor variables; must have same number of columns as  \code{x}.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param opt_par one of  \code{c("best","1se")}, chooses whether to select the
#'  best pair of  \code{nus} and  \code{nummods} according to CV-Meas, or the
#'  sparsest solution within one sd of that optimal CV-Meas;
#'  ignored when  \code{nummod} and  \code{nu}, or  \code{coef} are given
#' @param nummod number of models used to form coefficients; value with
#' minimal validation  \code{Meas} is used if not provided.
#' @param nu threshold level used to form coefficients; value with minimal
#'  validation  \code{Meas} is used if not provided.
#' @param coef optional; result of \code{\link{coef.spar.cv}}, can be used if
#'  \code{\link{coef.spar.cv}} has already been called.
#' @param ... further arguments passed to or from other methods
#' @return Vector of predictions
#' @export
predict.spar.cv <- function(object,
                            xnew,
                            type = c("response","link"),
                            avg_type = c("link","response"),
                            opt_par = c("best","1se"),
                            nummod = NULL,
                            nu = NULL,
                            coef = NULL, ...) {
  if (ncol(xnew)!=length(object$xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  if (is.null(coef)) {
    coef <- coef(object,opt_par,nummod,nu)
  }
  if (avg_type=="link") {
    if (type=="link") {
      res <- as.numeric(xnew%*%coef$beta + coef$intercept)
    } else {
      eta <- as.numeric(xnew%*%coef$beta + coef$intercept)
      res <- object$family$linkinv(eta)
    }
  } else {
    if (type=="link") {
      res <- as.numeric(xnew%*%coef$beta + coef$intercept)
    } else {
      # do diff averaging
      final_coef <- object$betas[,1:coef$nummod,drop=FALSE]
      final_coef[abs(final_coef)<coef$nu] <- 0

      preds <- sapply(1:coef$nummod,function(j){
        tmp_coef <- final_coef[object$xscale>0,j]
        beta <- numeric(length(object$xscale))
        beta[object$xscale>0] <- object$yscale*tmp_coef/(object$xscale[object$xscale>0])
        intercept <- object$ycenter + object$intercepts[j]  - sum(object$xcenter*beta)
        eta <- as.numeric(xnew%*%beta + coef$intercept)
        object$family$linkinv(eta)
      })
      res <- rowMeans(preds)
    }
  }
  return(res)
}

#' plot.spar.cv
#'
#' Plot errors or number of active variables over different thresholds or number of models of spar.cv result, or residuals vs fitted
#' @param x result of spar.cv function of class  \code{"spar.cv"}.
#' @param plot_type one of  \code{c("Val_Measure","Val_numAct","res-vs-fitted","coefs")}.
#' @param plot_along one of  \code{c("nu","nummod")}; ignored when  \code{plot_type="res-vs-fitted"}.
#' @param opt_par one of  \code{c("1se","best")}, chooses whether to select the
#'  best pair of  \code{nus} and  \code{nummods} according to CV-Meas, or the
#'  sparsest solution within one sd of that optimal CV-Meas;
#' ignored when  \code{nummod} and  \code{nu}, or  \code{coef} are given
#' @param nummod fixed value for  \code{nummod} when  \code{plot_along="nu"} for
#'  \code{plot_type="Val_Measure"} or  \code{"Val_numAct"};
#'  same as for \code{\link{predict.spar.cv}} when plot_type="res-vs-fitted".
#' @param nu fixed value for \eqn{\nu} when  \code{plot_along="nummod"}
#' for  \code{plot_type="Val_Measure"} or  \code{"Val_numAct"}; same as for \code{\link{predict.spar.cv}} when  \code{plot_type="res-vs-fitted"}.
#' @param xfit data used for predictions in  \code{"res-vs-fitted"}.
#' @param yfit data used for predictions in  \code{"res-vs-fitted"}.
#' @param opt_par one of  \code{c("best","1se")}, only needed for
#'  \code{plot_type="res-vs-fitted"} to set type of predictions, see \code{\link{predict.spar.cv}}.
#' @param prange optional vector of length 2 for  \code{"coefs"}-plot to give the limits of the predictors' plot range; defaults to  \code{c(1, p)}.
#' @param coef_order optional index vector of length p for \code{"coefs"}-plot to give the order of the predictors; defaults to  \code{1 : p}.
#' @param digits number of significant digits to be displayed in the axis; defaults to 2L.
#' @param ... further arguments passed to or from other methods
#' @return ggplot2 object
#' @import ggplot2
#' @export
plot.spar.cv <- function(x,
                         plot_type = c("Val_Measure","Val_numAct","res-vs-fitted","coefs"),
                         plot_along = c("nu","nummod"),
                         nummod = NULL,
                         nu = NULL,
                         xfit = NULL,
                         yfit = NULL,
                         opt_par = c("best","1se"),
                         prange = NULL,
                         coef_order = NULL, digits = 2, ...) {
  spar_res <- x
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  opt_par <- match.arg(opt_par)
  mynummod <- nummod
  my_val_sum <- spar_res$val_sum
  colnames(my_val_sum)[match(c("mMeas", "mNumAct"),colnames(my_val_sum))] <- c("Meas", "numAct")
  # my_val_sum <- dplyr::rename(spar_res$val_sum,
  #                             Meas="mMeas",
  #                             numAct="mNumAct")

  if (plot_type=="res-vs-fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res-vs-fitted plot!")
    }
    pred <- predict(spar_res,xfit,opt_par=opt_par,nummod=nummod,nu=nu)
    res <- ggplot2::ggplot(data = data.frame(fitted=pred,residuals=yfit-pred),
                           ggplot2::aes(x=.data$fitted,y=.data$residuals)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  } else if (plot_type=="Val_Measure") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- subset(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas < (tmp_df$Meas+tmp_df$sdMeas)[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x = .data$nnu,y = .data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        # ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),1),labels=round(my_val_sum$nu,3)) +
        ggplot2::scale_x_continuous(
          breaks=seq(1,nrow(tmp_df)),
          labels=formatC(tmp_df$nu,
                         format = "e", digits = digits)) +
        ggplot2::labs(x=expression(nu),y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nnu[ind_min],
                                            y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$Meas-.data$sdMeas,
                                          ymax=.data$Meas+.data$sdMeas),
                             alpha=0.2,linetype=2,show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                   color=2,show.legend = FALSE,
                   data=data.frame(x = c(tmp_df$nnu[ind_min],tmp_df$nnu[allowed_ind][ind_1se]),
                                   y = c(tmp_df$Meas[ind_min],tmp_df$Meas[allowed_ind][ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nnu[ind_min],
                          y = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          xend = tmp_df$nnu[allowed_ind][ind_1se],
                          yend = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          color=2,linetype=2)
    } else {
      if (is.null(nu)) {
        nu <- my_val_sum$nu[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- subset(my_val_sum,nu==nu)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nummod,y=.data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),list(txt=tmp_title,v=round(nu,3)))) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$Meas-.data$sdMeas,
                                          ymax=.data$Meas+.data$sdMeas),
                             alpha=0.2,linetype=2,show.legend = FALSE)+
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$Meas[ind_min],tmp_df$Meas[allowed_ind][ind_1se]))) +
        ggplot2::annotate("segment",x = tmp_df$nummod[ind_min],
                          y = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          xend = tmp_df$nummod[allowed_ind][ind_1se],
                          yend = tmp_df$Meas[ind_min] + tmp_df$sdMeas[ind_min],
                          color=2,linetype=2)
    }
  } else if (plot_type=="Val_numAct") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- my_val_sum$nummod[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- subset(my_val_sum,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nnu,y=.data$numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        # ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),1),labels=round(my_val_sum$nu,3)) +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(my_val_sum),2),
                                    labels=formatC(my_val_sum$nu[seq(1,nrow(my_val_sum),2)],
                                                   format = "e", digits = digits)) +
        ggplot2::labs(x=expression(nu)) +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nnu[ind_min],tmp_df$nnu[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numAct[ind_min],tmp_df$numAct[allowed_ind][ind_1se]))) +
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(nu)) {
        nu <- my_val_sum$nu[which.min(my_val_sum$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- subset(my_val_sum,nu==nu)
      ind_min <- which.min(tmp_df$Meas)

      allowed_ind <- tmp_df$Meas<tmp_df$Meas[ind_min]+tmp_df$sdMeas[ind_min]
      ind_1se <- which.min(tmp_df$numAct[allowed_ind])

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nummod,y=.data$numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
                            color=2,show.legend = FALSE,
                            data=data.frame(x = c(tmp_df$nummod[ind_min],tmp_df$nummod[allowed_ind][ind_1se]),
                                            y = c(tmp_df$numAct[ind_min],tmp_df$numAct[allowed_ind][ind_1se]))) +
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),list(txt=tmp_title,v=round(nu,3))))

    }
  } else if (plot_type=="coefs") {
    p <- nrow(spar_res$betas)
    nummod <- ncol(spar_res$betas)
    if (is.null(prange)) {
      prange <- c(1,p)
    }
    if (is.null(coef_order)) {
      coef_order <- 1:p
    }


    tmp_mat <- data.frame(t(apply(as.matrix(spar_res$betas)[coef_order,],1,
                                  function(row)row[order(abs(row),decreasing = TRUE)])),
                          predictor=1:p)
    colnames(tmp_mat) <- c(1:nummod,"predictor")
    tmp_df <- tidyr::pivot_longer(tmp_mat, cols = (1:nummod),
                                  names_to = "marginal model",
                                  values_to = "value")
    tmp_df$`marginal model` <- as.numeric(tmp_df$`marginal model`)

    mrange <- max(Matrix::rowSums(spar_res$betas != 0))

    res <- ggplot2::ggplot(tmp_df,ggplot2::aes(x=.data$predictor,
                                               y=.data$`marginal model`,
                                               fill=.data$value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2() +
      ggplot2::coord_cartesian(xlim=prange,ylim=c(1,mrange)) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.border = ggplot2::element_blank())

  } else {
    res <- NULL
  }
  return(res)
}


#' print.spar.cv
#'
#' Print summary of spar.cv result
#' @param x result of  \code{spar.cv()} function of class  \code{"spar.cv"}.
#' @param ... further arguments passed to or from other methods
#' @return text summary
#' @export
print.spar.cv <- function(x, ...) {
  spar_res <- x
  mycoef_best <- coef(spar_res,opt_par = "best")
  mycoef_1se <- coef(spar_res,opt_par = "1se")
  cat(sprintf("spar.cv object:\nSmallest CV-Meas %.1f reached for nummod=%d,
              nu=%s leading to %d / %d active predictors.\n",
              min(spar_res$val_sum$mMeas),mycoef_best$nummod,
              formatC(mycoef_best$nu,digits = 2,format = "e"),
              sum(mycoef_best$beta!=0),length(mycoef_best$beta)))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_best$beta[mycoef_best$beta!=0]))
  cat(sprintf("\nSparsest coefficient within one standard error of best CV-Meas
              reached for nummod=%d, nu=%s \nleading to %d / %d active
              predictors with CV-Meas %.1f.\n",
              mycoef_1se$nummod,
              formatC(mycoef_1se$nu,digits = 2,format = "e"),
              sum(mycoef_1se$beta!=0),length(mycoef_1se$beta),
              spar_res$val_sum$mMeas[spar_res$val_sum$nummod==mycoef_1se$nummod
                                     & spar_res$val_sum$nu==mycoef_1se$nu]))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(mycoef_1se$beta[mycoef_1se$beta!=0]))
}
