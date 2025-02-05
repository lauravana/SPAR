###########################################
### Main implementation of sparse projected averaged regression (SPAR)
##########################################

# fix errors for auc evaluation

#' Sparse Projected Averaged Regression
#'
#' Apply Sparse Projected Averaged Regression to High-dimensional Data
#' \insertCite{parzer2024lm}{spar} \insertCite{parzer2024glms}{spar}.
#' This function performs the procedure for given thresholds nu and numbers
#' of marginal models, and acts as a help-function for the full cross-validated
#' procedure [spar.cv].
#'
#' @param x n x p numeric matrix of predictor variables.
#' @param y quantitative response vector of length n.
#' @param family  a \link[stats]{family}  object used for the marginal generalized linear model,
#'        default \code{gaussian("identity")}.
#' @param model function creating a "\code{sparmodel}" object; defaults to \code{spar_glmnet()}.
#' @param rp function creating a "\code{randomprojection}" object.
#' @param screencoef function creating a "\code{screeningcoef}" object
#' @param xval optional matrix of predictor variables observations used for
#'        validation of threshold nu and number of models; \code{x} is used
#'        if not provided.
#' @param yval optional response observations used for validation of
#'        threshold nu and number of models; \code{y} is used if not provided.
#' @param nnu number of different threshold values \eqn{\nu} to consider for thresholding;
#'        ignored when nus are given; defaults to 20.
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
#' @param inds optional list of index-vectors corresponding to variables kept
#'  after screening in each marginal model of length \code{max(nummods)};
#'  dimensions need to fit those of RPMs.
#' @param RPMs optional list of projection matrices used in each
#' marginal model of length \code{max(nummods)}, diagonal elements will be
#'  overwritten with a coefficient only depending on the given \code{x} and \code{y}.
#' @param ... further arguments mainly to ensure back-compatibility
#' @returns object of class \code{"spar"} with elements
#' \itemize{
#'  \item \code{betas} p x \code{max(nummods)} sparse matrix of class
#'  \code{"\link[=dgCMatrix-class]{dgCMatrix}"} containing the
#'   standardized coefficients from each marginal model
#'  \item \code{intercepts} used in each marginal model
#'  \item \code{scr_coef} p-vector of coefficients used for screening for standardized predictors
#'  \item \code{inds} list of index-vectors corresponding to variables kept after screening in each marginal model of length max(nummods)
#'  \item \code{RPMs} list of projection matrices used in each marginal model of length \code{max(nummods)}
#'  \item \code{val_res} \code{data.frame} with validation results (validation measure
#'   and number of active variables) for each element of \code{nus} and \code{nummods}
#'  \item \code{val_set} logical flag, whether validation data were provided;
#'  if \code{FALSE}, training data were used for validation
#'  \item \code{nus} vector of \eqn{\nu}'s considered for thresholding
#'  \item \code{nummods} vector of numbers of marginal models considered for validation
#'  \item \code{ycenter} empirical mean of initial response vector
#'  \item \code{yscale} empirical standard deviation of initial response vector
#'  \item \code{xcenter} p-vector of empirical means of initial predictor variables
#'  \item \code{xscale} p-vector of empirical standard deviations of initial predictor variables
#'  \item \code{rp} an object of class "\code{randomprojection}"
#'  \item \code{screencoef} an object of class "\code{screeningcoef}"
#' }
#' @references{
#'   \insertRef{parzer2024lm}{spar}
#'
#'   \insertRef{parzer2024glms}{spar}
#'
#'   \insertRef{Clarkson2013LowRankApprox}{spar}
#'
#'   \insertRef{ACHLIOPTAS2003JL}{spar}
#' }
#' @examples
#' \dontrun{
#' data("example_data")
#' spar_res <- spar(example_data$x,example_data$y,
#' xval=example_data$xtest,yval=example_data$ytest,nummods=c(5,10,15,20,25,30))
#' spar_res
#' coefs <- coef(spar_res)
#' pred <- predict(spar_res,xnew=example_data$x)
#' plot(spar_res)
#' plot(spar_res,"Val_Meas","nummod")
#' plot(spar_res,"Val_numAct","nu")
#' plot(spar_res,"coefs",prange=c(1,400))}
#' @seealso [spar.cv],[coef.spar],[predict.spar],[plot.spar],[print.spar]
#' @export
#' @importFrom stats glm.fit coef fitted gaussian predict rnorm quantile residuals sd var cor glm
#' @importFrom Matrix Matrix solve crossprod tcrossprod rowMeans
#' @importFrom Rdpack reprompt
#' @importFrom rlang list2
#' @importFrom glmnet glmnet
#' @importFrom ROCR prediction performance
#'
spar <- function(x, y,
                 family = gaussian("identity"),
                 model = NULL,
                 rp = NULL,
                 screencoef = NULL,
                 xval = NULL, yval = NULL,
                 nnu = 20, nus = NULL,
                 nummods = c(20),
                 measure = c("deviance","mse","mae","class","1-auc"),
                 inds = NULL, RPMs = NULL,
                 parallel = c("no", "multicore", "snow"),
                 ...) {

  # Set up and checks ----
  measure <- match.arg(measure)
  stopifnot("Length of y does not fit nrow(x)." = length(y) == nrow(x))
  stopifnot("Response y must be numeric." = is.numeric(y))
  # Ensure back compatibility ----
  args <- list(...)
  arg_list <- check_and_set_args(args, family, model,
                                 screencoef, rp,  measure)
  model <- arg_list$model; rp <- arg_list$rp
  screencoef <- arg_list$screencoef; measure <- arg_list$measure

  # Call SPAR algorithm ----
  res <- spar_algorithm(x = x, y = y,
                        family = family,
                        model = model, rp = rp, screencoef = screencoef,
                        xval = xval, yval = yval,
                        nnu = nnu, nus = nus,
                        nummods = nummods,
                        measure = measure,
                        inds = inds, RPMs = RPMs)
  return(res)

}

spar_algorithm <- function(x, y,
                           family, model, rp, screencoef,
                           xval = NULL, yval = NULL,
                           nnu, nus,
                           nummods, measure,
                           inds = NULL, RPMs = NULL,
                           parallel){
  # Start SPAR algorithm
  p <- ncol(x)
  n <- nrow(x)
  # Scaling the x matrix ----
  xcenter <- colMeans(x)
  xscale  <- apply(x, 2, sd)

  if (is.null(inds) || is.null(RPMs)) {
    actual_p <- sum(xscale > 0)
    z <- scale(x[, xscale > 0],
               center = xcenter[xscale > 0],
               scale  = xscale[xscale > 0])
  } else {
    actual_p <- p
    xscale[xscale == 0] <- 1
    z <- scale(x, center = xcenter, scale = xscale)
  }

  # Scaling the y vector ----
  if (family$family=="gaussian" & family$link=="identity") {
    ycenter <- mean(y)
    yscale <- sd(y)
  } else {
    ycenter <- 0
    yscale  <- 1
  }
  yz <- scale(y,center = ycenter,scale = yscale)
  # Setup model ----
  if (is.null(model$control$family))  {
    if (is.null(attr(model, "family"))) {
      model$control$family <- family
    } else {
      model$control$family <- attr(model, "family")
    }
  }
  if (!is.null(model$update_sparmodel)) {
    model <- model$update_sparmodel(model)
  }
  # Setup screening ----
  if (is.null(attr(screencoef, "family"))) {
    attr(screencoef, "family") <- family
  }
  if (!is.null(attr(screencoef, "split_data")) &
      is.null(attr(screencoef, "split_data_prop"))) {
    attr(screencoef, "split_data_prop") <- 0.25
  }
  if (!is.null(attr(screencoef, "split_data_prop"))) {
    attr(screencoef, "split_data") <- TRUE
    scr_inds <- sample(n,
                       ceiling(n * attr(screencoef, "split_data_prop")))  # TODO need to parametrize this
    mar_inds <- seq_len(n)[-scr_inds]
  } else {
    mar_inds <- scr_inds <- seq_len(n)
  }

  if (is.null(attr(screencoef, "nscreen"))) {
    nscreen <- attr(screencoef, "nscreen") <- 2 * n
  } else {
    nscreen <- attr(screencoef, "nscreen")
  }
  mslow <- attr(rp, "mslow")
  if (is.null(mslow)) mslow <- ceiling(log(p))
  msup <- attr(rp, "msup")
  if (is.null(msup)) msup <- ceiling(n/2)
  stopifnot(mslow <= msup)
  stopifnot(msup <= nscreen)

  # Perform screening ----
  if (nscreen < p) {
    scr_coef <- get_screencoef(object = screencoef,
                               x = z[scr_inds,],
                               y = yz[scr_inds, ])
    attr(screencoef, "importance") <- scr_coef

    inc_probs <- abs(scr_coef)
    max_inc_probs <- max(inc_probs)
    inc_probs <- inc_probs/max_inc_probs
    attr(screencoef, "inc_prob") <- inc_probs
  } else {
    scr_coef <- NULL
    message("nscreen >= p. No screening performed.")
  }

  ## Update RP with data only at the beginning if possible, not in each RP! ----
  if (is.null(attr(rp, "family"))) {
    attr(rp, "family") <- family
  }
  if (!is.null(rp$update_data_fun)) {
    rp <- rp$update_data_fun(rp = rp,
                             x = z[scr_inds,],
                             y = yz[scr_inds, ],
                             screencoef = screencoef)
  }

  max_num_mod <- max(nummods)
  intercepts <- numeric(max_num_mod)
  betas_std <- Matrix(data=c(0),actual_p,max_num_mod,sparse = TRUE)

  drawRPMs <- FALSE
  if (is.null(RPMs)) {
    RPMs <- vector("list", length = max_num_mod)
    drawRPMs <- TRUE
    ms <- sample(seq(floor(mslow), ceiling(msup)),
                 max_num_mod, replace=TRUE)
  }

  drawinds <- FALSE
  if (is.null(inds)) {
    inds <- vector("list", length = max_num_mod)
    drawinds <- TRUE
  }

  # SPAR algorithm  ----
  for (i in seq_len(max_num_mod)) {
    ## Screening step  ----
    if (drawinds) {
      if (nscreen < p) {
        if (attr(screencoef, "type") == "fixed") {
          ind_use <- order(inc_probs, decreasing = TRUE)[seq_len(nscreen)]
        }
        if (attr(screencoef, "type") == "prob") {
          ind_use <- sample(actual_p, nscreen, prob=inc_probs)
        }
      } else {
        ind_use <- seq_len(actual_p)
      }
      inds[[i]] <- ind_use
    } else {
      ind_use <- inds[[i]]
    }
    p_use <- length(ind_use)

    ## RP step  ----
    if (drawRPMs) {
      m <- ms[i]
      if (p_use < m) {
        m <- p_use
        RPM <- Matrix::Matrix(diag(1, m),sparse=TRUE)
        RPMs[[i]] <- RPM
      } else {
        RPM    <- get_rp(rp, m = m, included_vector = ind_use)
        RPMs[[i]] <- RPM
      }
    } else {
      RPM <- RPMs[[i]]
      if (!is.null(rp$update_rpm_w_data)) {
        RPM <- rp$update_rpm_w_data(rpm = RPM, rp = rp,
                                    included_vector = ind_use)
      }
    }

    ## Marginal model ----
    znew <- Matrix::tcrossprod(z[mar_inds, ind_use], RPM)

    # if (family$family=="gaussian" & family$link=="identity") {
    #   mar_coef <- tryCatch(solve(crossprod(znew),
    #                              crossprod(znew,yz[mar_inds])),
    #                        error=function(error_message) {
    #                          return(solve(crossprod(znew)+0.01*diag(ncol(znew)),
    #                                       crossprod(znew,yz[mar_inds])))
    #                        })
    #   intercepts[i] <- 0
    #   betas_std[ind_use,i] <- Matrix::crossprod(RPM,mar_coef)
    # } else {
    #   glmnet_res <- glmnet(znew,y[mar_inds],
    #                        family = fit_family, alpha=0)
    #   mar_coef <- coef(glmnet_res, s = min(glmnet_res$lambda))
    #   intercepts[i] <- mar_coef[1]
    #   betas_std[ind_use,i] <- crossprod(RPM,mar_coef[-1])
    # }
    res <- model$model_fun(yz[mar_inds], znew, object = model)
    intercepts[i] <- res$intercept
    betas_std[ind_use,i] <- crossprod(RPM, res$gammas)

  }

  if (is.null(nus)) {
    if (nnu>1) {
      nus <- unname(c(0, quantile(abs(betas_std@x),
                           probs=seq_len(nnu-1)/(nnu-1))))
    } else {
      nus <- 0
    }
  } else {
    nnu <- length(nus)
  }

  ## Validation set
  val_res <- data.frame(nnu = NULL, nu = NULL,
                        nummod = NULL,numAct = NULL, Meas = NULL)
  if (!is.null(yval) && !is.null(xval)) {
    val_set <- TRUE
  } else {
    val_set <- FALSE
    yval <- y
    xval <- x
  }

  val.meas <- get_val_measure_function(measure, family)

  tabnummodres <- lapply(nummods,  function(nummod) {
    coef <- betas_std[,seq_len(nummod),drop=FALSE]
    abscoef <- abs(coef)
    tabres <- lapply(seq_len(nnu), function(l){
      thresh <- nus[l]
      tmp_coef <- coef
      tmp_coef[abscoef<thresh] <- 0

      avg_coef <- rowMeans(tmp_coef)
      tmp_beta <- numeric(p)
      tmp_beta[xscale>0] <- yscale*avg_coef/(xscale[xscale>0])
      tmp_intercept <- mean(intercepts[seq_len(nummod)]) +
        drop(ycenter - sum(xcenter*tmp_beta) )
      eta_hat <- (xval %*% tmp_beta) + tmp_intercept

      c(nnu = l,
        nu = unname(thresh),
        nummod = nummod,
        numAct = sum(tmp_beta!=0),
        Meas = val.meas(yval,eta_hat)
      )
    })
    out <- do.call("rbind", tabres)
    colnames(out) <- c("nnu","nu","nummod","numAct","Meas")
    out
  })
  val_res <- do.call("rbind.data.frame", tabnummodres)
  betas <- Matrix(data=c(0),p,max_num_mod,sparse = TRUE)
  betas[xscale>0,] <- betas_std


  ## Clean up
  family_str <- paste0(family$family, "(", family$link, ")")
  # attr(rp,"family") <- paste0(attr(rp,"family")$family, "(",
  #                         attr(rp,"family")$link, ")")
  # attr(screencoef,"family") <- paste0(attr(screencoef,"family")$family, "(",
  #                             attr(screencoef,"family")$link, ")")
  res <- list(betas = betas, intercepts = intercepts,
              scr_coef = scr_coef,
              inds = inds, RPMs = RPMs,
              val_res = val_res, val_set = val_set,
              nus = nus, nummods = nummods,
              ycenter = ycenter, yscale = yscale,
              xcenter = xcenter, xscale = xscale,
              family = family_str,
              measure = measure,
              rp = rp,
              screencoef = screencoef
  )

  attr(res,"class") <- "spar"

  return(res)
}
#' coef.spar
#'
#' Extract coefficients from spar object
#' @param object result of spar function of class "spar".
#' @param nummod number of models used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param nu threshold level used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param ... further arguments passed to or from other methods
#' @return List of coefficients with elements
#' \itemize{
#'  \item intercept
#'  \item beta
#'  \item nummod
#'  \item nu
#' }
#' @export

coef.spar <- function(object,
                      nummod = NULL,
                      nu = NULL, ...) {
  if (is.null(nummod) & is.null(nu)) {
    best_ind <- which.min(object$val_res$Meas)
    par <- object$val_res[best_ind,]
    nummod <- par$nummod
    nu <- par$nu
  } else if (is.null(nummod)) {
    if (!nu %in% object$val_res$nu) {
      stop("Nu needs to be among the previously fitted values when nummod is not provided!")
    }
    tmp_val_res <- object$val_res[object$val_res$nu==nu,]
    nummod <- tmp_val_res$nummod[which.min(tmp_val_res$Meas)]
  } else if (is.null(nu)) {
    if (!nummod %in% object$val_res$nummod) {
      stop("Number of models needs to be among the previously fitted values when nu is not provided!")
    }
    tmp_val_res <- object$val_res[object$val_res$nummod==nummod,]
    nu <- tmp_val_res$nu[which.min(tmp_val_res$Meas)]
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

#' predict.spar
#'
#' Predict responses for new predictors from spar object
#' @param object result of spar function of class  \code{"spar"}.
#' @param xnew matrix of new predictor variables; must have same number of columns as  \code{x}.
#' @param type the type of required predictions; either on response level (default) or on link level
#' @param avg_type type of averaging the marginal models; either on link (default) or on response level
#' @param nummod number of models used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param nu threshold level used to form coefficients; value with minimal validation Meas is used if not provided.
#' @param coef optional; result of  \code{coef.spar}, can be used if  \code{coef.spar()} has already been called.
#' @param ... further arguments passed to or from other methods
#' @return Vector of predictions
#' @export

predict.spar <- function(object,
                         xnew,
                         type = c("response","link"),
                         avg_type = c("link","response"),
                         nummod = NULL,
                         nu = NULL,
                         coef = NULL, ...) {
  if (ncol(xnew)!=length(object$xscale)) {
    stop("xnew must have same number of columns as initial x!")
  }
  type <- match.arg(type)
  avg_type <- match.arg(avg_type)
  if (is.null(coef)) {
    coef <- coef(object,nummod,nu)
  }
  object$family <- eval(parse(text = object$family))
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
      final_coef <- object$betas[object$xscale>0,1:coef$nummod,drop=FALSE]
      final_coef[abs(final_coef)<coef$nu] <- 0

      preds <- sapply(1:coef$nummod,function(j){
        tmp_coef <- final_coef[,j]
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

#' plot.spar
#'
#' Plot errors or number of active variables over different thresholds or number of models of spar result, or residuals vs fitted
#' @param x result of spar function of class  \code{"spar"}.
#' @param plot_type one of  \code{c("Val_Measure","Val_numAct","res-vs-fitted","coefs")}.
#' @param plot_along one of \code{c("nu","nummod")}; ignored when  \code{plot_type="res-vs-fitted"}.
#' @param nummod fixed value for number of models when  \code{plot_along="nu"}
#'  for  \code{plot_type="Val_Measure"} or  \code{"Val_numAct"}; same as for \code{\link{predict.spar}} when  \code{plot_type="res-vs-fitted"}.
#' @param nu fixed value for \eqn{\nu} when  \code{plot_along="nummod"} for
#'  \code{plot_type="Val_Measure"} or  \code{"Val_numAct"}; same as for \code{\link{predict.spar}} when  \code{plot_type="res-vs-fitted"}.
#' @param xfit data used for predictions in  \code{"res-vs-fitted"}.
#' @param yfit data used for predictions in  \code{"res-vs-fitted"}.
#' @param prange optional vector of length 2 for  \code{"coefs"}-plot to give
#'  the limits of the predictors' plot range; defaults to  \code{c(1, p)}.
#' @param coef_order optional index vector of length p for "coefs"-plot to give
#'  the order of the predictors; defaults to  \code{1 : p}.
#' @param digits number of significant digits to be displayed in the axis; defaults to 2L.
#' @param ... further arguments passed to or from other methods
#' @return ggplot2 object
#' @import ggplot2
#' @export
plot.spar <- function(x,
                      plot_type = c("Val_Measure","Val_numAct","res-vs-fitted","coefs"),
                      plot_along = c("nu","nummod"),
                      nummod = NULL,
                      nu = NULL,
                      xfit = NULL,
                      yfit = NULL,
                      prange = NULL,
                      coef_order = NULL,
                      digits = 2L, ...) {
  spar_res <- x
  plot_type <- match.arg(plot_type)
  plot_along <- match.arg(plot_along)
  mynummod <- nummod
  if (plot_type == "res-vs-fitted") {
    if (is.null(xfit) | is.null(yfit)) {
      stop("xfit and yfit need to be provided for res-vs-fitted plot!")
    }
    pred <- predict(spar_res, xfit, nummod, nu, type = "response")
    res <- ggplot2::ggplot(data = data.frame(fitted=.data$pred,
                                             residuals=.data$yfit-.data$pred),
                           ggplot2::aes(x=.data$fitted,y=.data$residuals)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=0.5)
  } else if (plot_type == "Val_Measure") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }

      tmp_df <- subset(spar_res$val_res,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x=.data$nnu,y=.data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        # ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),labels=round(spar_res$val_res$nu,3)) +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(tmp_df)),
                                    labels=formatC(tmp_df$nu,
                                                   format = "e", digits = digits)) +
        ggplot2::labs(x=expression(nu),y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nnu[ind_min],
                                            y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red") +
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(nu)) {
        nu <- spar_res$val_res$nu[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- subset(spar_res$val_res,nu==nu)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,
                             ggplot2::aes(x=.data$nummod,y=.data$Meas)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::labs(y=spar_res$measure) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nummod[ind_min],y=tmp_df$Meas[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),list(txt=tmp_title,v=round(nu,3))))
    }
  } else if (plot_type=="Val_numAct") {
    if (plot_along=="nu") {
      if (is.null(nummod)) {
        mynummod <- spar_res$val_res$nummod[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal nummod="
      } else {
        tmp_title <- "Fixed given nummod="
      }
      tmp_df <- subset(spar_res$val_res,nummod==mynummod)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nnu,y=.data$numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        # ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),1),labels=round(spar_res$val_res$nu,3)) +
        ggplot2::scale_x_continuous(breaks=seq(1,nrow(spar_res$val_res),2),
                                    labels=formatC(spar_res$val_res$nu[seq(1,nrow(spar_res$val_res),2)],
                                                   format = "e", digits = digits)) +
        ggplot2::labs(x=expression(nu)) +
        ggplot2::geom_point(data=data.frame(x=tmp_df$nnu[ind_min],y=tmp_df$numAct[ind_min]),
                            ggplot2::aes(x=.data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(paste0(tmp_title,mynummod))
    } else {
      if (is.null(nu)) {
        nu <- spar_res$val_res$nu[which.min(spar_res$val_res$Meas)]
        tmp_title <- "Fixed optimal "
      } else {
        tmp_title <- "Fixed given "
      }
      tmp_df <- subset(spar_res$val_res,nu==nu)
      ind_min <- which.min(tmp_df$Meas)

      res <- ggplot2::ggplot(data = tmp_df,ggplot2::aes(x=.data$nummod,y=.data$numAct)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(
          data=data.frame(x=tmp_df$nummod[ind_min],
                          y=tmp_df$numAct[ind_min]),
          ggplot2::aes(x = .data$x,y=.data$y),col="red")+
        ggplot2::ggtitle(substitute(paste(txt,nu,"=",v),
                                    list(txt=tmp_title,v=round(nu,3))))
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
    tmp_df <- tidyr::pivot_longer(tmp_mat,cols = (1:nummod),
                                  names_to = "marginal model",values_to = "value")

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

#' print.spar
#'
#' Print summary of spar result
#' @param x result of  \code{spar()} function of class  \code{"spar"}.
#' @param ... further arguments passed to or from other methods
#' @return text summary
#' @export
print.spar <- function(x, ...) {
  mycoef <- coef(x)
  beta <- mycoef$beta
  cat(sprintf("spar object:\nSmallest Validation Measure reached for nummod=%d,
              nu=%s leading to %d / %d active predictors.\n",
              mycoef$nummod, formatC(mycoef$nu,digits = 2,format = "e"),
              sum(beta!=0),length(beta)))
  cat("Summary of those non-zero coefficients:\n")
  print(summary(beta[beta!=0]))
}
