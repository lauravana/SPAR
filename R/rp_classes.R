#' Constructor function for building randomprojection objects
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param name character
#' @param generate_fun function for generating the random projection matrix. This
#' function should have with arguments \code{rp}, which is a randomprojection
#' object, \code{m}, the target dimension and a vector of indexes
#' \code{included_vector} which shows the column index of the original variables in the
#' \code{x} matrix to be projected using the random projection. This is needed
#' due to the fact that screening can be employed pre-projection.
#' @param update_data_fun function for updating the randomprojection object with
#' information from the data. This
#' function should have arguments \code{rp}, which is a randomprojection
#' object and `x` (the matrix of predictors)
#' and `y` (the vector of responses).
#' @param update_rpm_w_data function for updating the random projection matrix with data.
#' This can be used for the case where a list of random projection matrices is
#' provided by argument \code{RPMs}. In this case, the random structure is kept
#' fixed, but the data-dependent part gets updated with the provided data. Defaults
#' to NULL. If not provided, the values of the provided RPMs do not change.
#' @param control list of controls for random projection. Can include minimum and
#' maximum dimension for the projection defaults to
#' \code{list(mslow = NULL, msup = NULL)}
#' @return a function which in turn creates an object of class randomprojection
#'
#' @export
constructor_randomprojection <- function(name, generate_fun,
                                         update_data_fun = NULL,
                                         update_rpm_w_data = NULL,
                                         control = list()) {
  ## Checks
  stopifnot(names(formals(generate_fun)) %in% c("rp", "m", "included_vector"))
  if (!is.null(update_data_fun)) {
    stopifnot("Function update_data_fun should have arguments rp, x, y, screencoef."=names(formals(update_data_fun)) %in% c("rp", "x", "y", "screencoef"))
  }
  if (!is.null(update_rpm_w_data)) {
    stopifnot(
      "Function update_rpm_w_data should have arguments rpm,  rp, included_vector."=names(formals(update_rpm_w_data)) %in% c("rpm", "rp", "included_vector"))
  }
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                generate_fun = generate_fun,
                update_data_fun = update_data_fun,
                update_rpm_w_data = update_rpm_w_data,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    if (is.null(attr(out, "data"))) {
      attr(out, "data") <- ifelse(is.null(out$update_data_fun),
                                  FALSE, TRUE)
    }
    class(out) <- c("randomprojection")
    return(out)
  }
}

#' Function which works on all random projection objects
#' @param rp an object of class randomprojection
#' @param m integer goal dimension of the projection
#' @param included_vector a vector containing column index of the original variables in the
#' \code{x} matrix to be projected using the random projection
#'
#' @keywords internal
get_rp <- function(rp, m, included_vector) {
  RM <- rp$generate_fun(rp, m, included_vector)
  return(RM)
}

#'
#' Gaussian random projection matrix
#'
#' @param rp object of class  "\code{randomprojection}"
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return matrix with m rows and
#'  \code{length(included_vector)} columns sampled from the normal distribution.
#' @keywords internal
generate_gaussian <- function(rp, m, included_vector) {
  p <- length(included_vector)
  control_rnorm <-
    rp$control[names(rp$control)  %in% names(formals(rnorm))]
  vals <- do.call(function(...)
    rnorm(m * p, ...), control_rnorm)
  RM <- matrix(vals, nrow = m, ncol = p)
  return(RM)
}
#'
#' Creates an object class "\code{randomprojection}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_data_fun}, \code{update_rpm_w_data}
#' @return object of class "\code{randomprojection}" with is a list with elements name,
#' generate_fun, update_data_fun, control
#' @description
#' The entries of the matrix will be generated from
#' a normal distribution (mean 0 and standard deviation 1 by default).
#'
#' @export
#'
rp_gaussian <- constructor_randomprojection(
  "rp_gaussian",
  generate_fun = generate_gaussian
)

#'
#' Sparse random projection matrix
#'
#' @param rp object of class  "\code{randomprojection}"
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_sparse <- function(rp, m, included_vector) {
  p <- length(included_vector)
  psi <- attr(rp, "psi")
  if (is.null(psi)) psi <- 1
  if (psi > 1 | psi <= 0) stop("For a sparse rpm, psi should lie in interval (0,1].")
  v <- sample(c(-1, 0, 1), size = m * p,
              prob = c(psi/2, 1 - psi, psi/2), replace=TRUE)
  RM <- matrix(v/sqrt(psi), nrow = m, ncol = p)
  RM <- RM[rowSums(abs(RM)) > 0, ]
  RM <- Matrix(RM, sparse = TRUE)
  return(RM)
}

#'
#' Sparse random projection matrix
#'
#' Creates an object class "\code{randomprojection}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix. The possible argument is \code{psi} in (0,1] which determines
#' the level of sparsity in the matrix.
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_data_fun}, \code{update_rpm_w_data}
#' @return object of class "\code{randomprojection}"
#' @description
#' The sparse matrix used in \insertCite{ACHLIOPTAS2003JL}{spar} with entries equal to
#' \eqn{\Psi_{ij} = \pm 1/\sqrt{\psi}} with probability \eqn{\psi/2} and zero otherwise
#' for \eqn{\psi\in (0,1]}. Default is \code{psi = 1}.
#' @references{
#'   \insertRef{ACHLIOPTAS2003JL}{spar}
#' }
#' @export
rp_sparse <- constructor_randomprojection(
  "rp_sparse",
  generate_fun = generate_sparse
)


#'
#' Sparse embedding matrix
#'
#' @param m goal dimension, which will be randomly sampled in the SPAR algorithm
#' @param included_vector integer vector of column indices for the variables to be
#' included in the random projection. These indices are produced in the
#' screening step of the SPAR algorithm.
#' @return (possibly sparse) matrix with m rows and
#'  \code{length(included_vector)} columns.
#' @keywords internal
generate_cw <- function(rp, m, included_vector) {
  p <- length(included_vector)
  use_data <- attr(rp, "data")
  if (!use_data) {
    diagvals <- sample(c(-1, 1), p, replace = TRUE)
  } else {
    if (is.null(attr(rp, "diagvals")))
      stop("Must provide vector of coefficients for data-driven RP.")
    diagvals <- attr(rp, "diagvals")[included_vector]
  }
  goal_dims <- sample(m, p, replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in seq_len(m)) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims > goal_dim - counter] <-
        goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }
  RM <- Matrix(0, nrow = m - counter, ncol = p,sparse = TRUE)
  RM@i <- as.integer(goal_dims - 1)
  RM@p <- 0:p
  RM@x <- diagvals
  return(RM)
}

update_data_cw <- function(rp, x, y, screencoef) {
  n <- NROW(x)
  p <- NCOL(x)
  # dots <- list(...)
  if (attr(screencoef, "reuse_in_rp")) {
    scr_coef <- attr(screencoef, "importance")
    inc_probs <- attr(screencoef, "inc_prob")
  } else {
    if (is.null(rp$control$family)) {
      rp$control$family <- attr(rp, "family")
    }
    family <- rp$control$family

    if (family$family=="gaussian") {
      dev.ratio_cutoff <- 0.999
    } else {
      dev.ratio_cutoff <- 0.8
    }

    if (is.null(rp$control$alpha)) rp$control$alpha <-  0
    if (is.null(rp$control$lambda.min.ratio)) {
      tmp_sc <- apply(x, 2, function(col) sqrt(var(col)*(n-1)/n))
      x2 <- scale(x, center = colMeans(x), scale = tmp_sc)
      ytX <- crossprod(y, x2[,tmp_sc > 0])
      lam_max <- 1000 * max(abs(ytX))/n *
        family$mu.eta(family$linkfun(mean(y)))/
        family$variance(mean(y))
      rp$control$lambda.min.ratio <- min(0.01, 1e-4 / lam_max)
    }

    control_glmnet <- rp$control[names(rp$control)  %in% names(formals(glmnet))]
    glmnet_res <- do.call(function(...)
      glmnet(x=x, y=y, ...), control_glmnet)

    lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
    scr_coef <- coef(glmnet_res,s=lam)[-1]
    inc_probs <- abs(scr_coef)
    max_inc_probs <- max(inc_probs)
    attr(rp, "diagvals") <- scr_coef/max_inc_probs
  }
  attr(rp, "diagvals") <- scr_coef/max(inc_probs)
  return(rp)
}

update_rpm_w_data_cw <- function(rpm, rp, included_vector) {
  rpm@x <-  attr(rp, "diagvals")[included_vector]
  return(rpm)
}
#'
#' Sparse embedding matrix
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @param control list of arguments to be used in functions
#' \code{generate_fun}, \code{update_data_fun}, \code{update_rpm_w_data}
#' @return object of class randomprojection
#' @description
#' The entries of the matrix are generated based on \insertCite{Clarkson2013LowRankApprox}{spar}.
#' @references{
#'   \insertRef{Clarkson2013LowRankApprox}{spar}
#' }
#' @export
rp_cw <- constructor_randomprojection(
  "rp_cw",
  generate_fun = generate_cw,
  update_data_fun = update_data_cw,
  update_rpm_w_data = update_rpm_w_data_cw
)

#' print.randomprojection
#'
#' Print method for a "\code{randomprojection}" object
#' @param x object of class randomprojection
#' @param ... further arguments passed to or from other methods
#' @return text summary
#'
#' @export
print.randomprojection <- function(x, ...) {
  cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  cat("* Data-dependent:", attr(x,"data"), "\n")
  cat("* Lower bound on goal dimension m:",
      ifelse(is.null(attr(x, "mslow")),
             "not provided, will default to log(p).",
             attr(x, "mslow")), "\n")
  cat("* Upper bound on goal dimension m:",
      ifelse(is.null(attr(x, "msup")),
             "not provided, will default to n/2.",
             attr(x, "mslow")), "\n")
}
