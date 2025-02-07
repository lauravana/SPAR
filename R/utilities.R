check_and_set_args <- function(args, x, y, family, model, screencoef, rp,  measure) {
  ## Check screen coef
  if (is.null(screencoef)) {
    screencoef <- screen_glmnet(nscreen = ncol(x))
  }
  if (!is.null(args$nscreen)) attr(screencoef, "nscreen") <- args$nscreen
  ##  Check if the old argument name 'old_arg' is used
  if (!is.null(args$type.measure)) {
    if (!is.null(measure)) {
      warning("Both 'measure' and deprecated 'type.measure' were provided. Using 'measure'.")
    } else {
      # Assign the value from 'old_arg' to 'new_arg' if 'new_arg' is not provided
      measure <- args$type.measure
      warning("'type.measure' is deprecated. Please use 'measure' instead.")
    }
  }
  if (!is.null(args$type.rpm)) {
    if (!is.null(rp)) {
      warning("Both 'rp' and deprecated 'type.rpm' were provided. Using 'rp'.")
    } else {
      # Assign the value from 'old_arg' to 'new_arg' if 'new_arg' is not provided
      rp <- switch(args$type.rpm,
                   "cw"           = rp_cw(data = FALSE),
                   "cwdatadriven" = rp_cw(data = TRUE),
                   "gaussian"     = rp_gaussian(),
                   "sparse"       = rp_sparse(psi = 0.1),
                   stop("Provided 'type.rpm' not implemented."))
      warning("'type.rpm' is deprecated. Please use 'rp' instead.")
    }
  }
  if (!is.null(args$type.screening)) {
    if (!is.null(screencoef)) {
      warning("Both 'screencoef' and deprecated 'type.screening' were provided. Using 'screencoef'.")
    } else {
      # Assign the value from 'old_arg' to 'new_arg' if 'new_arg' is not provided
      screencoef <- switch(args$type.screening,
                           "ridge" = screen_glmnet(),
                           "marglik" = screen_marglik(),
                           "corr" = screen_cor(),
                           stop("Provided 'type.screening' not implemented."))
      warning("'type.screening' is deprecated. Please use 'screencoef' instead.")
    }
  }
  ## TODO
  if (is.null(rp)) rp <- rp_cw(data = TRUE)
  if (!is.null(args$mslow)) attr(rp, "mslow") <- args$mslow
  if (!is.null(args$msup))  attr(rp, "msup") <- args$msup

  if (is.null(model)) {
    if (family$family == "gaussian" && family$link == "identity") {
      model <- spar_glm()
    } else {
      model <- spar_glmnet()
    }
  }
  out <- list(model = model, rp = rp,
              screencoef = screencoef, measure = measure)
  return(out)
}


get_val_measure_function <- function(measure, family) {
  val.meas <- switch(measure,
                     "deviance" = function(yval, eta_hat) {
                       sum(family$dev.resids(yval, family$linkinv(eta_hat), 1))
                     },
                     "mse" = function(yval, eta_hat) {
                       mean((yval - family$linkinv(eta_hat))^2)
                     },
                     "mae" = function(yval, eta_hat) {
                       mean(abs(yval - family$linkinv(eta_hat)))
                     },
                     "class" = {
                       stopifnot(family$family == "binomial")
                       function(yval, eta_hat) {
                         mean(yval != round(family$linkinv(eta_hat)))
                       }
                     },
                     "1-auc" = {
                       stopifnot(family$family == "binomial")
                       function(yval, eta_hat) {
                         if (var(yval) == 0) {
                           NA
                         } else {
                           phat <- prediction(family$linkinv(eta_hat), yval)
                           1 - performance(phat, measure = "auc")@y.values[[1]]
                         }
                       }
                     },
                     stop("Invalid measure")
  )
  return(val.meas)
}
