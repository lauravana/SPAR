#' Constructor function for building sparmodel objects
#'
#' Creates an object of class "\code{sparmodel}" using arguments passed by user.
#' @param name character
#' @param model_fun function for estimating the marginal models which returns the
#     intercept and the vector of coefficients. This
#'    function should have arguments  and   \code{y} (vector of responses -- standardized
#'    for Gaussian family), \code{z} (the matrix of projected predictors) and a
#'    "\code{sparmodel}" \code{object}.
#' @param update_fun optional function for updating the "\code{sparmodel}" object. before the
#' start of the algorithm.
#' @return a function which in turn creates a function which in turn creates an
#'    object of class "\code{sparmodel}".
#' @description
#' The created function will return a object of class "\code{sparmodel}" which
#' constitutes of a list.
#' @export
constructor_sparmodel <- function(name, model_fun, update_fun = NULL) {
  ## Checks
  args_generate_fun <- formals(model_fun)
  stopifnot("Function model_fun should contain three arguments: y, z and an object
            of class \"sparmodel\"." =
              length(args_generate_fun) == 3)
  stopifnot("Function model_fun should contain argument 'y', the vector of responses." =
              "y" %in% names(args_generate_fun))
  stopifnot("Function model_fun should contain argument 'z', the matrix of reduced predictors." =
              "z" %in% names(args_generate_fun))
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                model_fun = model_fun,
                update_fun = update_fun,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    class(out) <- c("sparmodel")
    return(out)
  }
}

#' Model object for estimating penalized glms in marginal models in the ensemble
#'
#' Creates an object class "\code{sparmodel}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{sparmodel}" object
#' @param control list of controls to be passed to the model function
#' @return object of class "\code{sparmodel}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{model_fun}  for generating the screening coefficient.
#'   This function should have arguments \code{y}, vector of standardized responses,
#'   \code{z}, a matrix of projected predictors in each marginal model, and \code{object}, which is a "\code{sparmodel}" object. Returns a list with two elements: \code{gammas} which is the vector of regression coefficients for the projected predictors and \code{intercept} which is the intercept of the model.
#'  \item \code{update_fun}  optional function for updating the sparmodel object. before the
#' start of the algorithm.
#' }
#' @description
#' Relies on \link[glmnet]{glmnet}.
#'
#' @export
#'
spar_glmnet <- function(..., control = list()) {
  ## Set defaults
  if (is.null(control$alpha)) {
    control$alpha <- 0
  }
  out <-  list(name = "glmnet",
               model_fun = model_glmnet,
               update_fun = update_sparmodel_glmnet,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

ols_fun <- function(y, z) {
  solve(crossprod(z), crossprod(z,y))
}
ols_fun_corrected <- function(y, z) {
  solve(crossprod(z) + 0.01*diag(ncol(z)), crossprod(z,y))
}

update_sparmodel_glmnet <- function(object) {
  family <- object$control$family
  if (family$family=="gaussian" & family$link=="identity") {
    fit_family <- "gaussian"
  } else {
    if (family$family=="binomial" & family$link=="logit") {
      fit_family <- "binomial"
    } else if (family$family=="poisson" & family$link=="log") {
      fit_family <- "poisson"
    } else {
      fit_family <- family
    }
  }
  attr(object, "family") <- fit_family
  object
}

model_glmnet <- function(y, z, object) {
  ## y - vector of n responses
  ## z - matrix with n rows
  glmnet_res <- do.call(function(...) glmnet(x = z, y = y, ...),
                        object$control)
  mar_coef <- coef(glmnet_res, s = min(glmnet_res$lambda))
  intercept <- mar_coef[1]
  gammas <- mar_coef[-1]
  list(gammas = gammas, intercept = intercept)
}

#' Model object for estimating penalized glms in marginal models in the ensemble models in the ensemble
#'
#' Creates an object class "\code{sparmodel}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{sparmodel}" object
#' @param control list of controls to be passed to the model function
#' @return object of class "\code{sparmodel}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{model_fun}  for generating the screening coefficient.
#'   This function should have arguments \code{y}, vector of standardized responses,
#'   \code{z}, a matrix of projected predictors in each marginal model, and \code{object}, which is a "\code{sparmodel}" object. Returns a list with two elements: \code{gammas} which is the vector of regression coefficients for the projected predictors and \code{intercept} which is the intercept of the model
#' }
#' @description
#' Relies on \link[stats]{glm.fit}.
#'
#' @export
#'
spar_glm <- function(..., control = list()) {
  out <-  list(name = "glm",
               model_fun = model_glm,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

model_glm <- function(y, z, object) {
  ## y - vector of n responses
  ## z - matrix with n rows
  family <- object$control$family
  if (family$family=="gaussian" & family$link=="identity") {
    gammas <- tryCatch(ols_fun(y, z),
                       error=function(error_message) {
                         return(ols_fun_corrected(y, z))
                       })
    intercept <- 0
  } else {
    glm_res <- do.call(function(...) glm(y ~ z, ...),
                       object$control)
    intercept <- coef(glm_res)[1]
    gammas <- coef(glm_res)[-1]
  }

  list(gammas = gammas, intercept = intercept)
}

