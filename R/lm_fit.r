#' @title Linear Model Fitter
#' 
#' A basic linear model fitter.
#' 
#' @description
#' The "minimal" function does not do less work. Instead, it is minimal in the
#' return, in that it does not return unnecessary QR or miscellaneous, easily
#' inferred information.
#' 
#' @param x 
#' The input data matrix.
#' @param y 
#' The vector/matrix of independent variable(s).
#' @param intercept 
#' Should an intercept term be included?
#' 
#' @examples
#' \dontrun{
#' library(linmod2)
#' 
#' n <- 10
#' p <- 3
#' x <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' 
#' .lm_fit(x, y)
#' }
#' 
#' @rdname lm_fit
#' @name lm_fit
NULL

#' @rdname lm_fit
#' @export
.lm_fit_minimal <- function(x, y, intercept=FALSE)
{
  .Call(R_dot_lm_fit_minimal, x, y, intercept)
}

#' @rdname lm_fit
#' @export
.lm_fit <- function(x, y, intercept=FALSE)
{
  .Call(R_dot_lm_fit, x, y, intercept)
}
