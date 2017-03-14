#' Linear Model Fitter
#' 
#' A basic linear model fitter.
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
#' library(linmod)
#' 
#' n <- 10
#' p <- 3
#' x <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' 
#' .lm_fit(x, y)
#' }
#' 
#' @rdname dotlm_fit
#' @export
.lm_fit <- function(x, y, intercept=FALSE)
{
  .Call(R_lm_fit, x, y, intercept)
}
