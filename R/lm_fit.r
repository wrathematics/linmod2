#' Linear Model Fitter
#' 
#' @description
#' Basic linear model fitters.
#' 
#' @details
#' \code{lm_fit()} corresponds to R's \code{lm.fit()}.  The returned object is
#' mostly identical between the two functions.  There are a few key differences
#' worth noting.  First, the QR algorithm employed is different, and so the
#' \code{qr} matrix and \code{qraux} vector are unlikely to agree.
#' Additionally, when the data matrix \code{x} is assumed to be full rank and
#' it has more columns (predictors) than rows (observations), the returned
#' coefficients vector will not have \code{NA}'s after the first "rank" values,
#' as R's \code{lm.fit()} does.  Finally, we always return a matrix of
#' coefficients, effects, fitted values, and residuals, regardless of the number
#' of "right hand sides" (number of columns of the response \code{y}).  +R's
#' \code{lm.fit()} will return a vector when there is one right hand side, and
#' a matrix otherwise.
#' 
#' \code{.lm_fit()} corresponds to R's \code{.lm.fit()}.  The above caveats
#' also hold.
#' 
#' \code{.lm_fit_minimal()} does not correspond to an existing R function.  It
#' exists for personal reasons.  Note that it does not do less work than
#' \code{.lm_fit()}. It is "minimal" in the return, in that it does not return
#' unnecessary QR or miscellaneous, easily inferred information about the model
#' fit.  The return is subject to change at this time.
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
#' @return
#' It's complicated. The returns are mostly identical to their corresponding R
#' function counterparts.  However, there are a few notable differences.  See
#' the details section for more information..
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

#' @rdname lm_fit
#' @export
lm_fit <- function(x, y, intercept=FALSE)
{
  .Call(R_lm_fit, x, y, intercept)
}
