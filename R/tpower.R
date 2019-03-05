#' Compute a truncated power basis
#' @author Paul Eilers
#' @param x a vector on which the basis is calculated.
#' @param knots a vector of truncation points (knots) for piecewise polynomials.
#' @param p a scalar power for the basis, e.g. p=3 for cubic TPF.
#' @return a TPF basis of degree \code{p}.
#' @examples
#' library(JOPS)
#' # Basis  on grid
#' xg = seq(0, 4, length = 500)
#' Tg = cbind(tpower(xg,0, p=1), tpower(xg,1, p=1), tpower(xg,2,1))
#' matplot(xg, Tg, type='l',lty=c(1:6), lwd=2, xlab='x', ylab='',
#' main='Linear TPF basis (knots: 1, 2, 3)')
#'
#' @export

tpower <- function(x, knots, p) {
  (x - knots) ^ p * (x > knots)
}

