#' Circular B-spline basis
#' @description Constructs a circular B-spline basis on evenly spaced knots.
#'
#' @param x a vector of argument values, at which the B-spline basis functions
#' are to be evaluated.
#' @param xl the lower limit of the domain of x.
#' @param xr the upper limit of the domain of x.
#' @param nseg the number of evenly spaced segments between xl and xr.
#' @param bdeg the degree of the basis, usually 1, 2, or 3.
#' @return a matrix with number of rows=\code{length(x)} and
#' number of columns=\code{nseg}.
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M.
#' (2015). Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#'
#' @examples
#' library(JOPS)
# # Circular basis  on grid
#' xg = seq(0, 4, length = 500)
#' Bg = cbase(xg, 0, 4, nseg =3, bdeg = 3)
#' nb1 = ncol(Bg)
#' matplot(xg, (Bg), type='l',lty=c(1:6), lwd=2, xlab='x', ylab='')
#'
#' @export

cbase <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {
  # Construct circular B-spline basis Domain: xl to xr, number of
  # segments on domain: nseg, degree: bdeg
  # Wrap around to cyclic basis
  B0 = bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  n = ncol(B0) - bdeg
  cc = (1:bdeg) + n
  B = B0[, 1:n]
  B[, 1:bdeg] = B[, 1:bdeg] + B0[, cc]
  B
}

