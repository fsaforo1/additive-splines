#' B-spline basis on evenly spaced knots
#' @description Constructs a B-spline basis on evenly spaced knots, used
#' for P-splines.
#'
#' @param x a vector of argument values, at which the B-spline basis functions
#' are to be evaluated.
#' @param xl the lower limit of the domain of x.
#' @param xr the upper limit of the domain of x.
#' @param nseg the number of evenly spaced segments between xl and xr.
#' @param bdeg the degree of the basis, usually 1, 2, or 3.
#' @return a matrix with number of rows=\code{length(x)} and
#' number of columns=\code{nseg+bdeg}.
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
#' B-splines and penalties (with comments and rejoinder), Statistical Science,
#' 11: 89-121.
#'
#' @examples
#' library(JOPS)
# # Basis  on grid
#' xg = seq(0, 4, length = 500)
#' Bg = bbase(xg, 0, 4, nseg =3, bdeg = 3)
#' nb1 = ncol(Bg)
#' matplot(xg, (Bg), type='l',lty=c(1:6), lwd=2, xlab='x', ylab='')
#' @export

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {
  # Function for B-spline basis Evalated at x Left and right boundaries
  # in xl and xr Number of segments (between xl and xr) nseg Degree of
  # the splines bdeg (3 gives cubic splines)

  dx <- (xr - xl)/nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  B
}

