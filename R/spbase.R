#' Sarse B-spline basis on evenly spaced knots
#' @description Constructs a B-spline basis on evenly spaced knots, used
#' for P-splines. A sparse matrix, in the format of the packag \code{spam}, is returned.
#'
#' @param x a vector of argument values, at which the B-spline basis functions
#' are to be evaluated.
#' @param xl the lower limit of the domain of x.
#' @param xr the upper limit of the domain of x.
#' @param nseg the number of evenly spaced segments between xl and xr.
#' @param bdeg the degree of the basis, usually 1, 2, or 3.
#' @return a sparse matrix with number of rows=\code{length(x)} and
#' number of columns = \code{nseg + bdeg}.
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
#' B-splines and penalties (with comments and rejoinder), Statistical Science,
#' 11: 89-121.
#'
#' @examples
#' library(JOPS)
#' # Basis  on grid
#' xg = seq(0, 4, length = 500)
#' Bg = spbase(xg, 0, 4, nseg = 3, bdeg = 3)
#' nb1 = ncol(Bg)
#' matplot(xg, Bg, type = 'l', lty = c(1:6), lwd = 2, xlab= 'x', ylab='')
#' @export

spbase = function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {
  # Compute a sparse B-spline matrix in spam format
  require(spam)

  # Reduce x to first interval between knots
  m = length(x)
  dx = (xr - xl) / nseg
  ix <- floor((x - xl) / (1.0000001 * dx))
  xr = (x - xl) - ix * dx

  # Full basis for reduced x
  Br = bbase(xr, xl = 0, xr = 0 + dx, nseg = 1, bdeg = bdeg)

  # Compute proper rows, columns
  nr = ncol(Br)
  rw = rep(1:m, each = nr)
  cl = rep(1:nr, m) + rep(ix, each = nr)

  # Make the sparse matrix
  b = as.vector(t(Br))
  Bs = spam(list(i = rw, j = cl, b), nrow = m, ncol = (nseg + bdeg))

  return(Bs)
}
