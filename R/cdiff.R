#' Constructs a second order cyclical differencing matrix
#' @description A difference matrix used for cyclical penalities, where
#' the ends are penalized to connnect and extrapolation is performed on a cosine curve.
#'
#' @param n scalar that determines the dimension of the square differencing matrix.
#' @return a square difference penalty matrix with
#' \code{n} rows and \code{n} columns.
#' @author Paul Eilers
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' cdiff(5)
#'
#' @export

cdiff = function(n) {
  # Compute cyclic difference matrix
  D2 = matrix(0, n, n + 2)
  p = c(-1, 2 * cos(2 * pi/n), -1)
  for (k in 1:n) D2[k, (0:2) + k] = p
  D = D2[, 2:(n + 1)]
  D[, 1] = D[, 1] + D2[, n + 2]
  D[, n] = D[, n] + D2[, 1]
  D
}
