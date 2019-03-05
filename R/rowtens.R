#' Compute row tensor product of a matrix
#' @description Compute row tensor product of a matrix a
#' support function for 2-D histogram smoothing.
#' @author Paul Eilers
#'
#' @export

rowtens = function(X, Y=X){
  # Row-wise tensor products
  one = matrix(1, nrow = 1, ncol = ncol(X))
  kronecker(X, one) * kronecker(one, Y)
}
