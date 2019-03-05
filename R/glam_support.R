# Support functions for GLAM

rowtens = function(X) {
    # Row-wise tensor products
    one = matrix(1, nrow = 1, ncol = ncol(X))
    kronecker(X, one) * kronecker(one, X)
}
