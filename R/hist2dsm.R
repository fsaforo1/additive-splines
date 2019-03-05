#' Smoothing 2D histograms using P-splines
#'
#' @description Fits 2D smooth P-spline surface using
#' Poisson response layed out in a 2D histogram
#'
#' @param Y a matrix of 2D histogram counts.
#' @param nsegx the number of evenly spaced knots along x for Tensor product B-spline basis.
#' @param nsegy the number of evenly spaced knots along y for Tensor product B-spline basis.
#' @param ndeg the integer (e.g. 1, 2, 3) for the degree of the B-spline tensor basis.
#' @param lambdax the positive number for the tuning parameter along \code{x}.
#' @param lambday the positive number for the tuning parameter along \code{y}.
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the degree of the basis, usually 1, 2, or 3.
#' @param dx the integer (e.g. 1, 2, 3) for the order of the penalty along \code{x}.
#' @param dy the integer (e.g. 1, 2, 3) for the order of the penalty along \code{y}.
#' @param kappa a (small, positive) number for ridge tuning parameter, if needed.
#'
#' @return A list comprising of
#' \item{ED}{the effective dimension of the smooth 2D surface}
#' \item{Mu}{a matrix with \code{dim(Y)} consisting of
#' P-spline density estimates}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' S=hist2dsm(hist2d(runif(1000),runif(1000),c(100,100))$H, nsegx = 5,
#' nsegy = 5, bdeg = 3, lambdax = 10, lambday = 10, dx = 3, dy =3,
#' Mu.start = NULL, kappa = 1e-4)
#' image(seq(0,1,length=100),seq(0,1,length=100),S$Mu, xlab='x', ylab='y')
#'
#' @export

hist2dsm = function(Y, nsegx = 10, nsegy = nsegx, bdeg = 3,
                    lambdax = 10, lambday = lambdax, dx = 3, dy = dx,
                    Mu.start = NULL, kappa = 1e-4) {
  # 2-D P-spline histogram estimation
    nx = nrow(Y)
    ny = ncol(Y)
    # cat(nx, ny, '\n')

    Bx = bbase(1:nx, 0, nx + 1, nsegx, bdeg)
    By = bbase(1:ny, 0, ny + 1, nsegy, bdeg)
    nbx = ncol(Bx)
    nby = ncol(By)
    Tx = rowtens(Bx)
    Ty = rowtens(By)

    Dx = diff(diag(nbx), diff = dx)
    Dy = diff(diag(nby), diff = dy)
    Px = lambdax * t(Dx) %*% Dx
    Py = lambday * t(Dy) %*% Dy
    P = kronecker(Py, diag(nbx)) + kronecker(diag(nby), Px)
    P = P + kappa * diag(nrow(P))

    # Initialize
    Mu = Y + 0.01
    if (!is.null(Mu.start)) Mu = Mu.start
    Z = log(Mu)
    Z = Z - log(sum(exp(Z)) / sum(Y))

    # Iterate
    for (it in 1:30) {
      Mu = exp(Z)
      U = Y - Mu + Mu * Z
      Q = t(Tx) %*% Mu %*% Ty
      dim(Q) = c(nbx, nbx, nby, nby)
      Q = aperm(Q, c(1, 3, 2, 4))
      dim(Q) = c(nbx * nby, nbx * nby)
      r = t(Bx) %*% U %*% By
      dim(r) = c(nbx * nby, 1)
      A = solve(Q + P, r)
      a = A
      dim(A) = c(nbx, nby)
      Znew = Bx %*% A %*% t(By)
      dz = sum(abs(Z - Znew))
      if (dz < 1e-5) break
      Z = Znew
    }

    K = solve(Q + P, Q)
    ed = sum(diag(K))

    apen = t(a) %*% P %*% a
    return(list(ed = ed, Mu = Mu, Y = Y, apen = apen))
}


