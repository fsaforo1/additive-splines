#' Smoothing scattered Poisson data using P-splines.
#'
#' @description psPoisson is used to smooth scattered
#' Poisson data using P-splines usign a log link function.
#'
#' @param y the response vector, usually count data.
#' @param x the vector for the continuous regressor of \code{length(y)} and
#' the abcissae of fit.
#' @param xl the number for the min along \code{x}.
#' @param xr the number for the max along \code{x}.
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the number of the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the number of the order of the difference penalty, usually 1, 2, or 3. (defalult)
#' @param lambda the (positive) number for the tuning parameter for the penalty.
#' @param wts the vector of weights, default is 1, 0/1 allowed
#' @param show Set to T or F to display iteration history.
#'
#' @return
#' \item{pcoeff}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{muhat}{a vector of length \code{m} of estimated means.}
#' \item{dev}{Deviance of fit.}
#' \item{effdim}{effective dimension of fit.}
#' \item{aic}{AIC}
#' \item{wts}{a vector of weights.}
#' \item{nseg, bdeg, pord, lambda}{design parameters.}
#' \item{xgrid}{gridded x values for plotting.}
#' \item{ygrid}{gridded fitted mean values for plotting.}
#' \item{lgrid}{gridded lower 2se values for plotting.}
#' \item{ugrid}{gridded lower 2se values for plotting.}
#'
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' library(JOPS)
#' library(boot)

#' # Extract the data
#' Count=hist(coal$date, breaks=c(1851:1963),plot=F)$counts
#' Year=c(1851:1962)
#' xl=min(Year)
#' xr=max(Year)

#' # Poisson smoothing
#' nseg=20
#' bdeg=3
#' fit1=psPoisson(Year, Count, xl,xr,nseg,bdeg,pord=2,lambda=1)
#' names(fit1)
#' plot(fit1, xlab="Year", ylab="Count", se=2)
#' @export

psPoisson = function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                     pord = 2, lambda = 1, wts = NULL, show = F) {

  # Compute B-spline basis
  m = length(x)
  B = bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)

  # Construct penalty stuff
  n = dim(B)[2]
  P = sqrt(lambda) * diff(diag(n), diff = pord)
  nix = rep(0, n - pord)

  # Initialize
  z = log(y + 0.01)
  if(missing(wts)){wts=rep(1,m)}
  # Fit
  for (it in 1:50) {
    mu = exp(z)
    w = mu
    u = (y - mu) / w + z
    wtprod=c(wts*w,(nix+1))
    f = lsfit(rbind(B, P), c(u, nix), intercept=F, wt = wtprod)
    beta = f$coef
    znew = B %*% beta
    dz = max(abs(z - znew))
    z = znew
    if (dz < 1e-06)
      break
    if (show)
      print(c(it, dz))
  }

  if(it > 49) {
    warning(paste("Did NOT converge in 50 iterations"))}

  # Compute AIC
  dev = 2 * sum(y * log((y + 1e-09)/mu))
  h = hat(f$qr)[1:m]
  ed = sum(h)
  aic = dev + 2 * ed

  # Compute curve on grid
  xgrid = seq(xl, xr, length = 500)
  Bu = bbase(xgrid, xl=xl, xr=xr, nseg = nseg, bdeg = bdeg)
  zu = Bu %*% beta
  ygrid = exp(zu)

  # SE bands on a grid
    w = exp(mu)
    wtprod=as.vector(w)*as.vector(wts)
    Bayese = solve(t(wtprod*B) %*%  B +  t(P) %*% P)
    vareta = Bu %*% Bayese %*% t(Bu)
    seeta = 2*sqrt(diag(vareta))
    ugrid = exp(zu + seeta)
    lgrid = exp(zu - seeta)

  # Return list
  pp = list(xl=xl,xr=xr,aic = aic, x = x, y = y, B=B, P=P, muhat = mu, nseg = nseg, bdeg = bdeg,
            pord = pord, pcoef = beta, lambda = lambda,
            effdim = ed, dispersion = 1,
            family = "poisson", link = "log", wts = wts,
            xgrid=xgrid, ygrid=ygrid, ugrid=ugrid, lgrid=lgrid)
  class(pp) = "pspfit"
  return(pp)
}
