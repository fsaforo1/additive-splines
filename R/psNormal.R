#' Smoothing scattered (normal) data using P-splines
#'
#' @description psNormal is used to smooth scattered (normal) data using P-splines(with identity link function)
#'
#' @param y the response vector, usually continuous data.
#' @param x the vector for the continuous regressor of \code{length(y)} and the abcissae of fit.
#' @param wts the vector of weights, default is 1; 0/1 allowed.
#' @param xl the number for the min along \code{x}.
#' @param xr the number for the max along \code{x}.
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the number of the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the number of the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) number for the tuning parameter for the penalty.
#' @param show set to T or F to display iteration history.
#'
#' @return
#' \item{coeff}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{muhat}{a vector of length \code{m} of smooth estimated means.}
#' \item{wts}{a vector of length \code{m} of weights}
#' \item{dev}{deviance}
#' \item{effdim}{estimated effective dimension}
#' \item{aic}{AIC}
#' \item{ed.resid}{approximate residual df}
#' \item{sigma}{root MSE}
#' \item{cv}{a statistic providing the standard error of leave-one-out prediction.}
#' \item{nseg, bdeg, pord, lambda}{design parameters.}
#' \item{xgrid}{gridded x values for plotting.}
#' \item{ygrid}{gridded fitted mean values for plotting.}
#' \item{lgrid}{gridded lower 2se values for plotting.}
#' \item{ugrid}{gridded lower 2se values for plotting.}
#'
#' @author  Paul Eilers and Brian Marx
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' library(JOPS)
#' #Extract data
#' library(MASS)
#' #Get the data
#' data(mcycle)
#' x = mcycle$times
#' y = mcycle$accel
#' fit1 = psNormal(x, y, nseg = 20, bdeg = 3, pord = 2, lambda = .8)
#' names(fit1)
#' plot(fit1, se=2, xlab="time (ms)", ylab="accel")
#' @export
#'
psNormal = function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                    pord = 2, lambda = 1, wts=NULL) {
  m = length(x)
  B = bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)

  # Construct penalty stuff
  n = dim(B)[2]
  P = sqrt(lambda) * diff(diag(n), diff = pord)
  nix = rep(0, n - pord)


  # Fit
  if(missing(wts)){wts=rep(1,m)}
  f = lsfit(rbind(B, P), c(y, nix), intercept = F, wt=c(wts,(nix+1)))
  h = hat(f$qr)[1:m]
  beta = f$coef
  mu = B %*% beta

  # Cross-validation and dispersion
  r = (y - mu)/(1 - h)
  cv = sqrt(mean(r^2))
  ed = sum(h)
  sigma = sqrt(sum((y - mu)^2)/(m - ed))

  # Compute curve on grid
  xgrid = seq(xl, xr, length = 500)
  Bu = bbase(xgrid, xl=xl, xr=xr, nseg = nseg, bdeg = bdeg)
  zu = Bu %*% beta
  ygrid = zu

  # SE bands on a grid
  Bayese = solve(t(wts*B) %*%  B +  t(P) %*% P)
  vareta = sigma^2*(Bu %*% Bayese %*% t(Bu))
  seeta = 2*sqrt(diag(vareta))
  ugrid = zu + seeta
  lgrid = zu - seeta

  # Return list
  pp = list(x = x, y = y, B=B, P=P, muhat = mu, nseg = nseg, xl = xl,
            xr = xr, bdeg = bdeg, pord = pord, lambda = lambda,
            cv = cv, effdim = ed, ed.resid = m - ed, wts=wts,
            pcoeff = beta, family = "gaussian", link = "identity",
            sigma = sigma, xgrid=xgrid, ygrid=ygrid, ugrid=ugrid,
            lgrid=lgrid)
  class(pp) = "pspfit"
  return(pp)
}
