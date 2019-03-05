#' Smoothing scattered binomial data using P-splines.
#'
#' @description psBinomial is used to smooth scattered
#' binomial data using P-splines usign a logit link function.
#'
#' @param y the response vector, usually 0/1 data.
#' @param x the vector for the continuous regressor of \code{length(y)} and
#' the abcissae of fit.
#' @param xl the number for the min along \code{x}.
#' @param xr the number for the max along \code{x}.
#' @param nseg the number of evenly spaced segments between xl and xr.
#' @param bdeg the number of the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the number of the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) number for the tuning parameter for the penalty.
#' @param wts the vector of weights, default is 1, 0/1 allowed.
#' @param show Set to T or F to display iteration history.
#'
#' @return
#' \item{pcoeff}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{muhat}{a vector of length \code{m} of estimated means (probabilities).}
#' \item{dev}{deviance}
#' \item{effdim}{effective dimension of the smooth}
#' \item{aic}{AIC}
#' \item{wts}{a vector of preset weights.}
#' \item{nseg, bdeg, pord, lambda}{design parameters.}
#' \item{xgrid}{gridded x values for plotting.}
#' \item{ygrid}{gridded fitted probability values for plotting.}
#' \item{lgrid}{gridded lower 2se values for plotting.}
#' \item{ugrid}{gridded lower 2se values for plotting.}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.

#' @examples
#' library(JOPS)
#' #Extract data
#' library(rpart)
#' Kyphosis=kyphosis$Kyphosis
#' Age=kyphosis$Age
#' y = 1 * (Kyphosis == "present")  # make y 0/1
#' fit1 = psBinomial(Age, y,xl=min(Age), xr=max(Age), nseg=20,
#'                  bdeg=3, pord=2, lambda= 1)
#' names(fit1)
#' plot(fit1, xlab="Age", ylab='0/1', se=2)
#' @export

psBinomial = function(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
                      pord = 2, lambda = 1, wts = NULL, show = F) {

    # Compute B-spline basis
  m = length(x)
  B = bbase(x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg)
  # Construct penalty stuff
  n = dim(B)[2]
  P = sqrt(lambda) * diff(diag(n), diff = pord)
  nix = rep(0, n - pord)

  # Initialize
  znew = 0.5  #log(y + 0.01)
  z = rep(1, m)
  if(missing(wts)){wts=rep(1,m)}

  # Fit
  for (it in 1:50) {
    mu = exp(z)/(1 + exp(z))
    w = mu * (1 - mu)
    u = (y - mu)/w + z
    wtprod=c(wts*w,(nix+1))
    f = lsfit(rbind(B, P), c(u, nix), wt = wtprod, intercept = F)
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
  e = 1e-09
  dev = 2 * sum((y + e) * log((y + e)/mu) +
                  (1 - y + e) * log((1 - y + e)/(1 - mu)))
  h = hat(f$qr)[1:m]
  ed = sum(h)
  aic = dev + 2 * ed

  # Compute curve on grid
  xgrid = seq(xl, xr, length = 500)
  Bu = bbase(xgrid, xl=xl, xr=xr, nseg = nseg, bdeg = bdeg)
  zu = Bu %*% beta
  ygrid = exp(zu)/(1+exp(zu))

  # SE bands on a grid
  w = exp(mu)/((1+exp(mu))^2)
  wtprod=as.vector(w)*as.vector(wts)
  Bayese = solve(t(wtprod*B) %*%  B +  t(P) %*% P)
  vareta = Bu %*% Bayese %*% t(Bu)
  seeta = 2*sqrt(diag(vareta))
  ugrid = exp(zu + seeta)/(1+exp(zu+seeta))
  lgrid = exp(zu - seeta)/(1+exp(zu-seeta))

  # Return list
  pp = list(aic = aic, B=B, P=P, wts=wts, xl=xl, xr=xr, x = x,
            y = y, muhat = mu, nseg = nseg, bdeg = bdeg,
            pord = pord, pcoef = beta, lambda = lambda,
            effdim = ed, dispersion = 1,
            family = "binomial",  link = "logit",
            xgrid=xgrid, ygrid=ygrid, ugrid=ugrid, lgrid=lgrid)
 class(pp)=c("pspfit")
 return(pp)
}
