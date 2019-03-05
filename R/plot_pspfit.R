#' Plotting function for psNormal, psPoisson, psBinomial
#'
#' @description Plotting function for P-spline smooth (\code{class pspfit}), with se bands.
#'
#' @param ps the P-spline function, usually from psNormal, psPoisson, psBinomial.
#' @param se a scalar, e.g. \code{se}=2 to produce twice se bands, set \code{se} > 0 (or set \code{se}=0 to supress).
#' @param xlab label for the x-axis, e.g. "my x" (quotes required).
#' @param ylab label for the y-axis, e.g. "my y" (quotes required).
#' @param col color for points, e.g. "blue".
#' @param pch point character, e.g. pch=1 or pch="O".

#' @return
#' \item{Plot}{a plot of smoothed normal, Poisson, or binomial data, with or with se bands.}
#'
#' @author  Paul Eilers and Brian Marx
#' @references  Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' library(JOPS)
#' #Extract data
#' library(MASS)
#' # Get the data
#' data(mcycle)
#' x = mcycle$times
#' y = mcycle$accel
#' fit1 = psNormal(x, y, nseg = 20, bdeg = 3, pord = 2, lambda = .8)
#' plot(fit1, se=2, xlab="time (ms)", ylab="accel")
#' #'
#' @examples
#' library(JOPS)
#' library(boot)
#' # Extract the data
#' Count=hist(coal$date, breaks=c(1851:1963),plot=F)$counts
#' Year=c(1851:1962)
#' xl=min(Year)
#' xr=max(Year)
#'
#' # Poisson smoothing
#' nseg=20
#' bdeg=3
#' fit1=psPoisson(Year, Count, xl,xr,nseg,bdeg,pord=2,lambda=1)
#' names(fit1)
#' plot(fit1, xlab="Year", ylab="Count", se=2)
#'
#' @examples
#' library(JOPS)
#' #Extract data
#' library(rpart)
#' Kyphosis=kyphosis$Kyphosis
#' Age=kyphosis$Age
#' y = 1 * (Kyphosis == "present")  # make y 0/1
#' # Binomial smoothing
#' fit1 = psBinomial(Age, y,xl=min(Age), xr=max(Age), nseg=20,
#'                  bdeg=3, pord=2, lambda= 1)
#' names(fit1)
#' plot(fit1, xlab="Age", ylab='0/1', se=2)
#' @export
plot.pspfit = function(ps, se = 2, xlab="", ylab="", col='black', pch=1){
# Code block for psNormal
  if (ps$family == 'gaussian') {
    plot(ps$x, ps$y, main = "", xlab = xlab,
         ylab = ylab, pch=pch, col=col)
    # Compute curve on grid
    xl = ps$xl
    xr = ps$xr
    u = seq(ps$xl, ps$xr, length = 500)
    Bu = bbase(u, xl = ps$xl,xr = ps$xr,
               nseg = ps$nseg, bdeg = ps$bdeg)
    zu = Bu %*% ps$pcoeff
    lines(u, zu, col = 'blue')
    sefit = 0

    # Error bands (Bayesian estimate)
    ugrid2se = lgrid2se = NULL
    if (se > 0) {
      Covb = solve(t(ps$wts*ps$B) %*% ps$B + t(ps$P) %*% ps$P)
      Covz = ps$sigma ^ 2 * Bu %*% Covb %*% t(Bu)
      sefit = sqrt(diag(Covz))
      lgrid2se = zu - se * sefit
      ugrid2se = zu + se * sefit
      seb = se * sefit
      lines(u, zu + seb, lty = 2, col = 'red')
      lines(u, zu - seb, lty = 2, col = 'red')
    }
  }
# Code block for psPoisson
  if(ps$family=='poisson'){
    # Compute curve on grid
    xgrid = seq(ps$xl, ps$xr, length = 500)
    Bu = bbase(xgrid, xl = ps$xl, xr = ps$xr, nseg = ps$nseg, bdeg = ps$bdeg)
    zu = Bu %*% ps$pcoef
    ygrid = exp(zu)
    # Plot data and fit
    plot(ps$x, ps$y, xlab = xlab, ylab = ylab, col=col, pch=pch)
    lines(xgrid, exp(zu), col = "blue")

    # SE bands on a grid
    if(se>0){
      w = exp(ps$muhat)
      wtprod=as.vector(w)*as.vector(ps$wts)
      Bayese = solve(t(wtprod*ps$B) %*%  ps$B +  t(ps$P) %*% ps$P)
      vareta = Bu %*% Bayese %*% t(Bu)
      seeta = se*sqrt(diag(vareta))
      ugrid = exp(zu + seeta)
      lgrid = exp(zu - seeta)
      lines(xgrid,ugrid, col='red')
      lines(xgrid,lgrid, col='red')
    }
  }
# Code block for psBinomial
  if(ps$family=="binomial"){
    # Compute curve on grid
    xgrid = seq(ps$xl, ps$xr, length = 500)
    Bu = bbase(xgrid, xl = ps$xl, xr = ps$xr, nseg = ps$nseg, bdeg = ps$bdeg)
    zu = Bu %*% ps$pcoef
    ygrid = exp(zu)/(1 + exp(zu))
  # Plot data and fit
    plot(ps$x, ps$y, xlab = xlab, ylab = ylab,col=col,pch=pch)
    lines(xgrid, exp(zu)/(1 + exp(zu)), col = "blue")
  # SE bands on a grid
    if(se>0){
      w = exp(ps$muhat)/(1-exp(ps$muhat))^2
      wtprod=as.vector(w)*as.vector(ps$wts)
      Bayese = solve(t(wtprod*ps$B) %*%ps$B + t(ps$P) %*% ps$P)
      vareta = Bu %*% Bayese %*% t(Bu)
      seeta = se*sqrt(diag(vareta))
      ugrid = exp(zu + seeta)/(1+exp(zu +  seeta) )
      lgrid = exp(zu - seeta)/(1+exp(zu - seeta) )
      lines(xgrid,lgrid,col='red')
      lines(xgrid,ugrid,col='red')
    }
  }
    if (se < 0) {
    warning(paste("se should be nonnegative"))
  }
}
