#' Two-dimensional smoothing scattered (normal) data using P-splines.
#'
#' @description psNormal is used to smooth scattered
#' (normal) data using P-splines, with aniosotripic penalization of
#' tensor product B-splines.
#'
#' @param Data a matrix of 3 columns \code{x, y, z} of equal length;
#' the response is \code{z}.
#' @param Pars a matrix of 2 rows, where each the first (second) row
#' sets the P-spline paramters for \code{x (y)}: \code{min max nseg bdeg lambda pdeg}.
#' @param min(max) A scalar for the min(max) along x or y.
#' @param nseg the number of evenly spaced segments between min and max.
#' @param bdeg a number for the degree of the basis, usually 1, 2, or 3.
#' @param pord the number for the order of the difference penalty, usually 1, 2, or 3.
#' @param lambda the (positive) tuning parameter for the penalty.
#' @param XYpred a matrix with two columns \code{(x,y)} that give the coordinates
#' of (future) prediction; the default is the data locations.

#' @return
#' \item{coef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[1,3]+Pars[1,4])}
#' of (unfolded) estimated P-spline coefficients.}
#' \item{fit}{a vector of \code{length(y)} of smooth estimated means (at the \code{x,y} locations).}
#' \item{pred}{a vector of length \code{nrow(XYpred)} of (future) predictions.}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
#' Twenty years of P-splines, SORT, 39(2): 149-186.
#'
#' @examples
#' library(SemiPar)
#' library(fields)
#' library(spam)
#' library(JOPS)

#' # Get the data
#' data(ethanol)
#' m = nrow(ethanol)
#' x = ethanol$C
#' y = ethanol$E
#' z = ethanol$NOx

#' # Set parameters for domain
#' xlo <- 7
#' xhi <- 19
#' ylo <- 0.5
#' yhi <- 1.25

#' # Set P-spline parameters, fit and compute surface
#' xseg <- 10
#' xdeg <- 3
#' xpars <- c(xlo, xhi, xseg, xdeg, 3, 1)
#' yseg <- 10
#' ydeg <- 3
#' ypars <- c(ylo, yhi, yseg, ydeg, 3, 1)
#' Pars1=rbind(xpars,ypars)
#' resol=100
#' Xgrid=seq(xlo,xhi,length=resol)
#' Ygrid=seq(ylo,yhi,length=resol)
#' XYgrid=as.matrix(expand.grid(seq(xlo,xhi,length=resol),seq(ylo,yhi,length=resol)))
#' fit=psp2d(cbind(x,y,z),Pars1, XYpred=XYgrid)
#' image(Xgrid,Ygrid,matrix(fit$pred,resol,resol), xlab='C',ylab='E')
#' @export

psp2d <- function(Data, Pars, XYpred) {
  # Fitting with 2-D P-splines
  # Input:
  #   Data: 3 columns, giving x, y, z
  #   Pars: 2 rows with P-spline parameters: [min max nseg deg lambda pdeg]
  #   XYpred: 2 columns, giving x and y at points to predict
  # Output: list with elements
  #   coef: coeffcients
  #   fit: fitted values
  #   pred: predicted z
  #
  # Paul Eilers, 2000

  # Prepare bases
  p1 = Pars[1, ]
  p2 = Pars[2, ]
  Bx <- bbase(Data[, 1], p1[1], p1[2], p1[3], p1[4])
  By <- bbase(Data[, 2], p2[1], p2[2], p2[3], p2[4])
  nx <- ncol(Bx);
  ny <- ncol(By);

  # Compute tensor products
  B1 <- kronecker(Bx, t(rep(1, ny)));
  B2 <- kronecker(t(rep(1, nx)), By);
  B <- B1 * B2;

  # Construct penalty matrices
  dx <- Pars[1,6];
  Dx <- diag(nx); for (j in 1:dx) { Dx <-  diff(Dx)}
  lambdax <- Pars[1, 5];
  Px <- sqrt(lambdax) * kronecker(Dx, diag(ny));
  dy <- Pars[2, 6];
  Dy <- diag(ny); for (j in 1:dy) { Dy <-  diff(Dy)}
  lambday <- Pars[2, 5];
  Py <- sqrt(lambday) * kronecker(diag(nx), Dy);

  # Data augmentation and regression
  zx <- rep(0, ny * (nx - dx));
  zy <- rep(0, nx * (ny - dy));
  zplus <- c(Data[, 3], zx, zy);
  Bplus <- rbind(B, Px, Py);
  pcoef <- lsfit(Bplus, zplus, intercept = F)$coeff;
  pfit <- B %*% pcoef;

  # Prediction
  Bxp <- bbase(XYpred[, 1], p1[1], p1[2], p1[3], p1[4])
  Byp <- bbase(XYpred[, 2], p2[1], p2[2], p2[3], p2[4])
  B1p <- kronecker(Bxp, t(rep(1, ny)));
  B2p <- kronecker(t(rep(1, nx)), Byp);
  Bp <- B1p * B2p;
  zpred <- Bp %*% pcoef;

  P <- list(coef = pcoef, fit = pfit, pred = zpred, XYpred=XYpred)
  P
}
