#' @export
# psp2d.der Support function needed for VSISR (Marx, 2015)
"psp2d.der"<-function(Data, Pars, XYpred, Plots=F,dPlots=F, pPlots=T,iPlots=F, sPlots=F)
{
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


  Bx <- bbase(as.vector(Data[, 1]), Pars[1, 1 ],Pars[1,2],Pars[1,3],Pars[1,4])
  By <- bbase(as.vector(Data[, 2]), Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])

  m=nrow(Bx)
  nx <- ncol(Bx)
  ny <- ncol(By)# Compute tensor products
  B1 <- kronecker(Bx, t(rep(1, ny)))
  B2 <- kronecker(t(rep(1, nx)), By)
  B <- B1 * B2# Construct penalty matrices
  dx <- Pars[1, 6]
  Dx <- diag(nx)
  for(j in 1:dx) {
    Dx <- diff(Dx)
  }
  lambdax <- Pars[1, 5]
  Px <- sqrt(lambdax) * kronecker(Dx, diag(ny))
  dy <- Pars[2, 6]
  Dy <- diag(ny)
  for(j in 1:dy) {
    Dy <- diff(Dy)
  }
  lambday <- Pars[2, 5]
  Py <- sqrt(lambday) * kronecker(diag(nx), Dy)# Data augmentation and regression
  zx <- rep(0, ny * (nx - dx))
  zy <- rep(0, nx * (ny - dy))
  zplus <- c(Data[, 3], zx, zy)
  Bplus <- rbind(B, Px, Py)


  f <- lsfit(Bplus, zplus,intercept = F)
  pcoef. <- f$coef

  pcoef=solve(t(Bplus)%*%Bplus)%*%t(Bplus)%*%zplus
  #pcoef <- solve(Bplus, zplus)



  pfit <- B %*% pcoef# Prediction

  h <- hat(f$qr, intercept = F)[1:m]
  press.e <- f$residuals[1:m]/(1 - h)
  cv <- sqrt(sum((press.e)^2)/(m))



  R=100

  # Prepare bases for estimated surfaces
  eta.index. <- seq(from = Pars[1, 1], to = Pars[1, 2], length = R)
  oeta <- outer(rep(1, R), eta.index.)
  Bx. <- bbase(as.vector(oeta),Pars[1, 1 ],Pars[1,2],Pars[1,3],Pars[1,4])
  t.index. <- seq(from = Pars[2, 1], to = Pars[2, 2], length = R)
  ot <- outer(t.index., rep(1, R))
  By. <- bbase(as.vector(ot), Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
  # Compute tensor products for estimated alpha surface
  B1. <- kronecker(Bx., t(rep(1, ncol(By.))))
  B2. <- kronecker(t(rep(1, ncol(Bx.))), By.)
  B. <- B1. * B2.
  pfit.. <- B. %*% pcoef
  y.lab="Temperature (C)"
  x.lab="Linear predictor"
  z.lab= " "


  fit.matrix <- matrix(pfit.., R, R, byrow = T)
  if(Plots){
    if(pPlots){
      persp(oeta[1,  ], ot[, 1], fit.matrix, xlab = x.lab, ylab = y.lab, zlab = z.lab,
            theta = 30, phi = 30, expand = 0.5, col = "lightblue")}
    if(iPlots){
      image.plot(oeta[1,  ], ot[, 1], fit.matrix, xlab = x.lab, ylab = y.lab, main='Link surface')}
    if(sPlots){
      matplot(eta.index.,fit.matrix[,seq(1,R, length=6)],type='l',col=(1:10),
              lty=c(1:10), ylab=' ', xlab='Linear predictor', main='Link slices, by temp')
      abline(0,1,lty=2, col=1, lwd=2)}
  }


  Bxp <- bbase(XYpred[, 1], Pars[1, 1 ],Pars[1,2],Pars[1,3],Pars[1,4])
  Byp <- bbase(XYpred[, 2], Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
  B1p <- kronecker(Bxp, t(rep(1, ny)))
  B2p <- kronecker(t(rep(1, nx)), Byp)
  Bp <- B1p * B2p
  zpred <- Bp %*% pcoef


  ##########################
  pcoef.m=matrix(pcoef,nx,ny, byrow=T)

  #d.pcoef.m=t(diff(t(pcoef.m)))
  d.pcoef.m=(diff((pcoef.m)))

  Pars1.d=Pars[1,]
  Pars1.d[4]=Pars1.d[4]-1

  Bx.d <- bbase(Data[, 1], Pars1.d[1], Pars1.d[2], Pars1.d[3],Pars1.d[4])
  By <- bbase(Data[, 2],Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
  nx <- ncol(Bx.d)
  ny <- ncol(By)# Compute tensor products
  B1d <- kronecker(Bx.d, t(rep(1, ny)))
  B2d <- kronecker(t(rep(1, nx)), By)
  B.d <- B1d * B2d
  d.pcoef=c(t(d.pcoef.m))/((Pars[1,2]-Pars[1,1])/(1*Pars[1,3]))

  #d.pcoef=c((d.pcoef.m))
  d.fit=B.d%*%d.pcoef # /((Pars[1,2]-Pars[1,1])/(1*Pars[1,3]))

  Pars1.d=Pars[1,]
  Pars1.d[4]=Pars1.d[4]-1
  #R=100
  # Prepare bases for estimated DERIVATIVE surfaces
  eta.index. <- seq(from = Pars[1, 1], to = Pars[1, 2], length = R)
  oeta <- outer(rep(1, R), eta.index.)
  Bx. <- bbase(as.vector(oeta), Pars1.d[1], Pars1.d[2], Pars1.d[3],Pars1.d[4])
  t.index. <- seq(from = Pars[2, 1], to = Pars[2, 2], length = R)
  ot <- outer(t.index., rep(1, R))
  By. <- bbase(as.vector(ot), Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
  # Compute tensor products for estimated alpha surface
  B1. <- kronecker(Bx., t(rep(1, ncol(By.))))
  B2. <- kronecker(t(rep(1, ncol(Bx.))), By.)
  B. <- B1. * B2.
  pfit... <- B. %*% d.pcoef # /((Pars[1,2]-Pars[1,1])/(1*Pars[1,3]))

  y.lab="Temperature (C)"
  x.lab="Linear predictor"
  z.lab= " "


  dfit.matrix <- matrix(pfit..., R, R, byrow = T)
  if(dPlots){
    if(pPlots){
      persp(oeta[1,  ], ot[, 1], dfit.matrix, xlab = x.lab, ylab = y.lab, zlab = z.lab,
            theta = 30, phi = 30, col = "lightblue")}
    if(iPlots){
      image.plot(oeta[1,  ], ot[, 1], dfit.matrix, xlab = x.lab, ylab = y.lab)}
    if(sPlots){
      matplot(eta.index.,dfit.matrix[,seq(1,R, length=6)],type='l',col=(1:10),
              lty=c(1:10), ylab=' ', xlab='Linear predictor')}
  }

  #B.der<-bbase(x, obj$xmin, obj$xmax, obj$nseg, obj$bdeg-1)
  #    alpha.der<-diff(obj$pcoef)
  #    out<-as.vector(B.der%*%alpha.der)/((obj$xmax-obj$xmin)/obj$nseg)
  ##########################
  P <- list(coef = pcoef, fit = pfit,d.fit=d.fit, pred = zpred,cv=cv,
            fit.matrix.grid=fit.matrix, dfit.matrix.grid=dfit.matrix)
  P
}

#' @export
# Support function needed for SISR (Eilers, Li, Marx, 2009)
`pnormal.der` <-
  function(x, y, nseg, bdeg, pord, lambda, plot = F, se = F) #, xpred)
  {
    # Function pnormal: smooths scatterplot data with P-splines.
    # Input:
    #   x = abcissae of data
    #   y = response
    #   nseg = number of intervals for B-splines
    #   bdeg = degree of B-splines
    #   pord = order of difference penalty
    #   lambda = smoothness parameter
    #   plot = plot parameter (T of F)
    #   se = plot parameter (T or F)

    # Output: an object of class "pspfit" with the following fields
    #   bdeg = degree of B-splines
    #   cv = cross-validation sum of squares
    #   ed.resid = effective degrees of freedom residuals
    #   effdim = effective dimension P-spline model
    #   family = "gaussian" (like glm object)
    #   lambda = smoothing parameter
    #   link = "identity" (like glm object)
    #   muhat = expected values for y (at x)
    #   mse = standard deviation of errors
    #   nseg = number of B-spline segments on domain from xmin to xmax)
    #   pord = order of difference penalty
    #   x = x as input
    #   xgrid = x grid used for plotting curve
    #   xmin = left boundary of B-spline domain
    #   xmax = right boundary of B-spline domain
    #   y = y as input
    #   ygrid = computed curve on x grid

    #
    # Side effect: a plot of (x,y) and the estimated curve (if plot = T) with twice se bands (if se=T).

    #
    # Paul Eilers and Brian Marx, 2003 (c)
    #

    # Compute B-spline basis

    m <- length(x)
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.5 * (xr - xl)
    xmin <- xl - 0.5 * (xr - xl)
    B <- bbase(x, xmin, xmax, nseg, bdeg)

    # Construct penalty stuff
    n <- dim(B)[2]
    P <- sqrt(lambda) * ndiff(n, pord)
    nix <- rep(0, n - pord)

    # Fit
    if(lambda == 0) {
      f <- lsfit(B, y, intercept = F)
    }
    if(lambda > 0) {
      f <- lsfit(rbind(B, P), c(y, nix), intercept = F)
    }
    h <- hat(f$qr)[1:m]
    beta <- as.vector(f$coef)
    mu <- B %*% beta

    # Cross-validation and dispersion
    r <- (y - mu)/(1 - h)
    cv <- sqrt((sum(r^2))/m)
    s <- sqrt(sum((y - mu)^2)/(m - sum(h)))

    # Compute curve on grid
    u <- seq(xl-.1, xr+.1, length = 100)
    Bu <- bbase(u, xmin-.1, xmax+.1, nseg, bdeg)
    zu <- Bu %*% as.vector(f$coef)
    #Bu2<-bbase(xpred, xmin, xmax, nseg, bdeg)
    #fit.xpred<-Bu2 %*% as.vector(f$coef)

    #Derivative
    B.der <- bbase(x, xmin, xmax, nseg, bdeg-1)
    alpha.der <- diff(f$coef)
    der <- B.der%*%alpha.der
    #B.der2 <- bbase(xpred, xmin, xmax, nseg, bdeg-1)
    #der.xpred<-B.der2 %*% alpha.der

    # Compute derivative on grid
    #u <- seq(xl, xr, length = 100)
    #Bu.der <- bbase(u, xmin, xmax, nseg, bdeg-1)
    #zu.der <- Bu.der %*% as.vector(diff(f$coef))


    # Plot data and fit
    if(plot) {
      plot(x, y)
      lines(u, zu, col = 2)
      lines(u, zu.der, col = 4)
      if(se) {
        varf <- diag(Bu %*% solve(t(B) %*% B + t(P) %*% P) %*% t(Bu))
        sef <- s * sqrt(varf)
        upperu <- zu + 2 * sef
        loweru <- zu - 2 * sef
        lines(u, upperu, lty = 3, col = 3)
        lines(u, loweru, lty = 3, col = 3)
      }
    }

    # Return list
    pp <- list(x = x, y = y, muhat = mu, nseg = nseg, xmin = xmin, bdeg = bdeg, pord
               = pord, lambda = lambda, xgrid = u, ygrid = zu, cv = cv, effdim = sum(h
               ), ed.resid = m - sum(h), family = "gaussian", link = "identity", sqrt.mse =
                 s, pcoef=beta, xmin=xmin, xmax=xmax) #, fit.xpred=fit.xpred, der.xpred=der.xpred)
    class(pp) <- "pspfit"
    pp
  }

# Support function for SISR
`predict.pnormal` <-
  function(obj,x,der){
    #der can only be either 0 or 1
    if (der==0){
      bu<-bbase(x,obj$xmin,obj$xmax,obj$nseg,obj$bdeg)
      out<-as.vector(bu%*%obj$pcoef)
    }
    if (der==1){
      B.der<-bbase(x, obj$xmin, obj$xmax, obj$nseg, obj$bdeg-1)
      alpha.der<-diff(obj$pcoef)
      out<-as.vector(B.der%*%alpha.der)/((obj$xmax-obj$xmin)/obj$nseg)
    }
    out
  }

