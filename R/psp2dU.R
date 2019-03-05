#' Varying-coefficient penalized signal regression using P-splines.
#'
#' @description psp2dU is used to regress a (glm) response onto a two-dimensional
#'  signal such that the signal coefficients can vary over another covariate \code{t}.
#'  Aniosotripic penalization of tensor product B-splines produces a 2D coefficient surface that
#'  can be sliced at \code{t}.
#' @param y a response vector of length \code{m}, usually continuous, binary/bimomial or counts.
#' @param S a m x p Signal matrix of regressors.
#' @param S.index p-vector for index of Signal (e.g. wavelength)
#' @param t.var other (indexing) variable in Alpha surface (e.g. temperature)
#' @param Pars a matrix with 2 rows, each with P-spline parameters:
#' \code{min max nseg deg lambda pdeg}, for row and columns.
#' @param S the signal regressors, which are either "stacked" or "unfolded",
#' with dimensions (m x p1) x p2 (i.e m matrices each of p1xp2) or
#  m x (p1 x p2) (i.e. regressor matrix with m regressor rows, each with column)
#' length p1xp2, respectively.
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed.
#' @param m.binomial a vector of binomial trials having length(y). Default is 1 vector for binomial, NULL otherwise.
#' @param r.gamma a vector of gamma shape parameters. Default is 1 vector for gamma, NULL otherwise.
#' @param link link function (\code{identity, log, sqrt, logit, probit, cloglog, loglog, recipical}). Quotes are needed.
#' @param coef.plot set to T or F to display coefficient perspective surface.
#' @param image.plot set to T or F to display coefficient image surface.
#' @param x.lab(y.lab,z.lab) "character" labels for estimated coefficient surface
#' @param ridge.adj A ridge penalty tuning parameter (usually set to small value, e.g. 1e-8 to stabilize estimation).
#' @param se.bands set to T or F to produce se surfaces for plots.
#' @param wts the weight vector of \code{length(y)}. Deault is 1.
#' @param int set to T or F to include intercept term in linear predictor.
#' @param S.pred a matrix of signal inputs for  new predictions.
#' @param y.predicted a vector of responses from a cv data set (assoc. with \code{S.pred}).
#' @param t.pred: vector t covariates for new predictions
#' @param R resolution of alpha surface (e.g. 20 or 100).
#' @param x.lab,y.lab,z.lab character" labels for estimated alpha coefficient surface.

#' @return
#' \item{coef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[1,3]+Pars[1,4])}
#' of estimated P-spline coefficients for tensor surface.}
#' \item{summary.predicted}{inverse link prediction vectors, and +/- twice se bands.}
#' \item{deviance}{deviance of fit.}
#' \item{eff.df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC}
#' \item{df.resid}{approx. df resid.}
#' \item{cv}{leave-one-out standard error prediction (\code{normal, identity}).}
#' \item{cv.predicted}{standard error prediction for \code{y.predict} (normal, identity).}
#' \item{Pars}{nseg, bdeg, pord, lambda design parameters}
#' \item{Dispersion.parm}{estimate of dispersion, Dev/df.resid.}
#' \item{summary.predicted}{inverse link prediction vectors, and +/- twice se bands.}
#' \item{eta.predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press.mu}{leave-one-out prediction of mean (\code{normal, identity}).}
#' \item{bin.percent.correct}{percent correct classification based on 0.5 cut-off (\code{binomial}).}
#' @author Paul Eilers and Brian Marx
#' @references Eilers, P. H. C. and Marx, B. D. (2003). Multivariate calibration with temperature
#' interaction using two-dimensional penalized signal regression. Chemometrics
#' and Intellegent Laboratory Systems, 66, 159–174.
#'
#' @examples
#' library(fields)
#' library(JOPS)

#' psp2dU(y=rnorm(100),S.index=1:500,
#' S=matrix(rnorm(100*500),100,500),t.var=rnorm(100),
#' Pars=rbind(c(-3.5,3.5,10,3, 100, 3),c(-3.5,3.5,20,3,.1,2)),
#' se.bands=F, image.plot=T)
#' @export

"psp2dU"=
function(y, S, S.index, t.var, Pars, ridge.adj = 1e-6, S.pred = NULL, t.pred = NULL, y.predicted = NULL, R
 = 50, x.lab = "X1", y.lab = "X2", z.lab = "A.hat", coef.plot = T, image.plot = T, se.bands = T,
family = "gaussian", link = "default", m.binomial = NULL, wts = NULL, r.gamma = NULL, int = T)
{
# P-spline signal regression, with additional index t for 2-D coefficient surface
# Input:
#   y: response (of length m)
#   S: m x p Signal matrix
#   S.index: p-vector for index of Signal (e.g. wavelength)
#   t.var: other (indexing) variable in Alpha surface (e.g. temperature)
#   Pars: 2 rows with P-spline parameters: [min max nseg deg lambda pdeg]
#   ridge.adj: small ridge penalty to stabilize estimation, can be zero
#   S.pred, t.pred: (qxp signal, q-vector t) for q new predictions
#   y.predicted= a vector of responses from a cv data set (assoc. with cv M.pred)
#   R: resolution of alpha surface (e.g. 20 or 100)
#   x.lab, y.lab, z.lab: "character" labels for estimated alpha coefficient surface
#   coef.plot: T or F for perspective plot of coefficient surface
#   image.plot: T or F for image plot of coefficent surface
#   se.bands: T or F for twice standard error surface plots
#   family: "gaussian", "binomial", "poisson", "Gamma"
#   link: "logit", "probit", "log", "sqrt", "loglog", "cloglog", "identity", "inverse"
#   wts: non-negative weights (can be zero)
#   m.binomial: number of trials associated with binomial r.v. (can vary)
#   r.gamma: gamma scale parameter
#   int: T or F for intercept column of ones
#
#  Output: list with elements
#   coef: tensor product P-spline coefficients
#   summary.predicted: predicted value at new 2-D regressor locations (with 2 se bands)
#   cv: cross-validation statistic
#   aic, dev, df.residual
#   eff.dim: effective df of estimated 2-D estimated coefficient surface
#   perspective plot and image plot of estimated alpha coefficient surface
#
#  Support functions needed: bbase(), pspline2d.checker(), pspline.fitter()
#
#  Paul Eilers (2000) and Brian Marx (2001, 2002) (c)
#' @export
# Prepare bases for estimation
m <- length(y)
if(missing(wts)) {
wts <- rep(1, m)
}
if(missing(m.binomial)) {
m.binomial <- rep(1, m)
}
if(missing(r.gamma)) {
r.gamma <- rep(1, m)
}
parms <- pspline2d.checker(family, link, Pars[1, 4], Pars[2, 4], Pars[1, 6], Pars[2, 6], Pars[1,
3], Pars[2, 3], Pars[1, 5], Pars[2, 5], ridge.adj, wts)
family <- parms$family
link <- parms$link
ridge.adj <- parms$ridge.adj
wts <- parms$wts
Pars[1, 3:6] <- c(parms$ps.intervals1, parms$degree1, parms$lambda1, parms$order1)
Pars[2, 3:6] <- c(parms$ps.intervals2, parms$degree2, parms$lambda2, parms$order2)
S <- as.matrix(S)
Bx <- bbase(S.index, Pars[1, 1 ],Pars[1,2],Pars[1,3],Pars[1,4])
U <- S %*% Bx
By <- bbase(t.var, Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
n1 <- ncol(Bx)
n2 <- ncol(By)# Compute tensor products
SB1 <- kronecker(U, t(rep(1, n2)))

B1=kronecker(Bx,t(rep(1,n2)))
B2 <- kronecker(t(rep(1, n1)), By)

Q <- SB1 * B2#-----

d1 <- Pars[1, 6]
D1 <- diag(n1)
if(d1 != 0) {
for(j in 1:d1) {
D1 <- diff(D1)
}
}
lambda1 <- Pars[1, 5]
P1 <- sqrt(lambda1) * kronecker(D1, diag(n2))
d2 <- Pars[2, 6]
D2 <- diag(n2)
if(d2 != 0) {
for(j in 1:d2) {
D2 <- diff(D2)
}
}
lambda2 <- Pars[2, 5]
P2 <- sqrt(lambda2) * kronecker(diag(n1), D2)
Pen <- rbind(P1, P2)
p.ridge <- NULL
if(ridge.adj > 0) {
nix.ridge <- rep(0, n1 * n2)
p.ridge <- sqrt(ridge.adj) * diag(n1 * n2)
}
# Data augmentation and regression
z1 <- rep(0, n2 * (n1 - d1))
z2 <- rep(0, n1 * (n2 - d2))
n.col <- ncol(Q)
if(int) {
Q <- cbind(rep(1, nrow(Q)), Q)
Pen <- cbind(rep(0, nrow(Pen)), Pen)
if(ridge.adj > 0) {
p.ridge <- cbind(rep(0, nrow(p.ridge)), p.ridge)
}
}
ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b = Q, Pen, p.ridge, nix
 = c(z1, z2), nix.ridge = rep(0, n1 * n2), ridge.adj, wts)
mu <- ps.fit$mu
pcoef <- ps.fit$coef
bin.percent.correct <- NULL
if(family == "binomial") {
pcount <- 0
p.hat <- mu/m.binomial
for(ii in 1:m) {
if(p.hat[ii] > 0.5) {
count <- y[ii]
}
if(p.hat[ii] <= 0.5) {
count <- m.binomial[ii] - y[ii]
}
count <- pcount + count
pcount <- count
}
bin.percent.correct <- count/sum(m.binomial)
}
w <- ps.fit$w
e <- 1e-009
h <- hat(ps.fit$f$qr, intercept = F)[1:m]
trace <- eff.dim <- sum(h)
if(family == "binomial") {
dev <- 2 * sum((y + e) * log((y + e)/mu) + (m.binomial - y + e) * log((m.binomial - y +
e)/(m.binomial - mu)))
dispersion.parm <- 1
cv <- NULL
}
if(family == "poisson") {
dev <- 2 * sum(y * log(y + e) - y - y * log(mu) + mu)
dispersion.parm <- 1
cv <- NULL
}
if(family == "Gamma") {
dev <- -2 * sum(r.gamma * (log((y + e)/mu) - ((y - mu)/mu)))
ave.dev <- dev/m
dispersion.parm <- (ave.dev * (6 + ave.dev))/(6 + 2 * ave.dev)
cv <- NULL
}
cv <- press.mu <- press.e <- NULL
if(family == "gaussian") {
dev <- sum(ps.fit$f$residuals[1:m]^2)
dispersion.parm <- dev/(m - trace)
press.e <- ps.fit$f$residuals[1:m]/(1 - h)
cv <- sqrt(sum((press.e)^2)/(m))
press.mu <- y - press.e
}
aic <- dev + 2 * trace
w.aug <- c(w, (c(z1, z2) + 1))
A.hat <- "set coef.plot=T"
if(int) {
yint <- ps.fit$coef[1]
}
if(!int) {
yint <- NULL
}
if(coef.plot) {
# Prepare bases for estimated alpha surface
S.index. <- seq(from = Pars[1, 1], to = Pars[1, 2], length = R)
oS <- outer(rep(1, R), S.index.)
Bx. <- bbase(as.vector(oS), Pars[1, 1 ],Pars[1,2],Pars[1,3],Pars[1,4])
t.index. <- seq(from = Pars[2, 1], to = Pars[2, 2], length = R)
ot <- outer(t.index., rep(1, R))
By. <- bbase(as.vector(ot), Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
# Compute tensor products for estimated alpha surface
B1. <- kronecker(Bx., t(rep(1, n2)))
B2. <- kronecker(t(rep(1, n1)), By.)
B. <- B1. * B2.
if(int) {
A.hat <- B. %*% pcoef[2:(n.col + 1)]
}
if(!int) {
A.hat <- B. %*% pcoef
}
A.hatm <- matrix(A.hat, R, R, byrow = T)
#par(mfrow = c(1, 1))
persp(oS[1,  ], ot[, 1], A.hatm, xlab = x.lab, ylab = y.lab,zlab=z.lab,theta = 30,
phi = 30, expand = 0.5, col = "lightblue")
if(image.plot) {
image(oS[1,  ], ot[, 1], A.hatm, xlab = x.lab, ylab = y.lab,sub = "Coefficient surface")
#, nlevels = 7)
}
if(se.bands) {
half.meat <- sqrt(c(w)) * Q
meat <- t(half.meat) %*% half.meat
if(ridge.adj > 0) {
bread <- solve(meat + t(Pen) %*% Pen + t(p.ridge) %*% p.ridge)
}
if(ridge.adj == 0) {
bread <- solve(meat + t(Pen) %*% Pen)
}
half.sw <- half.meat %*% bread[, (1 + int):(n.col + int)]
var.c <- t(half.sw) %*% half.sw
half.lunch <- half.sw %*% t(B.)
ones <- 0 * y + 1
var.Ahat <- ones %*% (half.lunch * half.lunch)
stdev.Ahat <- sqrt(dispersion.parm) * t(sqrt(var.Ahat))
pivot <- 2 * stdev.Ahat
upper <- A.hat + pivot
lower <- A.hat - pivot
L.hatm <- matrix(lower, R, R, byrow = T)
U.hatm <- matrix(upper, R, R, byrow = T)
persp(oS[1,  ], ot[, 1], U.hatm, xlab = x.lab, ylab = y.lab, zlab =
"2 se Upper Surface")
if(image.plot) {
image(oS[1,  ], ot[, 1], U.hatm, xlab = x.lab, ylab = y.lab, sub =
  "2 se Upper Surface")#, nlevels = 7)
}
persp(oS[1,  ], ot[, 1], L.hatm, xlab = x.lab, ylab = y.lab, zlab =
"2 se Lower Surface")
if(image.plot) {
image(oS[1,  ], ot[, 1], L.hatm, xlab = x.lab, ylab = y.lab, sub =
  "2 se Lower Surface")#, nlevels = 7)
}
}
}
# Prediction
summary.predicted <- NULL
cv.predicted <- eta.predicted <- avediff.pred <- NULL
if(!missing(t.pred)) {
q <- length(t.pred)
Up <- as.matrix(S.pred) %*% Bx
Byp <- bbase(t.pred, Pars[2, 1 ],Pars[2,2],Pars[2,3],Pars[2,4])
B1p <- kronecker(Up, t(rep(1, n2)))
B2p <- kronecker(t(rep(1, n1)), Byp)
Qp <- B1p * B2p
if(!int) {
eta.predicted <- Qp %*% pcoef
#var.pred <- Qp %*% var.c %*% t(Qp)
}
if(int) {
#var.c <- t(bread) %*% t(half.meat) %*% half.meat %*% bread
one.xpred.b <- cbind(rep(1, q), Qp)
eta.predicted <- Qp %*% pcoef[2:length(pcoef)] + yint
#var.pred <- one.xpred.b %*% var.c %*% t(one.xpred.b)
}
#stdev.pred <- as.vector(sqrt(diag(var.pred)))
#stdev.pred <- sqrt(dispersion.parm) * stdev.pred
#pivot <- as.vector(2 * stdev.pred)
upper <- 1#eta.predicted + pivot
lower <-1# eta.predicted - pivot
summary.predicted <- cbind(lower, eta.predicted, upper)
if(!missing(y.predicted)) {
if(family == "gaussian") {
cv.predicted <- sqrt(sum((y.predicted - eta.predicted)^2)/(length(
  y.predicted)))
avediff.pred <- (sum(y.predicted - eta.predicted))/length(y.predicted)
}
}
bin.percent.correct <- NULL
if(link == "logit") {
summary.predicted <- 1/(1 + exp( - summary.predicted))
pcount <- 0
p.hat <- exp(eta.predicted)/(1 + exp(eta.predicted))
if(!missing(y.predicted)) {
for(ii in 1:length(eta.predicted)) {
  if(p.hat[ii] > 0.5) {
    count <- y.predicted[ii]
  }
  if(p.hat[ii] <= 0.5) {
    count <- 1 - y.predicted[ii]
  }
  count <- pcount + count
  pcount <- count
}
bin.percent.correct <- count/length(y.predicted)
}
}
if(link == "probit") {
summary.predicted <- apply(summary.predicted, c(1, 2), pnorm)
}
if(link == "cloglog") {
summary.predicted <- (1 - exp( - exp(summary.predicted)))
}
if(link == "loglog") {
summary.predicted <- exp( - exp( - summary.predicted))
}
if(link == "sqrt") {
summary.predicted <- summary.predicted^2
}
if(link == "log") {
summary.predicted <- exp(summary.predicted)
}
if(link == "recipical") {
summary.predd <- 1/(summary.predicted)
summary.predicted[, 1] <- summary.predd[, 3]
summary.predicted[, 3] <- summary.predd[, 1]
summary.predd <- NULL
}
summary.predicted <- as.matrix(summary.predicted)
dimnames(summary.predicted) <- list(NULL, c("-2std_Lower", "Predicted", "+2std_Upper"))
}
P <- list(coef = pcoef, Pars = Pars, cv = cv, eff.dim = eff.dim, yint = yint, int = int,
bin.percent.correct = bin.percent.correct, family, link, aic = aic, dev = dev, df.resid
 = m - trace, dispersion.parm = dispersion.parm, mu = mu, press.mu = press.mu,
summary.predicted = summary.predicted, cv.predicted = cv.predicted, eta.predicted =
eta.predicted, avediff.pred = avediff.pred, ridge.adj = ridge.adj,Q=Q,Bx=Bx,By=By)
P
}


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

