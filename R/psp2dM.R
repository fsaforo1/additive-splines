#' Two-dimensional penalized signal regression using P-splines.
#'
#' @description psp2dM is used to regress a (glm) response onto a two-dimensional
#'  signal or image, with aniosotripic penalization of
#' tensor product B-splines.
#' @param y a response vector of length \code{m}, usually continuous, binary/bimomial or counts.
#' @param p1 the row dimension of the image.
#' @param p2 the column dimension of the image.
#' @param M.type "stacked" or "unfolded"
#' @param M1.index an index of length \code{p1} for rows of regressor matrix; default is simple sequence.
#' @param M2.index an index of length \code{p2} for columns of regressor matrix.
#' @param Pars a matrix with 2 rows, each with P-spline parameters:
#' \code{min max nseg deg lambda pdeg}, for row and columns.
#' @param M The image regressors, which are either "stacked" or "unfolded",
#' with dimensions (m x p1) x p2 (i.e m matrices each of p1xp2) or
#  m x (p1 x p2) (i.e. regressor matrix with m regressor rows, each with column)
#' length p1xp2, respectively.
#' @param min(max) A scalar for the min(max) along x or y.
#'
#' @param coef.plot set to T or F to display coefficient perspective surface.
#' @param image.plot set to T or F to display coefficient image surface.
#' @param x.lab(y.lab,z.lab) "character" labels for estimated coefficient surface
#' @param ridge.adj A ridge penalty tuning parameter (usually set to small value, e.g. 1e-8 to stabilize estimation).
#' @param se.bands set to T or F to produce se surfaces for plots.
#' @param wts the weight vector of \code{length(y)}. Deault is 1.
#' @param int set to T or F to include intercept term in linear predictor.
#' @param nseg the number of evenly spaced segments between min and max.
#' @param bdeg a number for the degree of the basis, usually 1, 2, or 3.
#' @param pord the number for the order of the difference penalty, usually 1, 2, or 3.
#' @param lambda the (positive) tuning parameter for the penalty.
#' @param M.pred (stacked) \code{qp1 x p2} signal inputs  or (unfolded) \code{q x (p1 x p2)} signal
#' inputs for \code{q} new predictions.
#' @param  y.predicted a vector of responses from a cv data set (assoc. with \code{M.pred}).
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed.
#' @param m.binomial a vector of binomial trials having length(y). Default is 1 vector for binomial, NULL otherwise.
#' @param r.gamma a vector of gamma shape parameters. Default is 1 vector for gamma, NULL otherwise.
#' @param link link function (\code{identity, log, sqrt, logit, probit, cloglog, loglog, recipical}). Quotes are needed.




#' @return
#' \item{coef}{a vector of length \code{(Pars[1,3]+Pars[1,4])*(Pars[1,3]+Pars[1,4])}
#' of (unfolded) estimated P-spline coefficients for tensor surface.}
#' \item{summary.predicted}{inverse link prediction vectors, and +/- twice se bands.}
#' \item{deviance}{deviance of fit.}
#' \item{eff.df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC}
#' \item{df.resid}{approx. df resid.}
#' \item{summary.beta}{a matrix of dimension \code{px3}, containing beta, and +/- twice se bands for beta.}
#' \item{cv}{leave-one-out standard error prediction (\code{normal, identity}).}
#' \item{cv.predicted}{standard error prediction for \code{y.predict} (normal, identity).}
#' \item{Pars}{nseg, bdeg, pord, lambda design parameters}
#' \item{Dispersion.parm}{estimate of dispersion, Dev/df.resid.}
#' \item{summary.predicted}{inverse link prediction vectors, and +/- twice se bands.}
#' \item{eta.predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press.mu}{leave-one-out prediction of mean (\code{normal, identity}).}
#' \item{bin.percent.correct}{percent correct classification based on 0.5 cut-off (\code{binomial}).}
#' @author Paul Eilers and Brian Marx
#' @references Marx, B.D. and Eilers, P.H.C. (2005).
#' Multidimensional penalized signal regression, Technometrics, 47: 13-22.
#'
#' @examples
#' library(fields)
#' library(JOPS)

#' # Get the data
#' x0 = Sugar[,4:4000]
#' x0=x0-apply(x0,1,mean) #center Signal
#' y=as.vector(Sugar[,3]) #Response is Ash

#' # Inputs for two-dimensional signal regression
#' nseg=c(7,37)
#' pord=c(3,3)
#' min.=c(230,275)
#' max.=c(340,560)
#' M1.index=rev(c(340,325, 305,290,255,240,230))
#' M2.index=seq(from=275, to=560, by=.5)
#' p1=length(M1.index)
#' p2=length(M2.index)
#' int     <- T   # intercept in model

# Fit optimal model based on LOOCV
#' opt.lam=c(8858.6679, 428.1332) #Found via svcm
#' Pars.opt=rbind(c(min.[1],max.[1],nseg[1],3,opt.lam[1],pord[1]),
#'               c(min.[2],max.[2],nseg[2],3,opt.lam[2],pord[2]))
#' fit=fit.opt=psp2dM(y, x0, p1, p2, "unfolded", M1.index, M2.index, Pars.opt,
#'                   coef.plot=F, se=F, int=int, ridge.adj=1e-6)


# Plotting coefficient image
#' par(mfrow = c(1, 1))
#' set_window()
#' adj=25 # Creating mask of 1s and NAs for contriant
#' mask.=1*outer(fit$M2grid-adj, fit$M1grid,'>')
#' mask.[mask.==0]=NA
#' A.hatm=matrix(fit$bgrid%*%fit$coef[-1],length(fit$M1grid),
#'               length(fit$M2grid), byrow=T)
#' image.plot(fit$M2grid, fit$M1grid, t(A.hatm)*mask.,
#'           col = terrain.colors(100), ylab = 'Excitation (nm)',
#'           xlab = 'Emission (nm)')
#' @export

"psp2dM"<-
function(y, M, p1, p2, M.type = "stacked", M1.index = NULL, M2.index = NULL,
	Pars, ridge.adj = 1e-6, M.pred = M, y.predicted = NULL, x.lab = "X1",
	y.lab = "X2", z.lab = "A.hat", coef.plot = T, image.plot = T, se.bands
	 = T, family = "gaussian", link = "default", m.binomial = NULL, wts =
	NULL, r.gamma = NULL, int = F)
{
# P-spline regression with 2-D regressors
# Input:
#   y: response (of length m)
#   M (stacked): (mp1) x p2 matrix stacked with m regressor matrices each of p1xp2
#   M: (unfolded) m x (p1 x p2) regressor matrix with m regressor rows, each length p1xp2
#   M.type: "stacked" or "unfolded"
#   p1, p2: dimensions of regressor matrices
#   M1.index (M2.index): index for rows (cols) of regressor matrix; default is simple seq
#   Pars: 2 rows with P-spline parameters: [min max nseg deg lambda pdeg]
#   M.pred (stacked): qp1 x p2 inputs for q new predicitons
#   M.pred (unfolded): q x (p1 x p2) inputs for q new predictions
#   y.predicted= a vector of responses from a cv data set (assoc. with cv M.pred)
#   ridge.adj: small ridge penalty to stabilize estimation, can be zero
#   x.lab, y.lab, z.lab: "character" labels for estimated alpha coefficient surface
#   coef.plot: T or F for perspective plot of coefficient surface
#   image.plot: T or F for image plot of coefficient surface
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
# Paul Eilers (2000) and Brian Marx (2001, 2002) (c)
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
	parms <- pspline2d.checker(family, link, Pars[1, 4], Pars[2, 4],
	                           Pars[1,6], Pars[2, 6], Pars[1, 3],
	                           Pars[2, 3], Pars[1, 5], Pars[2, 5],
	                           ridge.adj, wts)
	family <- parms$family
	link <- parms$link
	ridge.adj <- parms$ridge.adj
	wts <- parms$wts
	Pars[1, 3:6] <- c(parms$ps.intervals1, parms$degree1, parms$lambda1,
	                  parms$order1)
	Pars[2, 3:6] <- c(parms$ps.intervals2, parms$degree2, parms$lambda2,
	                  parms$order2)

		M <- X <- as.matrix(M)
	if(M.type == "stacked") {
		p1. <- nrow(M)/m
		p2. <- ncol(M)
		if(p1 != p1. | p2 != p2.) {
			warning(paste(
				"recheck input p1 or p2: at least one dimension does not match"
				))
		}
		x <- as.vector(t(M))
		X <- matrix(x, m, p1 * p2, byrow = T)
	}
	if(missing(M1.index)) {
		M1.index <- 1:p1
	}
	if(missing(M2.index)) {
		M2.index <- 1:p2
	}
	oM1 <- outer(rep(1, p2), M1.index)
	B1 <- bbase(as.vector(oM1), Pars[1,1],Pars[1,2],Pars[1,3],Pars[1,4])
	oM2 <- outer(M2.index, rep(1, p1))
	B2 <- bbase(as.vector(oM2), Pars[2,1],Pars[2,2],Pars[2,3],Pars[2,4])
	n1 <- ncol(B1)
	n2 <- ncol(B2)	# Compute tensor products for estimated alpha surface
	B1. <- kronecker(B1, t(rep(1, n2)))
	B2. <- kronecker(t(rep(1, n1)), B2)
	B. <- B1. * B2.	# Construct penalty matrices

	#High resolution Grid Bases
	reso=200
	M1grid=seq(from=min(M1.index),to=max(M1.index),length=reso)
	M2grid=seq(from=min(M2.index),to=max(M2.index),length=reso)
	oM1g <- outer(rep(1, length(M2grid)), M1grid)
	B1grid <- bbase(as.vector(oM1g), Pars[1,1],Pars[1,2],Pars[1,3],Pars[1,4])
	oM2g <- outer(M2grid, rep(1, length(M1grid)))
	B2grid<- bbase(as.vector(oM2g), Pars[2,1],Pars[2,2],Pars[2,3],Pars[2,4])
	n1g <- ncol(B1grid)
	n2g <- ncol(B2grid)	# Compute tensor products for estimated alpha surface
	B1.g <- kronecker(B1grid, t(rep(1, n2g)))
	B2.g <- kronecker(t(rep(1, n1g)), B2grid)
	Bgrid <- B1.g * B2.g	# Construct penalty matrices
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
	Q <- X %*% B.
	n.col <- ncol(Q)
	if(int) {
		Q <- cbind(rep(1, nrow(Q)), Q)
		Pen <- cbind(rep(0, nrow(Pen)), Pen)
		if(ridge.adj > 0) {
			p.ridge <- cbind(rep(0, nrow(p.ridge)), p.ridge)
		}
	}
	ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b
		 = Q, Pen, p.ridge, nix = c(z1, z2), nix.ridge = rep(0, n1 * n2
		), ridge.adj, wts)
	mu <- ps.fit$mu
	pcoef <- ps.fit$coef
	bin.percent.correct <- bin.0percent <- bin.1percent <- NULL
	if(family == "binomial") {
		count1 <- count2 <- pcount1 <- pcount2 <- 0
		p.hat <- mu/m.binomial
		for(ii in 1:m) {
			if(p.hat[ii] > 0.5) {
				count1 <- y[ii]
				count1 <- pcount1 + count1
				pcount1 <- count1
			}
			if(p.hat[ii] <= 0.5) {
				count2 <- m.binomial[ii] - y[ii]
				count2 <- pcount2 + count2
				pcount2 <- count2
			}
		}
		bin.percent.correct <- (count1 + count2)/sum(m.binomial)
		bin.1percent <- count1/sum(m.binomial[p.hat > 0.5])
		bin.0percent <- count2/sum(m.binomial[p.hat <= 0.5])
	}
	w <- ps.fit$w
	e <- 1e-009
	h <- hat(ps.fit$f$qr, intercept = F)[1:m]
	trace <- eff.dim <- sum(h)
	if(family == "binomial") {
		dev <- 2 * sum((y + e) * log((y + e)/mu) + (m.binomial - y + e) *
			log((m.binomial - y + e)/(m.binomial - mu)))
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
	cv <- press.mu <- press.e <- var.c <- NULL
	if(family == "gaussian") {
		dev <- sum(ps.fit$f$residuals[1:m]^2)
		dispersion.parm <- dev/(m - trace)
		press.e <- ps.fit$f$residuals[1:m]/(1 - h)
		cv <- sqrt(sum((press.e)^2)/(m))
		press.mu <- y - press.e
	}
	aic <- dev + 2 * trace
	w.aug <- c(w, (c(z1, z2) + 1))
	if(int) {
		A.hat <- B. %*% pcoef[2:(n.col + 1)]
		yint <- ps.fit$coef[1]
	}
	if(!int) {
		yint <- NULL
		A.hat <- B. %*% pcoef
	}
	A.hatm <- matrix(A.hat, p1, p2, byrow = T)
	if(coef.plot) {
		i.1 <- M1.index
		in1 <- 1:length(M1.index)
		i.2 <- M2.index
		in2 <- 1:length(M2.index)
		if(length(M1.index) > 100) {
			i.1 <- round(seq(from = min(M1.index), to = max(
				M1.index), length = 100))
			in1 <- round(seq(from = 1, to = p1, length = 50))
		}
		if(length(M2.index) > 100) {
			i.2 <- round(seq(from = min(M2.index), to = max(
				M2.index), length = 100))
			in2 <- round(seq(from = 1, to = p2, length = 50))
		}
		persp(M2.index[in2], (M1.index[in1]), t(A.hatm[in1, in2]), xlab
			 = y.lab, ylab = x.lab, zlab = z.lab)
		if(image.plot) {
			image(M2.index[in2], M1.index[in1], t(A.hatm[in1, in2]),
				xlab = x.lab, ylab = y.lab, sub = "A hat")
		}
		if(se.bands) {
			half.meat <- sqrt(c(w)) * Q
			meat <- t(half.meat) %*% half.meat
			if(ridge.adj > 0) {
				bread <- solve(meat + t(Pen) %*% Pen + t(
				  p.ridge) %*% p.ridge)
			}
			if(ridge.adj == 0) {
				bread <- solve(meat + t(Pen) %*% Pen)
			}
			half.sw <- half.meat %*% bread[, (1 + int):(n.col + int
				)]
			var.c <- t(half.sw) %*% half.sw
			half.lunch <- half.sw %*% t(B.)
			ones <- 0 * y + 1
			var.Ahat <- ones %*% (half.lunch * half.lunch)
			stdev.Ahat <- sqrt(dispersion.parm) * t(sqrt(var.Ahat))
			pivot <- 2 * stdev.Ahat
			upper <- A.hat + pivot
			lower <- A.hat - pivot
			L.hatm <- matrix(lower, p1, p2, byrow = T)
			U.hatm <- matrix(upper, p1, p2, byrow = T)
			persp(M2.index[in2], M1.index[in1], t(U.hatm[in1, in2]),
				xlab = x.lab, ylab = y.lab, zlab =
				"2 se Upper Surface")
			if(image.plot) {
				image(M2.index[in2], M1.index[in1], t(U.hatm[
				  in1, in2]), xlab = x.lab, ylab = y.lab, sub
				   = "2 se Upper Surface")
			}
			persp(M2.index[in2], M1.index[in1], t(L.hatm[in1, in2]),
				xlab = x.lab, ylab = y.lab, zlab =
				"2 se Lower Surface")
			if(image.plot) {
				image(M2.index[in2], M1.index[in1], t(L.hatm[
				  in1, in2]), xlab = x.lab, ylab = y.lab, sub
				   = "2 se Lower Surface")
			}
		}
	}
	summary.predicted <- NULL
	cv.predicted <- eta.predicted <- avediff.pred <- NULL
	if(!missing(M.pred)) {
		half.meat <- sqrt(c(w)) * Q
		meat <- t(half.meat) %*% half.meat
		if(ridge.adj > 0) {
			bread <- solve(meat + t(Pen) %*% Pen + t(p.ridge) %*%
				p.ridge)
		}
		if(ridge.adj == 0) {
			bread <- solve(meat + t(Pen) %*% Pen)
		}
		half.sw <- half.meat %*% bread[, (1 + int):(n.col + int)]
		var.c <- t(half.sw) %*% half.sw
		q <- nrow(M.pred)
		X.p <- as.matrix(M.pred)
		if(M.type == "stacked") {
			M.pred <- as.matrix(M.pred)
			q <- nrow(M.pred)/p1
			X.p <- matrix(x, q, p1 * p2, byrow = T)
		}
		if(!int) {
			eta.predicted <- X.p %*% as.vector(A.hat)
			var.pred <- X.p %*% B. %*% var.c %*% t(X.p %*% B.)
		}
		if(int) {
			var.c <- t(bread) %*% t(half.meat) %*% half.meat %*%
				bread
			one.xpred.b <- cbind(rep(1, q), (X.p %*% B.))
			eta.predicted <- X.p %*% A.hat + yint
			var.pred <- one.xpred.b %*% var.c %*% t(one.xpred.b)
		}
		stdev.pred <- as.vector(sqrt(diag(var.pred)))
		stdev.pred <- sqrt(dispersion.parm) * stdev.pred
		pivot <- as.vector(2 * stdev.pred)
		upper <- eta.predicted + pivot
		lower <- eta.predicted - pivot
		summary.predicted <- cbind(lower, eta.predicted, upper)
		if(!missing(y.predicted)) {
			if(family == "gaussian") {
				cv.predicted <- sqrt(sum((y.predicted -
				  eta.predicted)^2)/(length(y.predicted)))
				avediff.pred <- (sum(y.predicted -
				  eta.predicted))/length(y.predicted)
			}
		}
		bin.percent.correct <- bin.0percent <- bin.1percent <- NULL
		if(link == "logit") {
			summary.predicted <- 1/(1 + exp( - summary.predicted))
			count1 <- count2 <- pcount1 <- pcount2 <- 0
			p.hat <- exp(eta.predicted)/(1 + exp(eta.predicted))
			if(!missing(y.predicted)) {
				y.predicted <- as.vector(y.predicted)
				for(ii in 1:length(eta.predicted)) {
				  if(p.hat[ii] > 0.5) {
				    count1 <- y.predicted[ii]
				    count1 <- pcount1 + count1
				    pcount1 <- count1
				  }
				  if(p.hat[ii] <= 0.5) {
				    count2 <- 1 - y.predicted[ii]
				    count2 <- pcount2 + count2
				    pcount2 <- count2
				  }
				}
				bin.percent.correct <- (count1 + count2)/length(
				  y.predicted)
				bin.1percent <- count1/sum(y.predicted)
				bin.0percent <- count2/(length(y.predicted) -
				  sum(y.predicted))
			}
		}
		if(link == "probit") {
			summary.predicted <- apply(summary.predicted, c(1, 2),
				pnorm)
		}
		if(link == "cloglog") {
			summary.predicted <- (1 - exp( - exp(summary.predicted)
				))
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
		dimnames(summary.predicted) <- list(NULL, c("-2std_Lower",
			"Predicted", "+2std_Upper"))
	}
	P <- list(coef = pcoef, Pars = Pars, yint = yint, int = int, family,
		link, dev = dev, aic = aic, bin.percent.correct =
		bin.percent.correct, bin.0 = bin.0percent, bin.1 = bin.1percent,
		df.resid = m - trace, dispersion.parm = dispersion.parm, mu =
		mu, press.mu = press.mu, summary.predicted = summary.predicted,
		eta.predicted = eta.predicted, avediff.pred = avediff.pred,
		ridge.adj = ridge.adj, cv = cv, cv.predicted = cv.predicted,
		eff.dim = eff.dim,Q=Q,b=B.,h=h,bgrid=Bgrid,M1grid=M1grid,M2grid=M2grid)
	P
}
