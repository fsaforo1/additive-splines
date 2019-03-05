#' Smooth signal (multivariate calibration) regression using P-splines.
#' @description Smooth signal (multivariate calibration) regression using P-splines.
#' Support Functions: \code{pspline.fitter()} and \code{pspline.checker()}.
#'
#' @author Brian Marx
#'
#' @param response a response vector, usually continuous, binomial or count data.
#' @param X  a matrix of continuous regressor with nrow(X)=length(y), often
#' a discrete digitization of a signal or histogram or time series.
#' @param x.index a vector to of length \code{ncol(X)=p}, associated with the
#' ordering index of the signal. Default is \code{1:ncol(X)}.
#' @param family the response distribution, e.g.
#' \code{"gaussian", "binomial", "poisson", "Gamma"} distribution. Quotes are needed.
#' @param m.binomial a vector of binomial trials having length(y). Default is 1 vector for binomial, NULL otherwise.
#' @param r.gamma a vector of gamma shape parameters. Default is 1 vector for gamma, NULL otherwise.
#' @param link link function (\code{identity, log, sqrt, logit, probit, cloglog, loglog, recipical}). Quotes are needed.
#' @param nseg the number of evenly spaced segments between \code{xl} and \code{xr}.
#' @param bdeg the degree of the basis, usually 1, 2, or 3 (defalult).
#' @param pord the order of the difference penalty, usually 1, 2, or 3 (defalult).
#' @param lambda the (positive) tuning parameter for the penalty.
#' @param coef.plot set to T or F to display default plots.
#' @param ridge.adj A ridge penalty tuning parameter (usually set to small value, e.g. 1e-8 to stabilize estimation).
#' @param se.bands set to T or F to produce se bands in plots.
#' @param wts the weight vector of \code{length(y)}. Deault is 1.
#' @param int set to T or F to include intercept term in linear predictor.
#' @param x.predicted a matrix of external signals to yield exernal prediction.
#' @param y.predicted a vector of responses associated
#' with \code{x.predicted} which are used to calculated standard error of external prediction. Default is NULL.

#' @return
#' \item{coeff}{a vector with \code{length(n)} of estimated P-spline coefficients.}
#' \item{mu}{a vector with \code{length(m)} of estimated means.}
#' \item{eta}{a vector of \code{length(m)} of estimated linear predictors.}
#' \item{b}{the B-spline basis (for the coefficients), with dimension \code{p}x\code{n}.}
#' \item{deviance}{deviance of fit.}
#' \item{eff.df}{the approximate effective dimension of fit.}
#' \item{aic}{AIC}
#' \item{df.resid}{approx. df resid.}
#' \item{summary.beta}{a matrix of dimension \code{px3}, containing beta, and +/- twice se bands for beta.}
#' \item{cv}{leave-one-out standard error prediction (\code{normal, identity}).}
#' \item{cv.predicted}{standard error prediction for \code{y.predict} (normal, identity).}
#' \item{nseg, bdeg, pord, lambda}{design parameters}
#' \item{Dispersion.parm}{estimate of dispersion, Dev/df.resid.}
#' \item{summary.predicted}{inverse link prediction vectors, and +/- twice se bands.}
#' \item{eta.predicted}{estimated linear predictor of \code{length(y)}.}
#' \item{press.mu}{leave-one-out prediction of mean (\code{normal, identity}).}
#' \item{bin.percent.correct}{percent correct classification based on 0.5 cut-off (\code{binomial}).}
#' @references  Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression for sampled signals and
#'        curves: A P-spline approach. Technometrics, 41(1): 1-13.
#'
#' @examples
#' library(JOPS)
#' # Get the data
#' library(fds)
#' data(nirc)
#' iindex=nirc$x
#' X=nirc$y
#' sel= 50:650 #1200 <= x & x<= 2400
#' X=X[sel,]
#' iindex=iindex[sel]
#' dX=diff(X)
#' diindex=iindex[-1]
#' y=as.vector(labc[1,1:40])
#' oout=23
#' dX=t(dX[,-oout])
#' y=y[-oout]
#' fit1=signal.fit(y,diindex,dX, nseg=25,lambda=.0001,coef.plot=T)
#' title(main='25 B-spline segments with tuning=0.0001')
#' names(fit1)
#'
#' @export
"signal.fit"<-
function(response, x.index, x.signal, m.binomial = NULL, nseg = 8,
	bdeg = 3, pord = 3, wts = NULL, link = "default", family =
	"gaussian", r.gamma = NULL, lambda = 0, y.predicted = NULL, x.predicted
	 = NULL, ridge.adj = 0, int = T, coef.plot = T, se.bands = T)
{
# Function signal.fit: smooths signal (multivariate calibration) beta's using P-splines.
# Input: x.index= abcissae of spectra (1:p).
# Input: x.signal= (n X p) explanatory variable signal matrix (p >> n possible).
# Input: response= response variable.
# Input: family=gaussian, binomial, poisson, Gamma distribution.
# Input: m.binomial=vector of binomial trials. Default is 1 vector.
# Input: r.gamma=vector of gamma shape parameters. Default is 1 vector.
# Input: link= link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical).
# Input: nseg= number of intervals for B-splines. Default=8.
# Input: bdeg= degree of B-splines. Default=3.
# Input: pord= order of difference penalty. Default=3.
# Input: lambda= smoothness regulalizing parameter ( >= 0). Default=0.
# Input: x.predicted=a matrix of row (original) signals for prediction and twice stderr limits.
# Input: y.predicted= a vector of responses from a cv data set (assoc. with cv x.predicted).
# Input: ridge.adj= a small positive constant to help stabilize linear dependencies among B-splines.
# Result: a plot of smoothed beta's and twice stderr.
# Output: A list: including, AIC= deviance + 2*trace(Hat), dispers.parm, etc.
#
# Support Functions: pspline.fitter() and pspline.checker()
#
# References:
# Marx, B.D. and Eilers, P.H.C. (1999). Generalized linear regression for sampled signals and
#        curves: A P-spline approach. Technometrics, 41(1): 1-13.
# Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties (with comments
#        and rejoinder). Statistical Science, 11(2): 89-121.
#
# (c) 1995 Paul Eilers & Brian Marx
#
	y <- response
	x <- x.index
	n <- length(y)
	if(missing(wts)) {
		wts <- rep(1, n)
	}
	parms <- pspline.checker(family, link, bdeg, pord,
	                         nseg,
		lambda, ridge.adj, wts)
	family <- parms$family
	link <- parms$link
	q <- parms$bdeg
	d <- parms$pord
	ridge.adj <- parms$ridge.adj
	lambda <- parms$lambda
	nseg <- parms$nseg
	wts <- parms$wts
	if(missing(m.binomial)) {
		m.binomial <- rep(1, n)
	}
	if(missing(r.gamma)) {
		r.gamma <- rep(1, n)
	}
	xl <- min(x)
	xr <- max(x)
	xmax <- xr + 0.01 * (xr - xl)
	xmin <- xl - 0.01 * (xr - xl)
	dx <- (xmax - xmin)/nseg
	knots <- seq(xmin - q * dx, xmax + q * dx, by = dx)
	b <- bbase(x, xl, xr, nseg,3)
	  #spline.des(knots, x, q + 1, 0 * x)$design
	n.col <- ncol(b)
	if(d < 0) {
		d <- min(3, (n.col - 1))
		warning(paste("penalty order cannot be negative: have used", d)
			)
	}
	if((d - n.col + 1) > 0) {
		d <- n.col - 1
		warning(paste("penalty order was too large: have used", d))
	}
	p.ridge <- NULL
	if(ridge.adj > 0) {
		nix.ridge <- rep(0, n.col)
		p.ridge <- sqrt(ridge.adj) * diag(rep(1, n.col))
	}
	p <- diag(n.col)
	if(d != 0) {
		for(j in 1:d) {
			p <- diff(p)
		}
	}
	p <- sqrt(lambda) * p
	nix <- rep(0, n.col - d)
	x.signal <- as.matrix(x.signal)
	xb <- x.signal %*% as.matrix(b)
	if(int) {
		xb <- cbind(rep(1, n), xb)
		p <- cbind(rep(0, nrow(p)), p)
		if(ridge.adj > 0) {
			p.ridge <- cbind(rep(0, nrow(p.ridge)), p.ridge)
		}
	}
	ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b
		 = xb, p, p.ridge, nix, nix.ridge, ridge.adj, wts)
	mu <- ps.fit$mu
	coef <- ps.fit$coef
	eta=ps.fit$eta
	bin.percent.correct <- NULL
	if(family == "binomial") {
		pcount <- 0
		p.hat <- mu/m.binomial
		for(ii in 1:n) {
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
	h <- hat(ps.fit$f$qr, intercept = F)[1:n]
	trace <- sum(h)
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
		ave.dev <- dev/n
		dispersion.parm <- (ave.dev * (6 + ave.dev))/(6 + 2 * ave.dev)
		cv <- NULL
	}

	cv <- press.mu <- press.e <- NULL

	if(family == "gaussian") {
	  dev <- sum(ps.fit$f$residuals[1:n]^2)
		dispersion.parm <- dev/(n - trace)
		press.e <- ps.fit$f$residuals[1:n]/(1 - h)
		cv <- sqrt(sum((press.e)^2)/(n))
		press.mu <- y - press.e
	}
	aic <- dev + 2 * trace
	w.aug <- c(w, (nix + 1))
	if(int) {
		beta <- b %*% (as.vector(ps.fit$f$coef)[2:(n.col + 1)])
		yint <- ps.fit$coef[1]
	}
	if(!int) {
		yint <- NULL
		beta <- b %*% as.vector(ps.fit$f$coef)
	}
	half.meat <- sqrt(c(w)) * xb
	meat <- t(half.meat) %*% half.meat
	if(ridge.adj > 0) {
		bread <- solve(meat + t(p) %*% p + t(p.ridge) %*% p.ridge)
	}
	if(ridge.adj == 0) {
		bread <- solve(meat + t(p) %*% p)
	}
	half.sw <- half.meat %*% bread[, (1 + int):(n.col + int)]
	var.c <- t(half.sw) %*% half.sw
	half.lunch <- half.sw %*% t(b)
	ones <- 0 * y + 1
	var.beta <- ones %*% (half.lunch * half.lunch)
	stdev.beta <- sqrt(dispersion.parm) * t(sqrt(var.beta))
	pivot <- 2 * stdev.beta
	upper <- beta + pivot
	lower <- beta - pivot
	summary.beta <- cbind(lower, beta, upper)
	if(coef.plot) {
		if(se.bands) {
			matplot(x.index, summary.beta, type = "l", lty = c(3, 1,
				3), col = rep(1, 3), xlab = "Coefficient Index",
				ylab = "P-spline Coefficient")
		}
		if(!se.bands) {
			plot(x.index, summary.beta[, 2], col = 1, type = "l",
				lty = 1, xlab = "Coefficient Index", ylab =
				"P-spline Coefficient")
		}
	}
	summary.predicted <- NULL
	cv.predicted <- eta.predicted <- NULL
	if(!missing(x.predicted)) {
		x.predicted <- as.matrix(x.predicted)
		if(!int) {
			if(ncol(x.predicted) > 1) {
				eta.predicted <- x.predicted %*% beta
				var.pred <- x.predicted %*% b %*% var.c %*% t(
				  x.predicted %*% b)
			}
			if(ncol(x.predicted) == 1) {
				eta.predicted <- t(x.predicted) %*% beta
				var.pred <- t(x.predicted) %*% b %*% var.c %*%
				  t(b) %*% x.predicted
			}
		}
		if(int) {
			var.c <- t(bread) %*% t(half.meat) %*% half.meat %*%
				bread
			dim.xp <- nrow(x.predicted)
			if(ncol(x.predicted) > 1) {
				one.xpred.b <- cbind(rep(1, dim.xp), (
				  x.predicted %*% b))
				eta.predicted <- x.predicted %*% beta + yint
				var.pred <- one.xpred.b %*% var.c %*% t(
				  one.xpred.b)
			}
			if(ncol(x.predicted) == 1) {
				one.xpred.b <- cbind(1, t(x.predicted) %*% b)
				eta.predicted <- t(x.predicted) %*% beta + yint
				var.pred <- (one.xpred.b) %*% var.c %*% t(
				  one.xpred.b)
			}
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
				bin.percent.correct <- count/length(y.predicted
				  )
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
	llist <- list()
	llist$b <- b
	llist$coef <- coef
	llist$y.intercept <- yint
	llist$int <- int
	llist$press.mu <- press.mu
	llist$bin.percent.correct <- bin.percent.correct
	llist$family <- family
	llist$link <- link
	llist$nseg <- nseg
	llist$pord <- d
	llist$bdeg <- q
	llist$lambda <- lambda
	llist$aic <- aic
	llist$deviance <- dev
	llist$eff.df <- trace - 1
	llist$df.resid <- n - trace + 1
	llist$bin.percent.correct <- bin.percent.correct
	llist$dispersion.param <- dispersion.parm
	llist$summary.predicted <- summary.predicted
	llist$eta.predicted <- eta.predicted
	llist$cv.predicted <- cv.predicted
	llist$cv <- cv
	llist$mu <- mu
	llist$eta=eta
	llist$summary.beta<-summary.beta
	llist
}

#' @export
"pspline.fitter"<-
  function(family, link, n.col, m.binomial, r.gamma, y, b, p, p.ridge, nix,
           nix.ridge, ridge.adj, wts, ...)
  {
    coef.est <- rep(1, ncol(b))
    if(family == "binomial") {
      mu <- (y + 0.5 * m.binomial)/2
    }
    if(family == "Gamma" || family == "poisson") {
      mu <- (y + 3)
    }
    if(family == "gaussian") {
      mu <- rep(mean(y), length(y))
    }
    it <- 0
    repeat {
      if(it == 0) {
        if(link == "identity") {
          eta <- mu
        }
        if(link == "log") {
          eta <- log(mu)
        }
        if(link == "sqrt") {
          eta <- sqrt(mu)
        }
        if(link == "logit") {
          eta <- log(mu/(m.binomial - mu))
        }
        if(link == "recipical") {
          eta <- 1/mu
        }
        if(link == "probit") {
          eta <- qnorm(mu/m.binomial)
        }
        if(link == "cloglog") {
          eta <- log( - log(1 - mu/m.binomial))
        }
        if(link == "loglog") {
          eta <-  - log( - log(mu/m.binomial))
        }
      }
      it <- it + 1
      if(it > 50)
        break
      if(link == "identity") {
        mu <- eta
        h.prime <- 1
      }
      if(link == "log") {
        mu <- exp(eta)
        h.prime <- mu
      }
      if(link == "sqrt") {
        mu <- eta^2
        h.prime <- 2 * eta
      }
      if(link == "logit") {
        mu <- m.binomial/(1 + exp( - eta))
        h.prime <- mu * (1 - mu/m.binomial)
      }
      if(link == "recipical") {
        mu <- 1/eta
        h.prime <-  - (mu^2)
      }
      if(link == "probit") {
        mu <- m.binomial * pnorm(eta)
        h.prime <- m.binomial * dnorm(eta)
      }
      if(link == "cloglog") {
        mu <- m.binomial * (1 - exp( - exp(eta)))
        h.prime <- (m.binomial) * exp(eta) * exp( - exp(eta))
      }
      if(link == "loglog") {
        mu <- m.binomial * exp( - exp( - eta))
        h.prime <- m.binomial * exp( - eta) * exp( - exp( - eta
        ))
      }
      if(family == "gaussian") {
        w <- rep(1, length(y))
      }
      if(family == "poisson") {
        w <- h.prime^2/mu
      }
      if(family == "binomial") {
        w <- h.prime^2/(mu * (1 - mu/m.binomial))
      }
      if(family == "Gamma") {
        w <- (r.gamma * h.prime^2)/mu^2
      }
      u <- (y - mu)/h.prime + eta
      if(ridge.adj > 0) {
        f <- lsfit(rbind(b, p, p.ridge), c(u, nix, nix.ridge),
                   wt = c(wts, nix + 1, nix.ridge + 1) * c(w, (nix +
                                                                 1), (nix.ridge + 1)), intercept = F)
      }
      if(ridge.adj == 0) {
        f <- lsfit(rbind(b, p), c(u, nix), wt = c(wts, nix + 1) *
                     c(w, (nix + 1)), intercept = F)
      }
      coef.old <- coef.est
      coef.est <- as.vector(f$coef)
      d.coef <- max(abs((coef.est - coef.old)/coef.old))
      if(d.coef < 1e-008)
        break
      print(c(it, d.coef))
      eta <- b %*% coef.est
    }
    if(it > 49) {
      warning(paste("parameter estimates did NOT converge in 50 iterations"
      ))
    }
    llist <- list(coef = coef.est, mu = mu, f = f, w = w * wts, eta=eta)
    return(llist)
  }

#' @export
"pspline.checker"<-
  function(family, link, bdeg, pord, nseg, lambda, ridge.adj, wts)
  {
    if(link == "default" && family == "gaussian") {
      link <- "identity"
    }
    if(link == "default" && family == "poisson") {
      link <- "log"
    }
    if(link == "default" && family == "binomial") {
      link <- "logit"
    }
    if(link == "default" && family == "Gamma") {
      link <- "log"
    }
    if(family != "binomial" && family != "gaussian" && family != "poisson" &&
       family != "Gamma") {
      warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"
      ))
    }
    if((family == "binomial") && (link != "logit" && link != "probit" &&
                                  link != "cloglog" && link != "loglog")) {
      warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"
      ))
    }
    if((family == "Gamma") && (link != "log" && link != "recipical" && link !=
                               "identity")) {
      warning(paste("Improper LINK option with family=Gamma. Choose: recipical, log, identity"
      ))
    }
    if((family == "poisson") && (link != "log" && link != "sqrt" && link !=
                                 "identity")) {
      warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"
      ))
    }
    if((family == "gaussian") && (link != "identity")) {
      warning(paste("Improper LINK option with family=gaussian. Choose: identity"
      ))
    }
    if(bdeg < 0) {
      bdeg <- 1
      warning(paste("bdeg must be non-neg integer: have used 1"))
    }
    if(pord < 0) {
      pord <- 0
      warning(paste("pord must be non-neg integer: have used 0"))
    }
    if(nseg < 2) {
      nseg <- 2
      warning(paste("nseg must be positive integer, > 1: have used 2"
      ))
    }
    if(lambda < 0) {
      lambda <- 0
      warning(paste("lambda cannot be negative: have used 0"))
    }
    if(ridge.adj < 0) {
      ridge.adj <- 0
      warning(paste("ridge.adj cannot be negative: have used 0"))
    }
    if(min(wts) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    llist <- list(family = family, link = link, bdeg = bdeg, pord =
                    pord, nseg=nseg, lambda = lambda, ridge.adj
                  = ridge.adj, wts = wts)
    return(llist)
  }

#' @export
"pspline2d.checker"<-
  function(family, link, degree1, degree2, order1, order2, ps.intervals1,
           ps.intervals2, lambda1, lambda2, ridge.adj, wts)
  {
    if(link == "default" && family == "gaussian") {
      link <- "identity"
    }
    if(link == "default" && family == "poisson") {
      link <- "log"
    }
    if(link == "default" && family == "binomial") {
      link <- "logit"
    }
    if(link == "default" && family == "Gamma") {
      link <- "log"
    }
    if(family != "binomial" && family != "gaussian" && family != "poisson" &&
       family != "Gamma") {
      warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"
      ))
    }
    if((family == "binomial") && (link != "logit" && link != "probit" &&
                                  link != "cloglog" && link != "loglog")) {
      warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"
      ))
    }
    if((family == "Gamma") && (link != "log" && link != "recipical" && link !=
                               "identity")) {
      warning(paste("Improper LINK option with family=Gamma. Choose: recipical, log, identity"
      ))
    }
    if((family == "poisson") && (link != "log" && link != "sqrt" && link !=
                                 "identity")) {
      warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"
      ))
    }
    if((family == "gaussian") && (link != "identity")) {
      warning(paste("Improper LINK option with family=gaussian. Choose: identity"
      ))
    }
    if(degree1 < 0) {
      degree1 <- 1
      warning(paste("degree1 must be non-neg integer: have used 1"))
    }
    if(order1 < 0) {
      order1 <- 0
      warning(paste("order1 must be non-neg integer: have used 0"))
    }
    if(ps.intervals1 < 2) {
      ps.intervals1 <- 2
      warning(paste("ps.intervals1 must be positive integer, > 1: have used 2"
      ))
    }
    if(lambda1 < 0) {
      lambda1 <- 0
      warning(paste("lambda1 cannot be negative: have used 0"))
    }
    if(degree2 < 0) {
      degree2 <- 1
      warning(paste("degree2 must be non-neg integer: have used 1"))
    }
    if(order2 < 0) {
      order2 <- 0
      warning(paste("order2 must be non-neg integer: have used 0"))
    }
    if(ps.intervals2 < 2) {
      ps.intervals2 <- 2
      warning(paste("ps.intervals2 must be positive integer, > 1: have used 2"
      ))
    }
    if(lambda2 < 0) {
      lambda2 <- 0
      warning(paste("lambda2 cannot be negative: have used 0"))
    }
    if(ridge.adj < 0) {
      ridge.adj <- 0
      warning(paste("ridge.adj cannot be negative: have used 0"))
    }
    if(min(wts) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    llist <- list(family = family, link = link, degree1 = degree1, order1
                  = order1, ps.intervals1 = ps.intervals1, lambda1 = lambda1,
                  degree2 = degree2, order2 = order2, ps.intervals2 =
                    ps.intervals2, lambda2 = lambda2, ridge.adj = ridge.adj, wts =
                    wts)
    return(llist)
  }
