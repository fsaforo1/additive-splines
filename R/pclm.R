#' Fittting composite link functions using P-splines
#'
#' @description pclm is used to fit penalized composite link functions.
#' It can be useful for density estimation with coarsely aggregated data
#' or digit prefence data.
#'
#'
#' @param y a vector, usually of counts to be fit.
#' @param C a composition matrix with number of rows equal to \code{length(y)} and
#' number columns equal to the latent distribution, \code{q=length(gamma)}.
#' @param X a matrix, ususally a B-spline basis built along the centers of
#' the bins, with the number of rows \code{q=length(gamma)}. The basis should have enough
#' columns to have sufficient flexibilty with \code{ncol(X)=n}
#'  A number that is the upper limit of \code{x}.
#' @param lambda the number for the (positive) tuning parameter of the penalty.
#' @param pord the number for the order of the difference penalty (default=2) usually 1, 2, or 3.
#' @param itmax  the maximum number of iterations (defualt=50).
#' @param show Set to T or F to display iteration history.
#'
#' @return
#' \item{coeff}{a vector of length \code{n} of estimated P-spline coefficients.}
#' \item{gamma}{a vector of length \code{q} for the latent coefficients.}
#' \item{mu}{a vector of length \code{m} of estimated means.}
#' \item{dev}{Deviance of fit}
#' \item{ed}{effective dimension}
#' \item{aic}{AIC}.
#' @author Paul Eilers and Jutta Gampe
#' @references
#' Eilers, P. H. C. (2007). III-posed problems with counts, the composite link
#' model and penalized likelihood. Statistical Modelling, 7(3), 239–254.
#' @examples
#' library(ggplot2)
#' library(JOPS)
#' cb <- c( 0, 20, 30, 40, 50, 60)
#' ce <- c(20, 30, 40, 50, 60, 70)
#' y <- c(79, 54, 19, 1, 1, 0)
#' m <- length(y)
#' n <- ce[m]
#' C <- matrix(0, m, n)
#' for (i in 1:m) C[i, cb[i]:ce[i]] <- 1
#' mids = (cb + ce) / 2 - 0.5
#' widths = ce - cb + 1
#' dens = y / widths / sum(y)
#' x=1:n
#' fit = pclm(y, C, X=bbase(x), lambda = 2, pord = 2, show = T)
#' gam2 = fit$gamma

#' #Plot with ggplot2
#' Data = data.frame(x = mids, y = dens)
#' Fit = data.frame(x = x - 0.5, y = gam2 / sum(gam2))
#' Dat = data.frame(cb, ce, dens)
#' plt = ggplot(Dat, aes(ymin = 0)) +
#'  geom_rect(aes(xmin = cb, xmax = ce, ymax = dens, fill = I("grey"),color=I('white'))) +
#'  geom_line(aes(x = x, y = y), size = 1, linetype = 1,data = Fit, col = I("blue")) +
#'  xlab("Lead concentration") + ylab("Density") +
#'  ggtitle("Lead in blood") +
#'  JOPS_theme()
#' set_window()
#' print(plt)
#'
#' @export


pclm <- function(y, C, X, lambda = 1, pord = 2, itmax = 50,  show = F) {
  # Joint work of Jutta Gampe and Paul Eilers, 2016

	# Set up penalty
	n <- dim(X)[2]
	D <- diff(diag(n), diff = pord)
	P <- lambda * t(D) %*% D

	# Start for b (uniform)
	b <- log(rep(sum(y) / n, n))

	# Iterate
	for (it in 1:50) {
		b.old  <- b
		eta <- X %*% b
		gamma <- exp(eta)
		mu  <- as.vector(C %*% gamma)
		M <- diag(1 / mu)
		U <- M %*% C %*% diag(c(gamma)) %*% X
		Q <- t(U) %*% diag(mu) %*% U
		z <- t(U) %*% (y - mu) + Q %*% b
		b <- solve(Q + P, z)
		diff <- max(abs(b - b.old))
		if (show)  cat(it, "  ", max(abs(b - b.old)), "\n")
		if (diff < 1e-7) break
	}

	if (it >= itmax) cat("No convergence after", itmax, "iterations\n")

	# Diagnostics
	H <- solve(Q + P) %*% Q
	ok <- y > 0
	dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
	ed  <- sum(diag(H))
	aic <- dev + 2 * ed

	# Prepare output
	fit = list(coeff = b, gamma = gamma, mu = mu, dev = dev, ed = ed, aic = aic)
	return(fit)
}

