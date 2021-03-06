% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psBinomial.R
\name{psBinomial}
\alias{psBinomial}
\title{Smoothing scattered binomial data using P-splines.}
\usage{
psBinomial(x, y, xl = min(x), xr = max(x), nseg = 10, bdeg = 3,
  pord = 2, lambda = 1, wts = NULL, show = F)
}
\arguments{
\item{x}{the vector for the continuous regressor of \code{length(y)} and
the abcissae of fit.}

\item{y}{the response vector, usually 0/1 data.}

\item{xl}{the number for the min along \code{x}.}

\item{xr}{the number for the max along \code{x}.}

\item{nseg}{the number of evenly spaced segments between xl and xr.}

\item{bdeg}{the number of the degree of the basis, usually 1, 2, or 3 (defalult).}

\item{pord}{the number of the order of the difference penalty, usually 1, 2, or 3 (defalult).}

\item{lambda}{the (positive) number for the tuning parameter for the penalty.}

\item{wts}{the vector of weights, default is 1, 0/1 allowed.}

\item{show}{Set to T or F to display iteration history.}
}
\value{
\item{pcoeff}{a vector of length \code{n} of estimated P-spline coefficients.}
\item{muhat}{a vector of length \code{m} of estimated means (probabilities).}
\item{dev}{deviance}
\item{effdim}{effective dimension of the smooth}
\item{aic}{AIC}
\item{wts}{a vector of preset weights.}
\item{nseg, bdeg, pord, lambda}{design parameters.}
\item{xgrid}{gridded x values for plotting.}
\item{ygrid}{gridded fitted probability values for plotting.}
\item{lgrid}{gridded lower 2se values for plotting.}
\item{ugrid}{gridded lower 2se values for plotting.}
}
\description{
psBinomial is used to smooth scattered
binomial data using P-splines usign a logit link function.
}
\examples{
library(JOPS)
#Extract data
library(rpart)
Kyphosis=kyphosis$Kyphosis
Age=kyphosis$Age
y = 1 * (Kyphosis == "present")  # make y 0/1
fit1 = psBinomial(Age, y,xl=min(Age), xr=max(Age), nseg=20,
                 bdeg=3, pord=2, lambda= 1)
names(fit1)
plot(fit1, xlab="Age", ylab='0/1', se=2)
}
\references{
Eilers, P.H.C., Marx, B.D., and Durban, M. (2015).
Twenty years of P-splines, SORT, 39(2): 149-186.
}
\author{
Paul Eilers and Brian Marx
}
