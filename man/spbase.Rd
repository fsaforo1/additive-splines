% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spbase.R
\name{spbase}
\alias{spbase}
\title{Sarse B-spline basis on evenly spaced knots}
\usage{
spbase(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3)
}
\arguments{
\item{x}{a vector of argument values, at which the B-spline basis functions
are to be evaluated.}

\item{xl}{the lower limit of the domain of x.}

\item{xr}{the upper limit of the domain of x.}

\item{nseg}{the number of evenly spaced segments between xl and xr.}

\item{bdeg}{the degree of the basis, usually 1, 2, or 3.}
}
\value{
a sparse matrix with number of rows=\code{length(x)} and
number of columns = \code{nseg + bdeg}.
}
\description{
Constructs a B-spline basis on evenly spaced knots, used
for P-splines. A sparse matrix, in the format of the packag \code{spam}, is returned.
}
\examples{
library(JOPS)
# Basis  on grid
xg = seq(0, 4, length = 500)
Bg = spbase(xg, 0, 4, nseg = 3, bdeg = 3)
nb1 = ncol(Bg)
matplot(xg, Bg, type = 'l', lty = c(1:6), lwd = 2, xlab= 'x', ylab='')
}
\references{
Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder), Statistical Science,
11: 89-121.
}
\author{
Paul Eilers and Brian Marx
}
