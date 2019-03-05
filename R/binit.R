#' Translated number vector to bin index.
#'
#' @description Translates number vector to bin index, given lower and
#' upper limits of the domain and number of bins. A support function
#' for histograms.
#'
#' @param x a numerical vector.
#' @param xl the lower limit of the domain.
#' @param xr the upper limit of the domain.
#' @param nbin the number of bins.
#' @return A list with components
#' \item{xbin}{a vector of \code{length(nbin)} with elements of the location (bin number) for each x.}
#' \item{xgrid}{a vector of length(x) representing the interior bin limits.}
#' \item{nbin}{the number of bins.}
#'
#'
#' @examples
#' binit(runif(20), 0, 1, 5)
#' @export


# Support functions for histogram smoothing

binit = function(x, xmin = min(x), xmax = max(x), nbin = 100) {
  # Put the x-values into bins
  dx <- (xmax - xmin) / (nbin - 1)
  xbin <- floor(1 + (x - xmin) / dx)
  xgrid <- xmin + (1:nbin - 0.5) * dx
  b = list(xbin = xbin, xgrid = xgrid, nbin = nbin)
  b }
