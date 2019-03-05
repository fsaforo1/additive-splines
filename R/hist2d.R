#' Compile a 2D histogram.
#'
#' @param x a numeric vector \code{x}.
#' @param y a numeric vector \code{y}.
#' @param xb a vector of length 2, providing the number of bins along x, y, resp. .
#' @param xlim default is min; each is a vector of length 2.
#' @param ylim default is max; each is a vector of length 2.

#' @return A list with components
#' \item{H}{a matrix of dimension \code{nb[1]} by \code{nb[2]} containing bin counts.}
#' \item{xgrid}{a vector of \code{length(nb[1])} representing the interior bin limits.}
#' \item{ygrid}{a vector of \code{length(nb[2])} representing the interior bin limits.}
#' \item{xbin}{a vector of \code{length(x)} with elements of the location (bin number) for each x.}
#' \item{ybin}{a vector of \code{length(y)} with elements of the location (bin number) for each x.}
#'
#' @examples
#' #' data(faithful)
#' x=faithful$eruptions
#' y=faithful$waiting
#' C=hist2d(x,y, c(50,50))
#' image(C$xgrid,C$ygrid, C$H, xlab='Eruptions', ylab='Wating')
#'
#' @export


# Support functions for 2-D histogram smoothing

hist2d = function(x, y, nb = c(100, 100),
                  xlim = c(min(x), max(x)), ylim = c(min(y), max(y))) {
  # Compute a 2-D histogram
  xb = binit(x, xmin = xlim[1], xmax = xlim[2], nbin = nb[1])
  yb = binit(y, xmin = ylim[1], xmax = ylim[2], nbin = nb[2])
  xx = xb$xbin
  yy = yb$xbin
  sel = 0 < xx & xx <= nb[1] & 0 < yy & yy <= nb[2]
  H = count2d(xx[sel], yy[sel], nb)
  h2d = list(H = H, xgrid = xb$xgrid, ygrid = yb$xgrid, xbin = xx, ybin = yy)
  h2d
}
