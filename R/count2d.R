#' Provides the counts in 2D histograms.
#'
#' @description A support function for \code{hist2d()}
#' that gives the counts in 2D bin.
#'
#'
#' @param xb a vector along \code{x} to be binned in 2D.
#' @param yb a vector along \code{y} to be binned in 2D.
#' @param nb a vector of length 2 that provides the number of bins for the 2D histogram on \code{x} and \code{y}.
#' @return a matrix of dimension \code{nb[1]} by \code{nb[2]} with counts.
#'
#'
#' @examples
#' count2d(1:4, sample(1:4), c(4,4))
#' @export


count2d <- function(xb, yb, nb) {
  # Fill nb[1] by nb[2] histogram; indexed by xb and yb
  H <- matrix(rep(0, prod(nb)), nb[1], nb[2])
  for (i in 1:length(xb)) H[xb[i], yb[i]] = H[xb[i], yb[i]] + 1
  H}
