#' Sugar Processing Data
#'
#' Sugar was sampled continuously during eight hours to make
#' a mean sample representative for one 'shift'
#' (eight hour period). Samples were taken during the
#' three months of operation (the so-called campaign) in
#' late autumn from a sugar plant in Scandinavia giving a
#' total of 268 samples. The sugar was sampled directly
#' from the final unit operation (centrifuge) of the process.
#'
#' @docType data
#'
#' @usage data(Sugar)
#'
#' @format A list consisting of the following:
#' \describe{
#' \item{\code{y}}{a 268 x 3 matrix of quality parameters: \code{date}, \code{color}, \code{ash}*1000}
#' \item{\code{X}}{fluoresence array, 268 (observations) x [571 (emission channels) x 7 (excitation channels)]}
#' \item{Lab}{Lab information}
#' \item{DimX}{array dimension for \code{X}}
#' \item{Yidx}{names (id) for \code{y}}
#' \item{time}{.}
#' \item{EmAx}{.}
#' \item{ExAx}{.}
#' \item{readmetime}{.}
#' \item{Lname}{.}
#' \item{LabNumber}{.}
#' \item{ProcNumber}{.}
#' \item{Proc}{.}
#' \item{DimLab}{.}
#' \item{DimProc}{.}}
#'
#' @keywords datasets
#'
#' @references R. Bro, Exploratory study of sugar production using fluorescence spectroscopy and multi-way analysis,
#' Chemom. Intell. Lab. Syst., 1999, (46), 133-147.

#' @source http://www.models.life.ku.dk/sugar_process
#'
#'
"Sugar"
