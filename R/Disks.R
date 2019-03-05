#' Prices of hard disk drives
#'
#' Prices and capacities of hard disk drives,
#' as advertised in a Dutch computer monthly in 1999.
#' Prices are given in Dutch guilders; the Euro did not yet exist.
#'
#' @docType data
#'
#' @usage data(Disks)
#'
#' @format A dataframe with six columns:
#' \describe{
#'   \item{\code{Year}}{1999-2000}
#'   \item{\code{Month13}}{month of data}
#'   \item{\code{Size}}{capacity in Gb}
#'   \item{\code{Buffer}}{buffer size (MB)}
#'   \item{\code{RPM}}{rotating speed}
#'   \item{\code{PriceDG}}{in Dutch Guilders, divide by 2.2 for Euro}}
#'
#' @keywords datasets
#'
#' @references
#'
#' @source
#'
#' @examples
#' plot(Size[Month==5], Price[Month==5], type='b',
#' xlab='Size', ylab='Price', main='Month of May 1999')

"Disks"
