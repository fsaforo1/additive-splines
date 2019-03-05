#' Environmental complaints in The Netherlands
#'
#' The Rijnomond Environmental Agency registers approximately 20,000
#' complaints annually from regional inhabitants regarding odors, dust,
#' noise, among other annoyances. The data are daily number of complaints (
#' Complaints$freq), ranging from 0-100. The frequency of the
#' each number of complaints is tallied (Complaints$count). The data
#' are from 1988.
#'
#' @docType data
#'
#' @usage data(Complaints)
#'
#'
#' @format A dataframe with two columns (\code{freq} and \code{count}):
#' \describe{
#'   \item{\code{freq}}{The daily number of complaints.}
#'   \item{\code{count}}{The number of days the specific complaint frequency occurred.}
#'   }
#'
#'
#'
#' @keywords datasets
#'
#' @references
#'
#' @source
#'
#' @examples
#' plot(Complaints$freq, Complaints$count, type='h',
#' xlab='Number of complaints per day', ylab='Frequency')

"Complaints"
