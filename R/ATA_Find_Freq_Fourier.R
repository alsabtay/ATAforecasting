#' Find Frequency Using Periodogram
#'
#' @param x an univariate time series
#'
#' @return frequency/cycle of the given data
#' 
#' @importFrom TSA periodogram
#' @importFrom utils head
#' 
#' @export
#'
find.freq.fourier <- function(x)
{
  pppx <- TSA::periodogram(x)
  dddx = data.frame(freq=pppx$freq, spec=pppx$spec)
  orderpppx = dddx[order(-dddx$spec),]
  top5X = utils::head(orderpppx, 5)
  freq_all <- 1/top5X$freq
  period <- sort(freq_all)
  period <- period[period < 367]
  return(period)
}
