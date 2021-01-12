#' Find Frequency Using Spectral Density Of A Time Series From AR Fit
#'
#' @param x an univariate time series
#'
#' @return frequency/cycle of the given time data
#' 
#' @importFrom stats spec.ar
#' 
#' @export
#'
find.freq <- function(x)
{
  n <- length(x)
  spec <- stats::spec.ar(c(x),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
}
