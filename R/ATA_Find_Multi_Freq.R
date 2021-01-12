#' Find Multi Frequency Using Spectral Density Of A Time Series From AR Fit
#'
#' @param x an univariate time series
#'
#' @return multi frequencies/cycles of the given data
#' 
#' @importFrom xts period.apply
#' 
#' @export
#'
find.multi.freq <- function(x)
{
	f = find.freq(x)
	if (is.na(f)){
		f <- 1
	}
	freqs=c(f)
	while(f>1){
		start=1 #also try start=f;
		x=xts::period.apply(x,seq(start,length(x),f),mean)
		f=find.freq(x)
		freqs=c(freqs,f)
		if (is.na(f)){
		  f <- 1
		}
	}
	if(length(freqs)==1){ return(freqs) }
	for(i in 2:length(freqs)){
		freqs[i]=freqs[i]*freqs[i-1]
	}
	freqs[1:(length(freqs)-1)]
}
