#' @export
find.freq.fourier <- function(x)
{
	pppx <- periodogram(x)
	dddx = data.frame(freq=pppx$freq, spec=pppx$spec)
	orderpppx = dddx[order(-dddx$spec),]
	top5X = head(orderpppx, 5)
	freq_all <- 1/top5X$freq
	period <- sort(freq_all)
	period <- period[period < 367]
    return(period)
}


#' @export
find.freq <- function(x)
{
    n <- length(x)
    spec <- spec.ar(c(x),plot=FALSE)
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


#' @export
find.freq.all <- function(x)
{  
  f=find.freq(x)
  if (is.na(f)){
	f <- 1
  }
  freqs=c(f) 
  while(f>1){
    start=1 #also try start=f;
    x=period.apply(x,seq(start,length(x),f),mean) 
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


