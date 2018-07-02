#' @export SeasonalityTest

SeasonalityTest <- function(input, ppy, attr_set)
{
	if (ppy==1){
		test_seasonal <- FALSE
	}else {
		s.tcrit <- attr_set$s.tcrit
		uroot.test <- attr_set$uroot.test
		uroot.type <- attr_set$uroot.type
		uroot.alpha <- attr_set$uroot.alpha
		uroot.maxd <- attr_set$uroot.maxd
		if (length(ppy)>1){
			ppy <- min(ppy)
		}
	  #Used to determine whether a time series is seasonal

		d <- forecast::ndiffs(input, alpha=uroot.alpha, test=uroot.test, type=uroot.type, max.d=uroot.maxd)
		if(d > 0) {
			input <- diff(input, differences=d, lag=1)
		}
		if (length(input)<3*ppy){
			test_seasonal <- FALSE
		}else {
			xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
			clim <- s.tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
			test_seasonal <- (abs(xacf[ppy]) > clim[ppy])
			if (is.na(test_seasonal)==TRUE){
				test_seasonal <- FALSE
			}
		}
	}
	return(test_seasonal)
}

#' @export find.freq.fourier

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


#' @export find.freq

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

#' @export find.freq.all

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