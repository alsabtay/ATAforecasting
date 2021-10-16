#' Confidence Interval function for the ATA Method
#'
#' @param object An \code{ATA} object is required.
#' @param ci.level Confidence level, for example: 90, 95 or 99.
#'
#' @return The confidence interval output for the ATA forecasts
#'
#' @importFrom stats qnorm qt sd ts tsp tsp<-
#'
#' @export
#'
ATA.CI <- function(object, ci.level = 95)
{
	ata.output <- object
		if (class(ata.output)!="ata"){
		return("The Input must be 'ata' object. Please use ATA function to produce 'ata' object. Calculation of Confidence Intervals of ATA Forecasts will terminate!")
	}
	ci.alpha <- 1 - (ci.level/100)
	length_resid <- length(ata.output$residuals[!is.na(ata.output$residuals)])
	if (length_resid<=30){
		ci.ZTvalue <- stats::qt(ci.alpha/2, df=length_resid)
	}else {
		ci.ZTvalue <- stats::qnorm(ci.alpha/2, lower.tail=FALSE)
	}
	std_resid <- sd(ata.output$residual, na.rm=TRUE)
	ci.value <- sqrt(1:ata.output$h) * abs(ci.ZTvalue) * std_resid
	tsp_F <- tsp(ata.output$forecast)
	forecast.lower <- ts(ata.output$forecast - ci.value, frequency = tsp_F[3], start = tsp_F[1])
	forecast.upper <- ts(ata.output$forecast + ci.value, frequency = tsp_F[3], start = tsp_F[1])
	my_list <- list("forecast"=ata.output$forecast, "forecast.lower"=forecast.lower, "forecast.upper"=forecast.upper)
	return(my_list)
}
