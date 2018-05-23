#' @export ATA.CI

ATA.CI <- function(ata.output, ci.level = 95)
{
	if (class(ata.output)!="ata"){
		return("The Input must be 'ata' object. Please use ATA function to produce 'ata' object. Calculation of Confidence Intervals of ATA Forecasts will terminate!")
	}
	ci.alpha <- 1 - (ci.level/100)
	length_resid <- length(ata.output$residual[!is.na(ata.output$residual)])
	if (length_resid<=30){
		ci.ZTvalue <- qt(ci.alpha/2, df=length_resid)
	}else {
		ci.ZTvalue <- qnorm(ci.alpha/2, lower.tail=FALSE)
	}
	std_resid <- sd(ata.output$residual, na.rm=TRUE)
	ci.value <- sqrt(1:ata.output$h) * ci.ZTvalue * std_resid
	tsp_F <- ata.output$forecast 
	forecast.lower <- ts(ata.output$forecast - ci.value, f = tsp_F[3], s = tsp_F[1]) 
	forecast.upper <- ts(ata.output$forecast + ci.value, f = tsp_F[3], s = tsp_F[1])
	my_list <- list("forecast"=ata.output$forecast, "forecast.lower"=forecast.lower, "forecast.upper"=forecast.upper)
	return(my_list) 
}