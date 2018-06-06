#' @export print.ata

print.ata <- function(x,...)
{
if (class(x)=="ata.forecast"){
	cat("Forecasts:","\n")
	print(x$forecast)
	cat("\n")
}else {
	cat(x$method,"\n\n")
	if (x$level.fixed==TRUE){
		cat("   level.fixed: TRUE","\n\n")
	}
	if (x$trend.fixed==TRUE){
		cat("   trend.fixed: TRUE","\n\n")
	}
	if(!is.null(x$transform.method)){
		cat(paste("   '",x$transform.method, "' transformation method was selected.","\n\n", sep=""))
	}
	if(!is.null(x$lambda)){
		cat("   Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
	}
	cat(paste("   model.type:",x$model.type, "\n\n"))
	if (x$is.season==FALSE){
		cat("   seasonal.model: no seasonality","\n\n")
	}else {
		cat(paste("   seasonal.model:",x$seasonal.model, "\n\n"))
	}
	if (x$is.season==TRUE){
		cat(paste("   seasonal.type:",x$seasonal.type, "\n\n"))	
	}
	cat(paste("   forecast horizon:",x$h, "\n\n"))
	cat(paste("   accuracy.type:",x$accuracy.type, "\n\n"))
	
	cat("In-Sample Accuracy Measures:","\n")
	stats <- c(x$accuracy$MAE$inSample$MAE, x$accuracy$MAE$inSample$MdAE, x$accuracy$MAE$inSample$stdDev.MAE, x$accuracy$MAE$inSample$skewness.MAE, x$accuracy$MAE$inSample$kurtosis.MAE)
    names(stats) <- c("MAE","MdAE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$stdDev.MSE, x$accuracy$MSE$inSample$skewness.MSE, x$accuracy$MSE$inSample$kurtosis.MSE)
    names(stats) <- c("MSE", "RMSE", "MdSE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$MPE$inSample$MPE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MPE$inSample$stdDev.MPE, x$accuracy$MPE$inSample$skewness.MPE, x$accuracy$MPE$inSample$kurtosis.MPE)
    names(stats) <- c("MPE","MdPE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")

	stats <- c(x$accuracy$MAPE$inSample$MAPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$MAPE$inSample$stdDev.MAPE, x$accuracy$MAPE$inSample$skewness.MAPE, x$accuracy$MAPE$inSample$kurtosis.MAPE)
    names(stats) <- c("MAPE","MdAPE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$sMAPE$inSample$sMdAPE, x$accuracy$sMAPE$inSample$stdDev.sMAPE, x$accuracy$sMAPE$inSample$skewness.sMAPE, x$accuracy$sMAPE$inSample$kurtosis.sMAPE)
    names(stats) <- c("sMAPE","sMdAPE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")

	stats <- c(x$accuracy$MASE$inSample$MASE, x$accuracy$MASE$inSample$MdASE, x$accuracy$MASE$inSample$stdDev.MASE, x$accuracy$MASE$inSample$skewness.MASE, x$accuracy$MASE$inSample$kurtosis.MASE)
    names(stats) <- c("MASE","MdASE","StdDev", "Skewness", "Kurtosis")
    cat("\n")
    print(stats)
	cat("\n")
	
	cat("Out-Sample Accuracy Measures:","\n\n")
	stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MAE$outSample$MdAE)
    names(stats) <- c("MAE","MdAE")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MSE$outSample$MdSE)
    names(stats) <- c("MSE", "RMSE", "MdSE")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$MPE$outSample$MPE, x$accuracy$MPE$outSample$MdPE)
    names(stats) <- c("MPE","MdPE")
    cat("\n")
    print(stats)
	cat("\n")

	stats <- c(x$accuracy$MAPE$outSample$MAPE, x$accuracy$MAPE$outSample$MdAPE)
    names(stats) <- c("MAPE","MdAPE")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("sMAPE","sMdAPE")
    cat("\n")
    print(stats)
	cat("\n")

	stats <- c(x$accuracy$MASE$outSample$MASE, x$accuracy$MASE$outSample$MdASE)
    names(stats) <- c("MASE","MdASE")
    cat("\n")
    print(stats)
	cat("\n")
	
	stats <- c(x$execution.time[1], x$execution.time[2], x$execution.time[3])
    names(stats) <- c("user","system","elapsed")
    cat("\n")
    print(stats)
	cat("\n")
	cat(paste("calculation.time:",x$calculation.time, "\n\n"))
}
}