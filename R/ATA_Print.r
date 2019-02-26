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
	stats <- c(x$accuracy$MAE$inSample$MAE, x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MPE$inSample$MPE, x$accuracy$MAPE$inSample$MAPE, x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$MASE$inSample$MASE, x$accuracy$OWA$inSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(stats)
	cat("\n")
	
	cat("In-Sample Accuracy Measures:","\n")
	stats <- c(x$accuracy$MAE$inSample$MdAE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$RMdSE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$sMAPE$inSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
	cat("\n")
	
	
	cat("Out-Sample Accuracy Measures:","\n\n")
	stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MPE$outSample$MPE, x$accuracy$MAPE$outSample$MAPE, x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$MASE$outSample$MASE, x$accuracy$OWA$outSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(stats)
	cat("\n")

	cat("Out-Sample Accuracy Measures:","\n\n")	
	stats <- c(x$accuracy$MAE$outSample$MdAE, x$accuracy$MSE$outSample$MdSE, x$accuracy$MSE$outSample$RMdSE, x$accuracy$MPE$outSample$MdPE, x$accuracy$MAPE$outSample$MdAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
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