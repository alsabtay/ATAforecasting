#' @export

AutoATA.Accuracy <- function(ata_opt, accryType){
	inSample <- ata_opt$actual
	in_sample_fit <- ata_opt$fitted
	in_sample <- ATA.Shift(as.numeric(inSample),1)
	ata.error <- in_sample - in_sample_fit
	ata.pe <- ata.error / in_sample * 100
	if (accryType=="MAE" | accryType=="MdAE" | accryType=="MASE"){
		ata.accuracy.insample <- abs(ata.error)
	}else if (accryType=="MSE" | accryType=="MdSE" | accryType=="RMSE"){
		ata.accuracy.insample <- ata.error^2
	}else if (accryType=="MPE" | accryType=="MdPE"){
		ata.accuracy.insample <- ata.pe
	}else if (accryType=="MAPE" | accryType=="MdAPE"){
		ata.accuracy.insample <- abs(ata.pe)
	}else if (accryType=="sMAPE" | accryType=="sMdAPE"){
		ata.accuracy.insample <- abs(in_sample - in_sample_fit)/(abs(in_sample) + abs(in_sample_fit)) * 200
	}else {
	}
	
	if (accryType=="MAE" | accryType=="MSE" | accryType=="MPE" | accryType=="MAPE" | accryType=="sMAPE"){
		m_accry <- mean(ata.accuracy.insample, na.rm=TRUE)
	}else if (accryType=="MdAE" | accryType=="MdSE" | accryType=="MdPE" | accryType=="MdAPE" | accryType=="sMdAPE"){
		m_accry <- median(ata.accuracy.insample, na.rm=TRUE)
	}else if (accryType=="MASE"){
		naive_benchmark <- mean(abs(inSample-mean(inSample, na.rm=TRUE)),na.rm=TRUE)
		m_accry <- mean(abs(insample.error/naive_benchmark), na.rm=TRUE)
	}else if (accryType=="RMSE") {
		m_accry <- sqrt(mean(ata.accuracy.insample, na.rm=TRUE))
	}else{
	}	
    return(m_accry)
}