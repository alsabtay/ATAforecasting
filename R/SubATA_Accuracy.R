#' @importFrom stats frequency median
SubATA.Accuracy <- function(ata_opt, accryType){
  inSample <- ata_opt$actual
  in_sample_fit <- ata_opt$fitted
  in_sample <- ATA.Shift(as.numeric(inSample),1)
  ata.error <- in_sample - in_sample_fit
  ata.pe <- ata.error / in_sample * 100
  if (accryType=="MAE" | accryType=="MdAE"){
    ata.accuracy.insample <- abs(ata.error)
  }else if (accryType=="MSE" | accryType=="MdSE" | accryType=="RMSE" | accryType=="RMdSE"){
    ata.accuracy.insample <- ata.error^2
  }else if (accryType=="MPE" | accryType=="MdPE"){
    ata.accuracy.insample <- ata.pe
  }else if (accryType=="MAPE" | accryType=="MdAPE"){
    ata.accuracy.insample <- abs(ata.pe)
  }else if (accryType=="sMAPE" | accryType=="sMdAPE"){
    ata.accuracy.insample <- abs(in_sample - in_sample_fit)/(abs(in_sample) + abs(in_sample_fit)) * 200
  }else if (accryType=="MASE"){
    ata.accuracy.insample <- inMASE(as.double(in_sample), as.double(in_sample_fit), as.integer(frequency(inSample)))
  }else if (accryType=="OWA"){
    preOWA_first <- abs(in_sample - in_sample_fit)/(abs(in_sample) + abs(in_sample_fit)) * 200
    preOWA_second <- abs(inMASE(as.double(in_sample), as.double(in_sample_fit), as.integer(frequency(inSample))))
    naiveAccry <- round(NaiveSD(as.double(in_sample), as.integer(frequency(inSample))),6)
  }else {
  }

  if (accryType=="MAE" | accryType=="MSE" | accryType=="MPE" | accryType=="MAPE" | accryType=="sMAPE"){
    m_accry <- mean(ata.accuracy.insample, na.rm=TRUE)
  }else if (accryType=="MdAE" | accryType=="MdSE" | accryType=="MdPE" | accryType=="MdAPE" | accryType=="sMdAPE"){
    m_accry <- median(ata.accuracy.insample, na.rm=TRUE)
  }else if (accryType=="MASE"){
    m_accry <- ata.accuracy.insample
  }else if (accryType=="RMSE") {
    m_accry <- sqrt(mean(ata.accuracy.insample, na.rm=TRUE))
  }else if (accryType=="RMdSE") {
    m_accry <- sqrt(median(ata.accuracy.insample, na.rm=TRUE))
  }else if (accryType=="OWA") {
    OWA_first<- mean(preOWA_first, na.rm=TRUE)
    OWA_second <- preOWA_second
    m_accry <- ((OWA_first/naiveAccry) + (OWA_second/naiveAccry)) / 2
  }else{
  }
  return(m_accry)
}
