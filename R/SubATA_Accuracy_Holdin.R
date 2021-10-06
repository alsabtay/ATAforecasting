#' @importFrom stats frequency median
#' @importFrom utils tail
SubATA.Accuracy.Holdin <- function(ata_opt, accryType, h){
  inSample <- ata_opt$actual
  HoldinSet = utils::tail(inSample, h);
  ata.error <- HoldinSet - ata_opt$forecast
  ata.pe <- ata.error / HoldinSet * 100
  if (accryType=="MAE" | accryType=="MdAE"){
    ata.accuracy.insample <- abs(ata.error)
  }else if (accryType=="MSE" | accryType=="MdSE" | accryType=="RMSE"){
    ata.accuracy.insample <- ata.error^2
  }else if (accryType=="MPE" | accryType=="MdPE"){
    ata.accuracy.insample <- ata.pe
  }else if (accryType=="MAPE" | accryType=="MdAPE"){
    ata.accuracy.insample <- abs(ata.pe)
  }else if (accryType=="sMAPE" | accryType=="sMdAPE"){
    ata.accuracy.insample <- abs(ata.error)/(abs(HoldinSet) + abs(ata_opt$forecast)) * 200
  }else if (accryType=="MASE"){
    ata.accuracy.insample <- outMASE(as.double(inSample), HoldinSet, as.double(ata_opt$forecast), as.integer(frequency(inSample)))
  }else if (accryType=="OWA"){
    preOWA_first <- abs(ata.error)/(abs(HoldinSet) + abs(ata_opt$forecast)) * 200
    preOWA_second <- abs(outMASE(as.double(inSample), HoldinSet, as.double(ata_opt$forecast), as.integer(frequency(inSample))))
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
  }else if (accryType=="OWA") {
    OWA_first<- mean(preOWA_first, na.rm=TRUE)
    OWA_second <- preOWA_second
    m_accry <- (OWA_first + OWA_second) / 2
  }else{
  }
  return(m_accry)
}
