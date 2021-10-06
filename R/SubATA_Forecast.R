#' @importFrom stats ts tsp tsp<-
SubATA.Forecast <- function(ata_output, hh=NULL, initialLevel)
{
  if (class(ata_output)!="ATA"){
    return("The Input must be 'ATA' object. Please use ATA(x) function to produce 'ATA' object. ATA Forecast was terminated!")
  }
  tsp_X <- tsp(ata_output$actual)
  X <- as.numeric(ata_output$actual)
  ph <- ata_output$p
  qh <- ata_output$q
  phih <- ata_output$phi
  modelType <- ata_output$model.type
  lenX <- length(X)
  if(is.null(hh)){
    hh <- ata_output$h
  }
  if (initialLevel==TRUE){
	Xobs <- mean(X)
  }else {
	Xobs <- X[lenX]
  }
  ata.forecast.fitted <- rep(NA, hh)
  if (modelType=="A"){
    coefph <- abs(ph/lenX)
    coefqh <- abs(qh/lenX)
    T_1 <- ata_output$trend[lenX-1]
    S_1 <- ata_output$level[lenX-1]
    ata_output$level[lenX] <- S <- coefph * Xobs + (1-coefph)*(S_1 + phih * T_1)
    ata_output$trend[lenX] <- T <- coefqh * (S-S_1) + (1-coefqh) * (phih * T_1)
    ata.forecast.fitted[1] <- S + (phih * T)
    phiTotal <- phih
    for (h in 2:hh){
      phiTotal <- phiTotal + (phih^h)
      ata.forecast.fitted[h] <- S + (phiTotal * T)
    }
  }
  if (modelType=="M"){
    coefph <- abs(ph/lenX)
    coefqh <- abs(qh/lenX)
    T_1 <- ata_output$trend[lenX-1]
    S_1 <- ata_output$level[lenX-1]
    ata_output$level[lenX] <- S <- coefph * Xobs + (1-coefph)* S_1 * (T_1^phih)
    ata_output$trend[lenX] <- T <- coefqh * (S/S_1) + (1-coefqh) * (T_1^phih)
    ata.forecast.fitted[1] <- S * (T^phih)
    phiTotal <- phih
    for (h in 2:hh){
      phiTotal <- phiTotal + (phih^h)
      ata.forecast.fitted[h] <- S * (T^phiTotal)
    }
  }
  ata.forecast.fitted <- ts(ata.forecast.fitted, frequency = tsp_X[3], start = tsp_X[2] + ifelse(tsp_X[3]>1, 1/tsp_X[3], 1))
  my_list <- ata_output
  my_list$forecast <- ata.forecast.fitted
  return(my_list)
}
