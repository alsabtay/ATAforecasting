#' @importFrom stats ts tsp tsp<-
SubATA.Forecast <- function(ata_output, hh=NULL, initialLevel)
{
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
    if (hh > 1){
      for (h in 2:hh){
        phiTotal <- phiTotal + (phih^h)
        ata.forecast.fitted[h] <- S + (phiTotal * T)
      }
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
    if (hh > 1){
      for (h in 2:hh){
        phiTotal <- phiTotal + (phih^h)
        ata.forecast.fitted[h] <- S * (T^phiTotal)
      }
    }
  }
  my_list <- ata_output
  my_list$forecast <- ata.forecast.fitted
  return(my_list)
}
