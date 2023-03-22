#' @importFrom stats ts tsp tsp<-
SubATA.Forecast <- function(ataModel, hh=NULL)
{
  X <- as.numeric(ataModel$actual)
  ph <- ataModel$p
  qh <- ataModel$q
  phih <- ataModel$phi
  modelType <- ataModel$model.type
  lenX <- length(X)
  if(is.null(hh)){
    hh <- ataModel$h
  }
	Xobs <- X[lenX]
  multistep.fitted <- rep(NA, hh)
  if (modelType=="A"){
    coefph <- abs(ph/lenX)
    coefqh <- abs(qh/lenX)
    T_1 <- ataModel$trend[lenX-1]
    S_1 <- ataModel$level[lenX-1]
    ataModel$level[lenX] <- S <- coefph * Xobs + (1-coefph)*(S_1 + phih * T_1)
    ataModel$trend[lenX] <- T <- coefqh * (S-S_1) + (1-coefqh) * (phih * T_1)
    multistep.fitted[1] <- S + (phih * T)
    phiTotal <- phih
    if (hh > 1){
      for (h in 2:hh){
        phiTotal <- phiTotal + (phih^h)
        multistep.fitted[h] <- S + (phiTotal * T)
      }
    }
  }
  if (modelType=="M"){
    coefph <- abs(ph/lenX)
    coefqh <- abs(qh/lenX)
    T_1 <- ataModel$trend[lenX-1]
    S_1 <- ataModel$level[lenX-1]
    ataModel$level[lenX] <- S <- coefph * Xobs + (1-coefph)* S_1 * (T_1^phih)
    ataModel$trend[lenX] <- T <- coefqh * (S/S_1) + (1-coefqh) * (T_1^phih)
    multistep.fitted[1] <- S * (T^phih)
    phiTotal <- phih
    if (hh > 1){
      for (h in 2:hh){
        phiTotal <- phiTotal + (phih^h)
        multistep.fitted[h] <- S * (T^phiTotal)
      }
    }
  }
  my_list <- ataModel
  my_list$forecast <- multistep.fitted
  my_list$onestep.forecast <- NA
  return(my_list)
}
