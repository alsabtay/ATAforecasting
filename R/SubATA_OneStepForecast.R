#' @importFrom stats ts tsp tsp<-
SubATA.OneStepForecast <- function(ataModel, outSample, hh=NULL)
{
  X <- as.numeric(ataModel$actual)
  pk <- ataModel$p
  qk <- ataModel$q
  phik <- ataModel$phi
  modelType <- ataModel$model.type
  lenX <- length(X)
  if(!is.null(outSample)){
    if(!is.na(outSample[1])){
      hh <- length(outSample)
      ataModel$h <- hh
      outflag <- TRUE
    }
  }
  if(is.null(hh)){
    hh <- ataModel$h
  }
  onestep.X <- rep(NA, lenX + hh)
  onestep.X[1:lenX] <- as.numeric(ataModel$actual)
  if(outflag){
    onestep.X[(lenX + 1):(lenX + hh)] <- outSample
  }
  onestep.S <- rep(NA, lenX + hh)
  onestep.S[1:lenX] <- as.numeric(ataModel$level)
  onestep.T <- rep(NA, lenX + hh)
  onestep.T[1:lenX] <- as.numeric(ataModel$trend)
  onestep.fitted <- rep(NA, lenX + hh)
  onestep.fitted[1:lenX] <- as.numeric(ataModel$fitted)
  onestep.coefp <- rep(NA, lenX + hh)
  onestep.coefp[1:lenX] <- as.numeric(ataModel$coefp)
  onestep.coefq <- rep(NA, lenX + hh)
  onestep.coefq[1:lenX] <- as.numeric(ataModel$coefq)

  T_1 <- onestep.T[lenX-1]
  S_1 <- onestep.S[lenX-1]

  for(i in lenX:(lenX + hh - 1)){
    Xobs = onestep.X[i]
    if (modelType=="A"){
      onestep.coefp[i] <- coefpk <- abs(pk/i)
      onestep.coefq[i] <- coefqk <- abs(qk/i)
      onestep.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * (S_1 + (phik * T_1)))
      onestep.T[i] <- T <- (coefqk * (S-S_1)) + ((1-coefqk) * phik * T_1)
      if(outflag){
        onestep.fitted[i+1] <- S + (phik * T)
      }else{
        onestep.fitted[i+1] <- onestep.X[i+1] <- S + (phik * T)
      }
      S_1 <- S
      T_1 <- T
    }
    if (modelType=="M"){
      onestep.coefp[i] <- coefpk <- abs(pk/i)
      onestep.coefq[i] <- coefqk <- abs(qk/i)
      onestep.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * S_1 * (T_1^phik))
      onestep.T[i] <- T <- (coefqk * (S/S_1)) + ((1-coefqk) * (T_1^phik))
      if(outflag){
        onestep.fitted[i+1] <- S * (T^phik)
      }else{
        onestep.fitted[i+1] <- onestep.X[i+1] <- S * (T^phik)
      }
      S_1 <- S
      T_1 <- T
    }
  }
  my_list <- ataModel
  if(outflag){
    my_list$forecast <- onestep.fitted[(lenX + 1):(lenX + hh)]
    my_list$onestep.forecast <- list("level" = onestep.S[lenX:(lenX + hh - 1)],
                                     "trend" = onestep.T[lenX:(lenX + hh - 1)],
                                     "coefp" = onestep.coefp[lenX:(lenX + hh - 1)],
                                     "coefq" = onestep.coefq[lenX:(lenX + hh - 1)])
  }else{
    my_list$onestep.forecast <- NA
  }
  return(my_list)
}
