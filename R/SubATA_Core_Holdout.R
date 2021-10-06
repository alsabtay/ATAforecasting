SubATA.Core.Holdout <- function(X, pk, qk, phik, mdlType, initialLevel, initialTrend, hh)
{
  X <- as.numeric(X)
  lenX <- length(X)
  ata.fitted <- rep(NA, lenX)
  ata.forecast <- rep(NA, hh)

  if (initialTrend==TRUE){
    if (mdlType=="A"){
      IT_0 <- X-ATA.Shift(X,1)
    }else {
      IT_0 <- X/ATA.Shift(X,1)
    }
  }
  for(i in 1:lenX-1){
    if (i==1)
    {
      Xlag = X[i]
      Xobs = X[i]
    }
    else
    {
      if (initialLevel==TRUE)
      {
        Xlag =  mean(X[1:i-1])
        Xobs =  mean(X[1:i])
      }
      else
      {
        Xlag = X[i-1]
        Xobs = X[i]
      }
    }

    if (mdlType=="A") {
      T_0 = Xobs - Xlag
    }else {
      T_0 = Xobs / Xlag
    }

    if (i == 1){
      if (mdlType=="A"){
        S <- Xobs
        T <- 0
        ata.fitted[i] <- S + (phik * T)
        S_1 <- S
        T_1 <- T
      }
      if (mdlType=="M"){
        S <- Xobs
        T <- 1
        ata.fitted[i] <- S * (T^phik)
        S_1 <- S
        T_1 <- T
      }
    }else if (i<=pk & i<=qk & pk>=qk){
      if (mdlType=="A"){
        S <- Xobs
        T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- S + (phik * T)
        S_1 <- S
        T_1 <- T
      }
      if (mdlType=="M"){
        S <- Xobs
        T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- S * (T^phik)
        S_1 <- S
        T_1 <- T
      }
    }else if (i<=pk & i>qk & pk>=qk){
      if (mdlType=="A"){
        coefqk <- abs(qk/i)
        S <- Xobs
        T <- coefqk * (S-S_1) + (1-coefqk) * phik * T_1
        ata.fitted[i] <- S + (phik * T)
        S_1 <- S
        T_1 <- T
      }
      if (mdlType=="M"){
        coefqk <- abs(qk/i)
        S <- Xobs
        T <- coefqk * (S/S_1) + (1-coefqk) * (T_1^phik)
        ata.fitted[i] <- S * (T^phik)
        S_1 <- S
        T_1 <- T
      }
    }else if (i>pk & i<=qk & pk>=qk){
      Xobs = X[i]
      if (mdlType=="A"){
        coefpk <- abs(pk/i)
        S <- coefpk * Xobs + (1-coefpk)*(S_1 + phik * T_1)
        T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- S + (phik * T)
        S_1 <- S
        T_1 <- T
      }
      if (mdlType=="M"){
        coefpk <- abs(pk/i)
        S <- coefpk * Xobs + (1-coefpk)* S_1 * (T_1^phik)
        T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- S * (T^phik)
        S_1 <- S
        T_1 <- T
      }
    }else if (i>pk & i>qk & pk>=qk){
      Xobs = X[i]
      if (mdlType=="A"){
        coefpk <- abs(pk/i)
        coefqk <- abs(qk/i)
        S <- coefpk * Xobs + (1-coefpk)*(S_1 + phik * T_1)
        T <- coefqk * (S-S_1) + (1-coefqk) * phik * T_1
        ata.fitted[i] <- S + (phik * T)
        S_1 <- S
        T_1 <- T
      }
      if (mdlType=="M"){
        coefpk <- abs(pk/i)
        coefqk <- abs(qk/i)
        S <- coefpk * Xobs + (1-coefpk)* S_1 * (T_1^phik)
        T <- coefqk * (S/S_1) + (1-coefqk) * (T_1^phik)
        ata.fitted[i] <- S * (T^phik)
        S_1 <- S
        T_1 <- T
      }
    }else {
      ata.fitted[i] <- NA
      S_1 <- NA
      T_1 <- NA
    }
  }

  if (initialLevel==TRUE){
    Xobs <- mean(X)
  }else {
    Xobs <- X[lenX]
  }
  if (mdlType=="A"){
    coefpk <- abs(pk/lenX)
    coefqk <- abs(qk/lenX)
    S <- coefpk * Xobs + (1-coefpk)*(S_1 + phik * T_1)
    T <- coefqk * (S-S_1) + (1-coefqk) * (phik * T_1)
    ata.forecast[1] <- S + (phik * T)
    phiTotal <- phik
    for (h in 2:hh){
      phiTotal <- phiTotal + (phik^h)
      ata.forecast[h] <- S + (phiTotal * T)
    }
  }
  if (mdlType=="M"){
    coefpk <- abs(pk/lenX)
    coefqk <- abs(qk/lenX)
    S <- coefpk * Xobs + (1-coefpk)* S_1 * (T_1^phik)
    T <- coefqk * (S/S_1) + (1-coefqk) * (T_1^phik)
    ata.forecast[1] <- S * (T^phik)
    phiTotal <- phik
    for (h in 2:hh){
      phiTotal <- phiTotal + (phik^h)
      ata.forecast[h] <- S * (T^phiTotal)
    }
  }
  my_list <- list("actual" = X, "fitted" = ata.fitted, "forecast" = ata.forecast)
  return(my_list)
}
