#' The core algorithm of the ATA Method
#'
#' @param X A numeric vector or time series.
#' @param pk Value of Level parameter.
#' @param qk Value of Trend parameter.
#' @param phik Value of Damping Trend parameter.
#' @param mdlType An one-character string identifying method using the framework terminology.
#' @param initialLevel If NULL, FALSE is default. If FALSE, ATA Method calculates the pth observation in \code{X} for level.
#' If TRUE, ATA Method calculates average of first p value in \code{X}for level.
#' @param initialTrend If NULL, FALSE is default. If FALSE, ATA Method calculates the qth observation in \code{X(T)-X(T-1)} for trend.
#' If TRUE, ATA Method calculates average of first q value in \code{X(T)-X(T-1)} for trend.
#' @param nmse If `accuracy.type == "AMSE"`, `nmse` provides the number of steps for average multistep MSE (`2<=nmse<=30`).
#'
#' @return Returns an object of class "\code{ATA}"
#'
#' @importFrom forecast msts
#' @importFrom stats as.ts ts tsp tsp<-
#'
#' @export
ATA.Core <- function(X, pk, qk, phik, mdlType, initialLevel, initialTrend, nmse)
{
  tspX <- tsp(X)
  lenX <- length(X)
  if ("msts" %in% class(X)) {
       X_msts <- attributes(X)$msts
       if (any(X_msts >= lenX / 2)) {
         X_msts <- X_msts[X_msts < lenX / 2]
       }
       X_msts <- sort(X_msts, decreasing = FALSE)
  }else if ("ts" %in% class(X)) {
       X_ts <- frequency(X)
  }else {
       X_ts <- 1L
  }
  X <- as.numeric(X)
  ata.S <- rep(NA, lenX)
  ata.T <- rep(NA, lenX)
  ata.error <- rep(NA, lenX)
  ata.fitted <- rep(NA, lenX)
  ata.coefp <- rep(NA, lenX)
  ata.coefq <- rep(NA, lenX)
  FC <- matrix(NA,nrow=lenX, ncol=nmse)
  S_1 <- NA
  T_1 <- NA

  if (initialTrend==TRUE){
    if (mdlType=="A"){
      IT_0 <- X-ATA.Shift(X,1)
    }else {
      IT_0 <- X/ATA.Shift(X,1)
    }
  }

  for(i in 1:(lenX-1)){
    Xh = X[i+1]
    if (i==1) {
      Xlag = X[i]
      Xobs = X[i]
    }else {
      if (initialLevel==TRUE){
        Xlag =  mean(X[1:i-1])
        Xobs =  mean(X[1:i])
      }else {
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
        ata.coefp[i] <- NA
        ata.coefq[i] <- NA
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- 0
        ata.fitted[i] <- FF <- S + (phik * T)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S + (phiTotal * T)
        }
       }
      if (mdlType=="M"){
        ata.coefp[i] <- NA
        ata.coefq[i] <- NA
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- 1
        ata.fitted[i] <- FF <- S * (T^phik)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S * (T^phiTotal)
        }
      }
    }else if (i<=pk & i<=qk & pk>=qk){
      if (mdlType=="A"){
        ata.coefp[i] <- NA
        ata.coefq[i] <- NA
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- FF <- S + (phik * T)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S + (phiTotal * T)
        }
      }
      if (mdlType=="M"){
        ata.coefp[i] <- NA
        ata.coefq[i] <- NA
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- FF <- S * (T^phik)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S * (T^phiTotal)
        }
      }
    }else if (i<=pk & i>qk & pk>=qk){
      if (mdlType=="A"){
        ata.coefp[i] <- NA
        ata.coefq[i] <- coefqk <- abs(qk/i)
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- (coefqk * (S-S_1)) + ((1-coefqk) * phik * T_1)
        ata.fitted[i] <- FF <- S + (phik * T)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for(j in 2:nmse) {
            phiTotal = phiTotal + (phik^j)
            FC[i,j] = S + (phiTotal * T)
        }
      }
      if (mdlType=="M"){
        ata.coefp[i] <- NA
        ata.coefq[i] <- coefqk <- abs(qk/i)
        ata.S[i] <- S <- Xobs
        ata.T[i] <- T <- (coefqk * (S/S_1)) + ((1-coefqk) * (T_1^phik))
        ata.fitted[i] <- FF <- S * (T^phik)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S * (T^phiTotal)
        }
      }
    }else if (i>pk & i<=qk & pk>=qk){
      Xobs = X[i]
      if (mdlType=="A"){
        ata.coefp[i] <- coefpk <- abs(pk/i)
        ata.coefq[i] <- NA
        ata.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * (S_1 + (phik * T_1)))
        ata.T[i] <- T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- FF <- S + (phik * T)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S + (phiTotal * T)
        }
      }
      if (mdlType=="M"){
        ata.coefp[i] <- coefpk <- abs(pk/i)
        ata.coefq[i] <- NA
        ata.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * S_1 * (T_1^phik))
        ata.T[i] <- T <- ifelse(initialTrend==TRUE, mean(IT_0[1:i]),T_0)
        ata.fitted[i] <- FF <- S * (T^phik)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S * (T^phiTotal)
        }
      }
    }else if (i>pk & i>qk & pk>=qk){
      Xobs = X[i]
      if (mdlType=="A"){
        ata.coefp[i] <- coefpk <- abs(pk/i)
        ata.coefq[i] <- coefqk <- abs(qk/i)
        ata.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * (S_1 + (phik * T_1)))
        ata.T[i] <- T <- (coefqk * (S-S_1)) + ((1-coefqk) * phik * T_1)
        ata.fitted[i] <- FF <- S + (phik * T)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S + (phiTotal * T)
        }
      }
      if (mdlType=="M"){
        ata.coefp[i] <- coefpk <- abs(pk/i)
        ata.coefq[i] <- coefqk <- abs(qk/i)
        ata.S[i] <- S <- (coefpk * Xobs) + ((1-coefpk) * S_1 * (T_1^phik))
        ata.T[i] <- T <- (coefqk * (S/S_1)) + ((1-coefqk) * (T_1^phik))
        ata.fitted[i] <- FF <- S * (T^phik)
        ata.error[i] <- Xh - FF
        S_1 <- S
        T_1 <- T
        FC[i,1] <- FF
        phiTotal <- phik
        for (j in 2:nmse) {
            phiTotal <- phiTotal + (phik^j)
            FC[i,j] <- S * (T^phiTotal)
        }
      }
    }else {
      ata.coefp[i] <- NA
      ata.coefq[i] <- NA
      ata.S[i] <- NA
      ata.T[i] <- NA
      ata.fitted[i] <- NA
      ata.error[i] <- NA
      S_1 <- NA
      T_1 <- NA
    }
  }
  ata.fitted <- ATA.Shift(ata.fitted,-1)
  ata.error <- ATA.Shift(ata.error,-1)
  FC <- ATA.Shift_Mat(FC, "down", 1)
  if ("msts" %in% class(X)) {
         X <- forecast::msts(X, start = tspX[1], seasonal.periods = X_msts)
         ata.fitted <- forecast::msts(ata.fitted, start = tspX[1], seasonal.periods = X_msts)
         ata.error <- forecast::msts(ata.error, start = tspX[1], seasonal.periods = X_msts)
         ata.S <- forecast::msts(ata.S, start = tspX[1], seasonal.periods = X_msts)
         ata.T <- forecast::msts(ata.T, start = tspX[1], seasonal.periods = X_msts)
         ata.coefp <- forecast::msts(ata.coefp, start = tspX[1], seasonal.periods = X_msts)
         ata.coefq <- forecast::msts(ata.coefq, start = tspX[1], seasonal.periods = X_msts)
  }else {
       X <- ts(X, frequency = X_ts, start = tspX[1])
       ata.fitted <- ts(ata.fitted, frequency = X_ts, start = tspX[1])
       ata.error <- ts(ata.error, frequency = X_ts, start = tspX[1])
       ata.S <- ts(ata.S, frequency = X_ts, start = tspX[1])
       ata.T <- ts(ata.T, frequency = X_ts, start = tspX[1])
       ata.coefp <- ts(ata.coefp, frequency = X_ts, start = tspX[1])
       ata.coefq <- ts(ata.coefq, frequency = X_ts, start = tspX[1])
  }
  my_list <- list("actual" = X, "fitted" = ata.fitted , "level" = ata.S, "trend" = ata.T, "residuals" = ata.error, "coefp" = ata.coefp, "coefq" = ata.coefq,
                  "p" = as.integer(pk), "q" = as.integer(qk), "phi" = signif(phik,6), "model.type" = mdlType, "amse.fc" = FC)
  attr(my_list, 'class') <- "ata"
  return(my_list)
}
