#' Transformation Techniques for ATA Method
#'
#' @description These functions are modified version of \code{trafo::trafos} written by Lily Medina, Piedad Castro, Ann-Kristin Kreutzmann, Natalia Rojas-Perilla \code{trafo} package.
#' Please review manual and vignette documents of latest \code{trafo} package. The \code{ATA.Transform} function works with many different types of inputs.
#'
#' @param X a numeric vector or time series of class \code{ts} or \code{msts} for in-sample.
#' @param tMethod Box-Cox power transformation family is consist of "BoxCox", "BoxCox Shift", "Sqrt", "Sqrt Shift", "Reciprocal", "Log", "Log Shift", "NegLog", "Modulus", "Bickel-Doksum", "Manly", "Dual", "Yeo-Johnson", "GPower", "GLog" in ATAforecasting package.
#' @param tLambda Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.
#' @param tShift Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated.
#' @param bcMethod Choose method to be used in calculating lambda. "loglik" is default. Other method is "guerrero" (Guerrero, V.M. (1993)).
#' @param bcLower Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0.
#' @param bcUpper Upper limit for possible lambda values. The upper value is limited by 5. Default value is 1.
#'
#' @return A list object consists of transformation parameters and transformed data.
#' \code{ATA.Transform} is a list containing at least the following elements:
#' \itemize{
#'		 \item{trfmX}   : Transformed data
#'		 \item{tLambda} : Box-Cox power transformation family parameter
#'		 \item{tShift}  : Box-Cox power transformation family shifting parameter
#'}
#'
#' @importFrom forecast BoxCox BoxCox.lambda InvBoxCox
#'
#' @export
ATA.Transform <- function(X
                          , tMethod = c("BoxCox", "BoxCox Shift", "Sqrt", "Sqrt Shift", "Reciprocal", "Log", "Log Shift", "NegLog",
                                        "Modulus", "Bickel-Doksum", "Manly", "Dual", "Yeo-Johnson", "GPower", "GLog")
                          , tLambda
                          , tShift = 0
                          , bcMethod = c("loglik", "guerrero")
                          , bcLower = 0
                          , bcUpper = 1)
{
  if (is.null(tMethod)){
    trfmX <- X
  }else {
    if (is.null(tLambda)){
      bcMethod <- match.arg(bcMethod)
      if (bcMethod == "guerrero") {
        tLambda <- forecast::BoxCox.lambda(X, method = "guerrero", lower=bcLower, upper=bcUpper)
      } else {
        tLambda <- forecast::BoxCox.lambda(X, method = "loglik", lower=bcLower, upper=bcUpper)
      }
    }
    if (tMethod=="BoxCox"){
      trfmX <- forecast::BoxCox(X,tLambda)
    }else if (tMethod=="BoxCox Shift"){
      tZ <- box_cox_shift(X, lambda = tLambda, shift = tShift)
      trfmX <- tZ$tX
      tShift <- tZ$shift
    }else if (tMethod=="Sqrt"){
      trfmX <- sqrt(X)
      tLambda <- NA
    }else if (tMethod=="Sqrt Shift"){
      tZ <- Sqrt_shift(X, tShift)
      trfmX <- tZ$tX
      tShift <- tZ$shift
      tLambda <- NA
    }else if (tMethod=="Reciprocal"){
      tLambda <- -1
      trfmX <- forecast::BoxCox(X,tLambda)
    }else if (tMethod=="Log"){
      trfmX <- box_cox(X, lambda = 0)
      tLambda <- 0
    }else if (tMethod=="Log Shift"){
      tZ <- box_cox_shift(X, lambda = 0, shift = tShift)
      trfmX <- tZ$tX
      tShift <- tZ$shift
      tLambda <- 0
    }else if (tMethod=="NegLog"){
      trfmX <- Neg_Log(X)
    }else if (tMethod=="Modulus"){
      trfmX <- Modulus(X, tLambda)
    }else if (tMethod=="Bickel-Doksum"){
      trfmX <- Bickel_Doksum(X, tLambda)
    }else if (tMethod=="Manly"){
      trfmX <- Manly(X, tLambda)
    }else if (tMethod=="Dual"){
      trfmX <- Dual(X, tLambda)
    }else if (tMethod=="Yeo-Johnson"){
      trfmX <- Yeo_Johnson(X, tLambda)
    }else if (tMethod=="GPower"){
      trfmX <- GPower(X, tLambda)
    }else if (tMethod=="GLog"){
      trfmX <- GLog(X)
    }else {
      trfmX <- X
      tLambda <- NA
    }
  }
  my_list <- list("trfmX" = trfmX, "tLambda" = tLambda, "tShift" = tShift)
  return(my_list)
}

ATA.BackTransform <- function(X, tMethod, tLambda, tShift, tbiasadj=FALSE, tfvar=NULL){
  if (is.null(tMethod)){
    trfmX <- X
  }else if (tMethod == "BoxCox") {
    trfmX <- forecast::InvBoxCox(X, lambda=tLambda, biasadj=tbiasadj, fvar=tfvar)
  }else if (tMethod=="BoxCox Shift"){
    trfmX <- box_cox_shift_back(X, lambda = tLambda, shift = tShift)
  }else if (tMethod=="Sqrt"){
    trfmX <- X^2
  }else if (tMethod=="Sqrt Shift"){
    trfmX <- Sqrt_shift_back(X, shift = tShift)
  }else if (tMethod=="Reciprocal"){
    trfmX <- box_cox_back(X, lambda = -1)
  }else if (tMethod=="Log"){
    trfmX <- box_cox_back(X, lambda = 0)
  }else if (tMethod=="Log Shift"){
    trfmX <- box_cox_shift_back(X, lambda = 0)
  }else if (tMethod=="NegLog"){
    trfmX <- Neg_Log_back(X)
  }else if (tMethod=="Modulus"){
    trfmX <- Modulus_back(X, tLambda)
  }else if (tMethod=="Bickel-Doksum"){
    trfmX <- Bickel_Doksum_back(X, tLambda)
  }else if (tMethod=="Manly"){
    trfmX <- Manly_back(X, tLambda)
  }else if (tMethod=="Dual"){
    trfmX <- Dual_back(X, tLambda)
  }else if (tMethod=="Yeo-Johnson"){
    trfmX <- Yeo_Johnson_back(X, tLambda)
  }else if (tMethod=="GPower"){
    trfmX <- GPower_back(X, tLambda)
  }else if (tMethod=="GLog"){
    trfmX <- GLog_back(X)
  }else {
    trfmX <- X
  }
  return(trfmX)
}

# Box Cox ----------------------------------------------------------------------
# Transformation: Box Cox
box_cox <- function(X, lambda) {
  lambda_cases <- function(X, lambda) {
    lambda_absolute <- abs(lambda)
    if (lambda_absolute <= 1e-12) {  #case lambda=0
      yt <- log(X)
    }else {
      yt <- ((X)^lambda - 1) / lambda
    }
    return(yt)
  }
  zt <- lambda_cases(X = X, lambda = lambda)
  return(zt)
}

# Back transformation: Box Cox
box_cox_back <- function(X, lambda) {
  lambda_cases_back <- function(X, lambda){
    if (abs(lambda) <= 1e-12) {   #case lambda=0
      yt <-  exp(X)
    }else {
      yt <- (lambda * X + 1)^(1 / lambda)
    }
    return(yt)
  }
  zt <- lambda_cases_back(X = X, lambda = lambda)
  return(zt)
}

#  Transformation: Box Cox shift
box_cox_shift <- function(X, lambda = lambda, shift = 0) {
  with_shift <- function(X, shift) {
    min_X <- min(X)
    if (min_X <= 0) {
      shift_new <- shift + abs(min(X)) + 1
    }else {
      shift_new <- shift
    }
    return(shift_new)
  }
  # Shift parameter
  shift_new <- with_shift(X = X, shift = shift)
  lambda_cases <- function(X, lambda, shift) {
    lambda_absolute <- abs(lambda)
    if (lambda_absolute <= 1e-12) {  #case lambda=0
      yt <- log(X + shift)
    }else {
      yt <- ((X + shift)^lambda - 1) / lambda
    }
    return(yt)
  }
  zt <- lambda_cases(X = X, lambda = lambda, shift = shift_new)
  return(list("tX" = zt, "shift" = shift_new))
}

# Back transformation: Box Cox Shift
box_cox_shift_back <- function(X, lambda, shift = 0) {
  lambda_cases_back <- function(X, lambda, shift){
    if (abs(lambda) <= 1e-12) {   #case lambda=0
      yt <-  exp(X) - shift
    }else {
      yt <- (lambda * X + 1)^(1 / lambda) - shift
    }
    return(yt)
  }
  zt <- lambda_cases_back(X = X, lambda = lambda, shift = shift)
  return(zt)
}

# The Modulus transformation ----------------------------------------------------------------------

#  Transformation: Modulus
Modulus <- function(X, lambda) {
  u <- abs(X) + 1L
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <-  sign(X)*log(u)
  }else {
    yt <- sign(X)*(u^lambda - 1L)/lambda
  }
  return(yt)
}

# Back transformation: Modulus
Modulus_back <- function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {
    yt <- sign(X) * (exp(abs(X)) - 1)
  }else {
    yt <- sign(X) * ((abs(X)*lambda + 1)^(1/lambda) - 1)
  }
  return(yt)
}

# The Bickel-Doksum transformation ----------------------------------------------------------------------

#  Transformation: Bick-Doksum
Bickel_Doksum <-  function(X, lambda) {
  if (lambda > 1e-12){
    yt <- (abs(X)^lambda * sign(X) - 1)/lambda
  }else {
    stop("lambda must be positive for the Bickel-Doksum transformation")
  }
  return(yt)
}

# Back transformation: Bick-Doksum
Bickel_Doksum_back <- function(X, lambda) {
  positivos <- which(X >= 0)
  X[positivos] <- (lambda * X[positivos] + 1)^(1 / lambda)
  negativos <- which(X < 0)
  X[negativos] <- (-1) * ((-1) * (lambda * X[negativos] + 1))^(1 / lambda)
  return(X)
}

# The Manly transformation ----------------------------------------------------------------------

# Transformation: Manly
Manly <-  function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <-  X
  }else {
    yt <- (exp(X*lambda) - 1L)/lambda
  }
  return(yt)
}

# Back transformation: Manly
Manly_back <- function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <- X
  }else {
    yt <- log(lambda * X + 1) / lambda
  }
  return(yt)
}

# The Dual transformation ----------------------------------------------------------------------

# Transformation: Dual
Dual <-  function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <-  log(X)
  }else if (lambda > 1e-12){
    yt <- (X^(lambda) - X^(-lambda))/(2 * lambda)
  }else {
    stop("lambda can not be negative for the dual transformation")
  }
  return(yt)
}

# Back transformation: Dual
Dual_back <- function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if(lambda_absolute <= 1e-12) {
    yt <- exp(X)
  }else {
    yt <- (lambda * X + sqrt(lambda^2 * X^2 + 1))^(1/lambda)
  }
  return(yt)
}

# The Yeo-Johnson transformation ----------------------------------------------------------------------

# Transformation: Yeo-Johnson
Yeo_Johnson <-  function(X, lambda) {
  n <- length(X)
  yt <- rep(NA, n)
  negativos <- which(X < 0)
  positivos <- which(X >= 0)
  if(abs(lambda) <= 1e-12) {
    yt[positivos] <- log(X[positivos] + 1)
  }else {
    yt[positivos] <- ((X[positivos] + 1)^lambda - 1)/lambda
  }
  if(abs(lambda - 2) <= 1e-12) {
    yt[negativos] <- -log(-X[negativos] + 1)
  }else {
    yt[negativos] <- -((-X[negativos] + 1)^(2-lambda) - 1)/(2-lambda)
  }
  return(yt)
}

# Back transformation: Yeo-Johnson
Yeo_Johnson_back <- function(X, lambda) {
  negativos <- which(X < 0)
  positivos <- which(X >= 0)
  lambda_absolute <- abs(lambda)
  if (lambda != 0) {
    X[positivos] <- ((X[positivos] * lambda + 1)^(1 / lambda)) - 1
  }
  if (lambda_absolute <= 1e-12) {
    X[positivos] <- exp(X[positivos]) - 1
  }
  if (lambda != 2) {
    X[negativos] <- (-1) * ((X[negativos] * (lambda - 2) + 1)^(1/(2 - lambda)) - 1)
  }
  if (lambda_absolute == 2) {
    X[negativos] <- (-1) * (exp(-X[negativos]) - 1)
  }
  return(X)
}

#  Transformation: neg_log
Neg_Log <- function(X) {
  u <- abs(X) + 1L
  yt <-  sign(X)*log(u)
  return(yt)
}

# Back transformation: neg_log
Neg_Log_back <- function(X) {
  yt <- sign(X) * (exp(abs(X)) - 1)
  return(yt)
}

# Transformation: Squared Root shift
Sqrt_shift <- function(X, shift) {
  with_shift <-  function(X, shift) {
    min_X <- min(X)
    if (min_X <= 0) {
      shift_new <- shift + abs(min_X) + 1
    } else {
      shift_new <- shift
    }
    return(shift_new)
  }
  # Shift parameter
  shift_new <- with_shift(X = X, shift = shift)
  sqrt_ata <- function(X, shift) {
    yt <- sqrt(X + shift)
    return(yt)
  }
  zt <- sqrt_ata(X = X, shift = shift_new)
  return(list("tX" = zt, "shift" = shift_new))
}

# Back transformation: Squared Root shift
Sqrt_shift_back <- function(X, shift) {
  yt <-  X^2 - shift
  return(yt)
}

# Transformation: GPower
GPower <-  function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <-  log(X + sqrt(X^2 + 1))
  } else if (lambda_absolute > 1e-12) {
    yt <- ((X + sqrt(X^2 + 1))^lambda - 1)/lambda
  }
  return(yt)
}

# Back transformation: GPower
GPower_back <- function(X, lambda) {
  lambda_absolute <- abs(lambda)
  if (lambda_absolute <= 1e-12) {  #case lambda=0
    yt <- (-(1 - exp(X*2))) / (2 * exp(X))
  } else if (lambda_absolute > 1e-12) {
    A <- (X * lambda + 1)^(1 / lambda)
    yt <- (-(1 - A^2)) / (2*A)
  }
  return(yt)
}

# Transformation: GLog
GLog <- function(X) {
  yt <-  log(X + sqrt(X^2 + 1))
  return(yt)
}

# Back-transformation: GLog
GLog_back <- function(X) {
  yt <- (-(1 - exp(X*2))) / (2 * exp(X))
  return(yt)
}
