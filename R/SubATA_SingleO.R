#' @importFrom forecast msts
#' @importFrom stats frequency ts tsp tsp<- var
SubATA.SingleO <- function(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                            level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method, lambda, shift, orig.X,
                            OutSample, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, partition.h, holdout.adjustedP, holdin, nmse)
{
  tspX <- tsp(X)
  firstTspX <- tsp(orig.X)
  if (seasonal.model=="none"){
    is.season <- FALSE
    seasonal.type <- "A"
  }else {
    if (seasonal.test==FALSE){
      if (max(s.frequency)==1){
        is.season <- FALSE
        seasonal.type <- "A"
      }else {
        is.season <- TRUE
      }
    }else {
      is.season <- ATA.Seasonality(X, s.frequency, seas_attr_set)
    }
  }
  if (is.season==TRUE){
    if (seasonal.model!="decomp" & seasonal.type=="M"){
      out.transform <- ATA.Transform(X, tMethod = "Box_Cox", tLambda = 0, tShift = 0)  # lambda = 0 for multiplicative model
	  X <- forecast::msts(out.transform$trfmX, start = tspX[1], seasonal.periods = s.frequency)
      seas.type <- "A"
      seas.lambda <- out.transform$tLambda
	  seas.shift <- out.transform$tShift
      seas.transform <- "Box_Cox"
    }else {
      seas.lambda <- NULL
      seas.transform <- NULL
      seas.type <- seasonal.type
	  seas.shift <- 0
    }
  }else {
    seasonal.model <- "none"
    seas.type <- seasonal.type <- "A"
    seas.lambda <- NULL
	seas.shift <- 0
	seas.transform <- NULL
  }
  ata.seasonal.component <- ATA.Decomposition(X, s.model=seasonal.model, s.type=seas.type, s.frequency=s.frequency, seas_attr_set=seas_attr_set)
  AdjInSample <- ATA.BackTransform(X=ata.seasonal.component$AdjustedX, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  SeasonalIndex <- ATA.BackTransform(X=ata.seasonal.component$SeasIndex, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  SeasonalActual <- ATA.BackTransform(X=ata.seasonal.component$SeasActual, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  if (seasonal.model=="x13" | seasonal.model=="x11"){
    seasonal.type <- ata.seasonal.component$SeasType
  }
  ChgX <- ATA.Transform(AdjInSample, tMethod = transform.method, tLambda = lambda, tShift = shift, bcMethod = boxcox_attr_set$bcMethod, bcLower = boxcox_attr_set$bcLower, bcUpper = boxcox_attr_set$bcUpper)
  AdjInSample <- ChgX$trfmX
  lambda <- ChgX$tLambda
  shift <- ChgX$tShift
  if (is.season==FALSE & seasonal.type=="A"){
    OS_SIValue <- rep(0,times=h)
  }else if (is.season==FALSE & seasonal.type=="M"){
    OS_SIValue <- rep(1,times=h)
  }else if (is.season==TRUE){
    OS_SIValue <- rep(NA,times=h)
    for (k in 1:h){
      OS_SIValue[k] <- SeasonalIndex[freqYh[k]]
    }
  }else{
  }
  if (holdout == TRUE){
    holdout_part <- ifelse(partition.h > 0 & partition.h < 1, floor(length(AdjInSample) * partition.h), partition.h)
    HoldOutLen <- length(AdjInSample) - holdout_part
    InsampleLen <- length(AdjInSample)
    HoldoutSet <- ts(AdjInSample[(HoldOutLen+1):InsampleLen], frequency = tspX[3], start = tspX[2] - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1))
    DeSeas <- ts(AdjInSample[1:HoldOutLen], frequency = tspX[3], start = tspX[1])
  }else {
    DeSeas <- AdjInSample
    HoldoutSet <- NA
  }
  ata.output <- SubATA.Damped(DeSeas, pb = parP, qb = parQ, model.Type = model.type, accuracy.Type = accuracy.type, level.fix = level.fixed, trend.fix = trend.fixed, trend.Search = trend.search, phiStart = start.phi, phiEnd = end.phi, phiSize = size.phi, initialLevel = initial.level, initialTrend = initial.trend, orig_X = AdjInSample, Holdout = holdout, HoldoutSet = HoldoutSet, Adjusted_P = holdout.adjustedP, h = h, Holdin = holdin, nmse = nmse)
  ata.output$h <- h
  ata.output <- SubATA.Forecast(ata.output, hh=h, initialLevel = initial.level)
  ata.output$actual <- orig.X
  fit.ata <- ATA.BackTransform(X=ata.output$fitted, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
  forecast.ata <- ATA.BackTransform(X=ata.output$forecast, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
  if (holdout == TRUE){
	holdout.ata <- ATA.BackTransform(X=ata.output$holdout.forecast, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
  }
  ata.output$level <- ts(ATA.BackTransform(X=ata.output$level, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))), frequency = firstTspX[3], start=firstTspX[1])
  ata.output$trend <- ts(ATA.BackTransform(X=ata.output$trend, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))), frequency = firstTspX[3], start=firstTspX[1])
  if(seasonal.type == "A"){
    ATA.fitted <- fit.ata + SeasonalActual
    ATA.forecast <- forecast.ata + OS_SIValue
    if (holdout == TRUE){
      ata.output$holdout.forecast <- ts(holdout.ata + SeasonalActual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
    }
  }else {
    ATA.fitted <- fit.ata * SeasonalActual
    ATA.forecast <- forecast.ata * OS_SIValue
    if (holdout == TRUE){
      ata.output$holdout.forecast <- ts(holdout.ata * SeasonalActual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
    }
  }
  ata.output$fitted <- ts(ATA.fitted, frequency = firstTspX[3], start=firstTspX[1])
  if (negative.forecast==TRUE){
    ata.output$forecast <- ts(ATA.forecast, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    ata.output$forecast <- ts(ATA.forecast, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }
  ata.output$residuals <- ata.output$actual - ata.output$fitted
  if (holdout == TRUE){
	ata.output$holdout.training <- ts(ata.output$actual[1:HoldOutLen], frequency = firstTspX[3], start = firstTspX[1])
    ata.output$holdout.validation <- ts(ata.output$actual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = tsp(ata.output$holdout.training)[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }
  my_list <- ata.output
  my_list$out.sample <- ts(OutSample, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  if (level.fixed==TRUE){
    method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
  }else if (trend.fixed==TRUE){
    method <- paste("ATA(", my_list$p, ",1," ,my_list$phi, ")", sep="")
  }else if (trend.search==TRUE){
    method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
  }else {
    method <- paste("ATA(", my_list$p, "," ,my_list$q, ",", my_list$phi, ")", sep="")
  }
  my_list$initial.level <- initial.level
  my_list$initial.trend <- initial.trend
  my_list$level.fixed <- level.fixed
  my_list$trend.fixed <- trend.fixed
  my_list$trend.search <- trend.search
  my_list$transform.method <- transform.method
  my_list$lambda <- lambda
  my_list$shift <- shift
  my_list$bcLower <- boxcox_attr_set$bcLower
  my_list$bcUpper <- boxcox_attr_set$bcUpper
  my_list$bcBiasAdj <- boxcox_attr_set$bcBiasAdj
  my_list$accuracy.type <- accuracy.type
  my_list$nmse <- nmse
  my_list$is.season <- is.season
  my_list$seasonal.model <- seasonal.model
  my_list$seasonal.type <- seasonal.type
  if(my_list$q==0){
    trend_mthd <- "N"
  }else if (my_list$q!=0 & my_list$phi!=1 & my_list$phi>0){
    trend_mthd <- paste(my_list$model.type, "d", sep="")
  }else{
    trend_mthd <- my_list$model.type
  }
  if(my_list$seasonal.model == "none"){
    seas_mthd <- "N"
  }else{
    seas_mthd <- my_list$seasonal.type
  }
  method <- paste(method, " (A,", trend_mthd, ",", seas_mthd, ")", sep="")
  my_list$method <- method
  my_list$par.specs <- list("p" = my_list$p, "q" = my_list$q, "phi" = my_list$phi,
                              "trend" = trend_mthd,
                              "seasonal" = seas_mthd,
                              "initial_level" = ifelse(my_list$initial.level==FALSE, NA, TRUE),
                              "initial_trend" = ifelse(my_list$initial.trend==FALSE, NA, TRUE))
  accuracy.ata <- ATA.Accuracy(my_list, OutSample)
  my_list$accuracy <- accuracy.ata
  my_list$seasonal.period <- s.frequency
  my_list$seasonal.index <- SeasonalIndex
  my_list$seasonal <- ts(SeasonalActual, frequency = firstTspX[3], start=firstTspX[1])
  my_list$seasonal.adjusted <- ts(ATA.BackTransform(X=AdjInSample, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))), frequency = firstTspX[3], start=firstTspX[1])
  ci.output <- ATA.CI(object = my_list, ci.level = ci.level)
  my_list$ci.level <- ci.level
  if (negative.forecast==TRUE){
    my_list$forecast.lower <- ci.output$forecast.lower
    my_list$forecast.upper <- ci.output$forecast.upper
  }else {
    ci_low <- ci.output$forecast.lower
    ci_up <- ci.output$forecast.upper
    ci_low[ci_low<0] <- 0
    ci_up[ci_up<0] <- 0
    my_list$forecast.lower <- ci_low
    my_list$forecast.upper <- ci_up
  }
  return(my_list)
}
