#' @importFrom forecast msts
#' @importFrom stats end start ts tsp tsp<- var
SubATA_Single_After <- function(train_set, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                            level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method, lambda, shift, main_set,
                            test_set, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, hold_set_size, holdout.adjustedP, holdin, nmse, onestep, holdout.onestep)
{
  tspX <- tsp(train_set)
  firstTspX <- tsp(main_set)
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
      is.season <- ATA.Seasonality(train_set, s.frequency, seas_attr_set)
    }
  }
  if (is.season==TRUE){
    if (seasonal.model!="decomp" & seasonal.type=="M"){
      out.transform <- ATA.Transform(train_set, tMethod = "Box_Cox", tLambda = 0, tShift = 0)  # lambda = 0 for multiplicative model
	    seas_train_set <- forecast::msts(out.transform$trfmX, start = start(train_set), seasonal.periods = s.frequency)
      seas.type <- "A"
      seas.lambda <- out.transform$tLambda
	    seas.shift <- out.transform$tShift
      seas.transform <- "Box_Cox"
    }else {
      seas_train_set <- train_set
      seas.lambda <- NULL
      seas.transform <- NULL
      seas.type <- seasonal.type
	    seas.shift <- 0
    }
  }else {
    seas_train_set <- train_set
    seasonal.model <- "none"
    seas.type <- seasonal.type <- "A"
    seas.lambda <- NULL
	  seas.shift <- 0
	  seas.transform <- NULL
  }
  ata.seasonal.component <- ATA.Decomposition(seas_train_set, s.model=seasonal.model, s.type=seas.type, s.frequency=s.frequency, seas_attr_set=seas_attr_set)
  seasadj_train_set <- ATA.BackTransform(X=ata.seasonal.component$AdjustedX, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  SeasonalIndex <- ATA.BackTransform(X=ata.seasonal.component$SeasIndex, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  SeasonalActual <- ATA.BackTransform(X=ata.seasonal.component$SeasActual, tMethod=seas.transform, tLambda=seas.lambda, tShift = seas.shift)
  if (seasonal.model=="x13" | seasonal.model=="x11"){
    seasonal.type <- ata.seasonal.component$SeasType
  }
  ChgX <- ATA.Transform(seasadj_train_set, tMethod = transform.method, tLambda = lambda, tShift = shift, bcMethod = boxcox_attr_set$bcMethod, bcLower = boxcox_attr_set$bcLower, bcUpper = boxcox_attr_set$bcUpper)
  seasadj_train_set <- ChgX$trfmX
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
    holdout_part <- ifelse(hold_set_size > 0 & hold_set_size < 1, floor(length(seasadj_train_set) * hold_set_size), hold_set_size)
    valid_len <- length(seasadj_train_set) - holdout_part
    train_len <- length(seasadj_train_set)
    train_set_mat <- forecast::msts(seasadj_train_set[1:valid_len], start = start(train_set), seasonal.periods = s.frequency)
    validation_set <- forecast::msts(seasadj_train_set[(valid_len+1):train_len], start = end(train_set_mat) - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1), seasonal.periods = s.frequency)
  }else {
    train_set_mat <- seasadj_train_set
    validation_set <- NA
  }
  ata.output <- SubATA.Damped(train_set_mat, pb = parP, qb = parQ, model.Type = model.type, accuracy.Type = accuracy.type, level.fix = level.fixed, trend.fix = trend.fixed,
                              trend.Search = trend.search, phiStart = start.phi, phiEnd = end.phi, phiSize = size.phi, initialLevel = initial.level, initialTrend = initial.trend,
                              main_set = seasadj_train_set, Holdout = holdout, HoldoutSet = validation_set, Adjusted_P = holdout.adjustedP, h = h, Holdin = holdin, nmse = nmse,
                              seas_periods = s.frequency, holdout_onestep = holdout.onestep)
  ata.output$h <- h
  if (onestep == FALSE){
    ata.output <- SubATA.Forecast(ata.output, hh=h)
  }else {
    ata.output <- SubATA.OneStepForecast(ata.output, test_set, hh=h)
  }
  ata.output$actual <- main_set
  fit_ata <- ATA.BackTransform(X=ata.output$fitted, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
  forecast_ata <- ATA.BackTransform(X=ata.output$forecast, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
  ata.output$level <- forecast::msts(ATA.BackTransform(X=ata.output$level, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))),
                                     start = start(main_set), seasonal.periods = s.frequency)
  ata.output$trend <- forecast::msts(ATA.BackTransform(X=ata.output$trend, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))),
                                     start = start(main_set), seasonal.periods = s.frequency)
  if(seasonal.type == "A"){
    ATA.fitted <- fit_ata + SeasonalActual
    ATA.forecast <- forecast_ata + OS_SIValue
    if (holdout == TRUE){
      holdout.ata <- ATA.BackTransform(X=ata.output$holdout.forecast, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
      ata.output$holdout.forecast <- forecast::msts(holdout.ata + SeasonalActual[(valid_len+1):train_len], start = end(train_set_mat) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
    }
  }else {
    ATA.fitted <- fit_ata * SeasonalActual
    ATA.forecast <- forecast_ata * OS_SIValue
    if (holdout == TRUE){
      holdout.ata <- ATA.BackTransform(X=ata.output$holdout.forecast, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals)))
      ata.output$holdout.forecast <- forecast::msts(holdout.ata * SeasonalActual[(valid_len+1):train_len], start = start(train_set_mat) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
    }
  }
  ata.output$fitted <- forecast::msts(ATA.fitted, start = start(main_set), seasonal.periods = s.frequency)
  if (negative.forecast==TRUE){
    ata.output$forecast <- forecast::msts(ATA.forecast, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    ata.output$forecast <- forecast::msts(ATA.forecast, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
  }
  ata.output$residuals <- ata.output$actual - ata.output$fitted
  if (holdout == TRUE){
	  ata.output$holdout.training <- forecast::msts(ata.output$actual[1:valid_len], start = start(main_set), seasonal.periods = s.frequency)
    ata.output$holdout.validation <- forecast::msts(ata.output$actual[(valid_len+1):train_len], start = end(ata.output$holdout.training) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
  }
  my_list <- ata.output
  my_list$out.sample <- forecast::msts(test_set, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = s.frequency)
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
  my_list$seasonal.period <- s.frequency
  my_list$seasonal.index <- SeasonalIndex
  my_list$seasonal <- forecast::msts(SeasonalActual, start = start(main_set), seasonal.periods = s.frequency)
  my_list$seasonal.adjusted <- forecast::msts(ATA.BackTransform(X=seasadj_train_set, tMethod=transform.method, tLambda=lambda, tShift=shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ata.output$residuals))),
                                              start = start(main_set), seasonal.periods = s.frequency)
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
  my_list$par.specs <- list("p" = my_list$p, "q" = my_list$q, "phi" = my_list$phi,
                              "trend" = trend_mthd,
                              "seasonal" = seas_mthd,
                              "period" = s.frequency,
                              "decomp_model" = ifelse(seas_mthd == "N", NA, my_list$seasonal.model),
                              "initial_level" = ifelse(my_list$initial.level=="none", NA, TRUE),
                              "initial_trend" = ifelse(my_list$initial.trend=="none", NA, TRUE))
  accuracy_ata <- ATA.Accuracy(my_list, test_set, print.out = FALSE)
  my_list$accuracy <- accuracy_ata
  my_list$onestep <- onestep
  return(my_list)
}
