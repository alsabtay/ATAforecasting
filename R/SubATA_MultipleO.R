#' @importFrom forecast msts
#' @importFrom stats frequency ts tsp tsp<- var
SubATA.MultipleO <- function(ts_input, pb, qb, model.type, seasonal.Test, seasonal.Model, seasonal.Type, seasonal.Frequency, h, accuracy.Type,
                              level.Fix, trend.Fix, trend.Search, phiStart, phiEnd, phiSize, initialLevel, initialTrend, transform.Method, Lambda, Shift, orig_X,
                              OutSample, seas_attr_set, freqYh, ci.Level, negative.Forecast, boxcox_attr_set, Holdout, partition_h, Adjusted_P, Holdin, nmse)
{
  tspX <- tsp(ts_input)
  firstTspX <- tsp(orig_X)
  if (is.null(seasonal.Test)){
    is.season <- ATA.Seasonality(ts_input, seasonal.Frequency, seas_attr_set)
  }else if (seasonal.Test==TRUE){
    is.season <- ATA.Seasonality(ts_input, seasonal.Frequency, seas_attr_set)
  }else {
    if (max(seasonal.Frequency)==1){
      is.season <- FALSE
      seasonal.Type <- "A"
    }else {
      is.season <- TRUE
    }
  }
  if (is.null(seasonal.Model)){
    if (is.season==TRUE & length(seasonal.Frequency)>1){
      seas.model <- c("stl","stR","tbats")
    }else if (is.season==TRUE & length(seasonal.Frequency)==1){
      seas.model <- c("decomp","stl","stR","tbats")
    }else {
      if (is.season==FALSE | min(seasonal.Frequency)==1){
        seas.model <- "none"
        seasonal.Type <- "A"
      }else {
        if (seasonal.Frequency!=12){
          seas.model <- c("decomp","stl", "stlplus", "stR", "tbats")
        }else {
          seas.model <- c("decomp","stl", "stlplus", "stR", "tbats", "x13", "x11")
        }
      }
    }
  }else {
    if (is.season==FALSE | min(seasonal.Frequency)==1){
      seas.model <- "none"
      seasonal.Type <- "A"
    }else if (length(seasonal.Frequency)>1){
      seas.model <- seasonal.Model[!(seasonal.Model %in% "decomp")]
    }else {
      seas.model <- seasonal.Model
    }
  }
  ifelse(is.null(seasonal.Type), seas.type <- c("A","M"), seas.type <- seasonal.Type)
  model.Type <- ifelse(is.null(model.type), "B", model.type)
  max_smo <- length(seas.model)
  if (length(seas.type)==1){
    max_st <- 1
  }else {
    max_st <- 2
  }
  if (is.season==TRUE){
    DeSeas <- rep(NA,length(ts_input))
    DeSI <- rep(NA,max(seasonal.Frequency))
    DeSA <- rep(NA,length(ts_input))
    TA_0 <- rep(NA,length(ts_input))
    TM_0 <- rep(NA,length(ts_input))
    typeName <- as.data.frame("omit")
    for (smo in 1:max_smo){
      for (st in 1:max_st){
        if (seas.model[smo]!="none"){
          org.seas.Type <- seas.type[st]
          if (seas.model[smo]!="decomp" & seas.type[st]=="M"){
            out.transform <- ATA.Transform(ts_input, tMethod = "Box_Cox", tLambda = 0, tShift = 0)  # lambda = 0 for multiplicative model
            X <- out.transform$trfmX
			      seas.Type <- "A"
            seas.Model <- seas.model[smo]
            seas.Lambda <- out.transform$tLambda
			      seas.Shift <- out.transform$tShift
            seas.Transform <- "Box_Cox"
          }else {
            X <- ts_input
            seas.Type <- seas.type[st]
            seas.Model <- seas.model[smo]
            seas.Lambda <- NULL
			      seas.Shift <- 0
            seas.Transform <- NULL
          }
        }else {
          X <- ts_input
          seas.Type <- "A"
          seas.Model <- "none"
          seas.Lambda <- NULL
		      seas.Shift <- 0
          seas.Transform <- NULL
        }
        ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
        AdjX <- ATA.BackTransform(X=ata.seasonal.component$AdjustedX, tMethod=seas.Transform, tLambda=seas.Lambda, tShift = seas.Shift)
        AdjSI <- ATA.BackTransform(X=ata.seasonal.component$SeasIndex, tMethod=seas.Transform, tLambda=seas.Lambda, tShift = seas.Shift)
        AdjSA <- ATA.BackTransform(X=ata.seasonal.component$SeasActual, tMethod=seas.Transform, tLambda=seas.Lambda, tShift = seas.Shift)
        if (seas.Model=="x13" | seas.Model=="x11"){
          seas.Type <- ata.seasonal.component$SeasType
        }
        ChgX <- ATA.Transform(AdjX, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, bcMethod = boxcox_attr_set$bcMethod, bcLower = boxcox_attr_set$bcLower, bcUpper = boxcox_attr_set$bcUpper)
        AdjX <- ChgX$trfmX
        Lambda <- ChgX$tLambda
        Shift <- ChgX$tShift
        DeSeas <- as.matrix.data.frame(cbind(DeSeas, as.numeric(AdjX)))
        DeSI <- as.matrix.data.frame(cbind(DeSI, as.numeric(AdjSI)))
        DeSA <- as.matrix.data.frame(cbind(DeSA, as.numeric(AdjSA)))
        typeName <- cbind(typeName, seas.Type)
        TA_0 <- cbind(TA_0, as.double(AdjX-ATA.Shift(AdjX,1)))
        TM_0 <- cbind(TM_0, as.double(AdjX/ATA.Shift(AdjX,1)))
      }
    }
    orig_DeSeas <- DeSeas <- DeSeas[,-1]
    DeSI <- DeSI[,-1]
    DeSA <- DeSA[,-1]
    TA_0 <- TA_0[,-1]
    TM_0 <- TM_0[,-1]
    if (Holdout == TRUE){
      holdout_part <- ifelse(partition_h > 0 & partition_h < 1, floor(length(ts_input) * partition_h), partition_h)
      HoldOutLen <- length(ts_input) - holdout_part
      InsampleLen <- length(ts_input)
      HoldoutSet <- ts(DeSeas[(HoldOutLen+1):InsampleLen,], frequency = tspX[3], start = tspX[2] - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1))
      DeSeas <- ts(DeSeas[1:HoldOutLen,], frequency = tspX[3], start = tspX[1])
      output <- SubATAHoldout(as.matrix.data.frame(DeSeas)
                               , as.integer(ifelse(pb=="opt", -1, pb))
                               , as.integer(ifelse(qb=="opt", -1, qb))
                               , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                               , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16))
                               , as.integer(ifelse(level.Fix, 1, 0))
                               , as.integer(ifelse(trend.Fix, 1, 0))
                               , as.integer(ifelse(trend.Search, 1, 0))
                               , as.double(phiStart)
                               , as.double(phiEnd)
                               , as.double(phiSize)
                               , as.integer(ifelse(initialLevel, 1, 0))
                               , as.integer(ifelse(initialTrend, 1, 0))
                               , as.matrix.data.frame(TA_0)
                               , as.matrix.data.frame(TM_0)
                               , as.integer(sapply(seas.model, switch, "none"=0,"decomp"=1,"stl"=2,"stlplus"=3,"stR"=4,"tbats"=5,"x13"=6,"x11"=7))
                               , as.integer(sapply(seas.type, switch, "A"=0,"M"=1))
                               , as.integer(max_smo)
                               , as.integer(max_st)
                               , as.integer(frequency(ts_input))
                               , as.matrix.data.frame(HoldoutSet))

    }else if (Holdin == TRUE){
      HoldoutSet <- NA
      output <- SubATAHoldhin(as.matrix.data.frame(DeSeas)
                               , as.integer(ifelse(pb=="opt", -1, pb))
                               , as.integer(ifelse(qb=="opt", -1, qb))
                               , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                               , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16))
                               , as.integer(ifelse(level.Fix, 1, 0))
                               , as.integer(ifelse(trend.Fix, 1, 0))
                               , as.integer(ifelse(trend.Search, 1, 0))
                               , as.double(phiStart)
                               , as.double(phiEnd)
                               , as.double(phiSize)
                               , as.integer(ifelse(initialLevel, 1, 0))
                               , as.integer(ifelse(initialTrend, 1, 0))
                               , as.matrix.data.frame(TA_0)
                               , as.matrix.data.frame(TM_0)
                               , as.integer(sapply(seas.model, switch, "none"=0,"decomp"=1,"stl"=2,"stlplus"=3,"stR"=4,"tbats"=5,"x13"=6,"x11"=7))
                               , as.integer(sapply(seas.type, switch, "A"=0,"M"=1))
                               , as.integer(max_smo)
                               , as.integer(max_st)
                               , as.integer(frequency(ts_input))
                               , as.integer(h)
                               , as.integer(nmse))
    }else {
      HoldoutSet <- NA
      output <- SubATA(as.matrix.data.frame(DeSeas)
                        , as.integer(ifelse(pb=="opt", -1, pb))
                        , as.integer(ifelse(qb=="opt", -1, qb))
                        , as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
                        , as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13,"AMSE"=14,"lik"=15,"sigma"=16))
                        , as.integer(ifelse(level.Fix, 1, 0))
                        , as.integer(ifelse(trend.Fix, 1, 0))
                        , as.integer(ifelse(trend.Search, 1, 0))
                        , as.double(phiStart)
                        , as.double(phiEnd)
                        , as.double(phiSize)
                        , as.integer(ifelse(initialLevel, 1, 0))
                        , as.integer(ifelse(initialTrend, 1, 0))
                        , as.matrix.data.frame(TA_0)
                        , as.matrix.data.frame(TM_0)
                        , as.integer(sapply(seas.model, switch, "none"=0,"decomp"=1,"stl"=2,"stlplus"=3,"stR"=4,"tbats"=5,"x13"=6,"x11"=7))
                        , as.integer(sapply(seas.type, switch, "A"=0,"M"=1))
                        , as.integer(max_smo)
                        , as.integer(max_st)
                        , as.integer(frequency(ts_input))
                        , as.integer(nmse))
    }
    #output[1] = d_opt_p
    #output[2] = d_opt_q
    #output[3] = d_opt_phi
    #output[4] = d_opt_mo
    #output[5] = LastIXSMO
    #output[6] = LastIXST
    #output[7] = mod_clmn
    #output[8] = holdout.accuracy

    AdjInput <- forecast::msts(as.numeric(orig_DeSeas[,output[7]]), start=firstTspX[1], seasonal.periods = seasonal.Frequency)
    SeasonalActual <- forecast::msts(as.numeric(DeSA[,output[7]]), start=firstTspX[1], seasonal.periods = seasonal.Frequency)
    SeasonalIndex <- as.numeric(DeSI[,output[7]])
    if (is.season==FALSE & output[6]==0){
      OS_SIValue <- rep(0,times=h)
    }else if (is.season==FALSE & output[6]==1){
      OS_SIValue <- rep(1,times=h)
    }else if (is.season==TRUE){
      OS_SIValue <- rep(NA,times=h)
      for (k in 1:h){
        OS_SIValue[k] <- SeasonalIndex[freqYh[k]]
      }
    }else{
    }
    ifelse(Holdout==TRUE & Adjusted_P==TRUE, new_pk <- round((output[1] * length(ts_input))/ length(DeSeas[,output[7]])), new_pk <- output[1])
    ATA.last <- ATA.Core(AdjInput, pk = new_pk, qk = output[2], phik = output[3], mdlType = ifelse(output[4]==1,"A","M"), initialLevel = initialLevel, initialTrend = initialTrend, nmse = nmse)
    ATA.last$holdout <- Holdout
    ATA.last$holdin <- Holdin
    if(Holdout==TRUE){
      ATA.last$holdout.accuracy <- output[8]
      ATA.last$holdout.forecast <- ATAHoldoutForecast(as.double(DeSeas[,output[7]])
                                                      , as.integer(output[1])
                                                      , as.integer(output[2])
                                                      , as.double(output[3])
                                                      , as.integer(output[4])
                                                      , as.integer(ifelse(initialLevel, 1, 0))
                                                      , as.integer(ifelse(initialTrend, 1, 0))
                                                      , as.double(TA_0)
                                                      , as.double(TM_0)
                                                      , as.integer(frequency(ts_input))
                                                      , as.integer(length(HoldoutSet)))
    }
  }else {
    X <- ts_input
    seas.Type <- "A"
    OS_SIValue <- rep(0,times=h)
    seas.Model <- "none"
    seas.Lambda <- NULL
	  seas.Shift <- 0
    seas.Transform <- NULL
    ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
    SeasonalActual <- ata.seasonal.component$SeasActual
    SeasonalIndex <- ata.seasonal.component$SeasIndex
    ChgX <- ATA.Transform(ata.seasonal.component$AdjustedX, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, bcMethod = boxcox_attr_set$bcMethod, bcLower = boxcox_attr_set$bcLower, bcUpper = boxcox_attr_set$bcUpper)
    X <- AdjInput <- ChgX$trfmX
    Lambda <- ChgX$tLambda
	  Shift <- ChgX$tShift
    if (Holdout == TRUE){
      holdout_part <- ifelse(partition_h > 0 & partition_h < 1, floor(length(ts_input) * partition_h), partition_h)
      HoldOutLen <- length(ts_input) - holdout_part
      InsampleLen <- length(ts_input)
      HoldoutSet <- ts(X[(HoldOutLen+1):InsampleLen], frequency = tspX[3], start = tspX[2] - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1))
      X <- ts(X[1:HoldOutLen], frequency = tspX[3], start = tspX[1])
    }else {
      HoldoutSet <- NA
    }
    ATA.last <- SubATA.Damped(X, pb = pb, qb = qb, model.Type = model.Type, accuracy.Type = accuracy.Type, level.fix = level.Fix, trend.fix = trend.Fix, trend.Search = trend.Search, phiStart = phiStart, phiEnd = phiEnd, phiSize = phiSize, initialLevel = initialLevel, initialTrend = initialTrend, orig_X = ts_input, Holdout = Holdout, HoldoutSet = HoldoutSet, Adjusted_P = Adjusted_P, h = h, Holdin = Holdin, nmse = nmse)
  }
  ATA.last$h <- h
  ATA.last <- SubATA.Forecast(ATA.last, hh=h, initialLevel = initialLevel)
  ATA.last$actual <- orig_X
  fit.ata <- ATA.BackTransform(X=ATA.last$fitted, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
  forecast.ata <- ATA.BackTransform(X=ATA.last$forecast, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
  fc.amse <- ATA.BackTransform(X=ATA.last$amse.fc, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
  ATA.last$level <- ts(ATA.BackTransform(X=ATA.last$level, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))), frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  ATA.last$trend <- ts(ATA.BackTransform(X=ATA.last$trend, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))), frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  houldout.ata <- ATA.BackTransform(X = ATA.last$holdout.forecast, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
  crit_a <- ifelse(is.season==TRUE, ifelse(output[6]==0,"A","M"), seas.Type)
  crit_a <- ifelse(is.season==FALSE, "A", crit_a)
  if(crit_a=="A"){
    ATA.fitted <- fit.ata + SeasonalActual
    ATA.forecast <- forecast.ata + OS_SIValue
    if (Holdout == TRUE){
      ATA.last$holdout.forecast <- ts(houldout.ata + SeasonalActual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
    }
  }else {
    ATA.fitted <- fit.ata * SeasonalActual
    ATA.forecast <- forecast.ata * OS_SIValue
    if (Holdout == TRUE){
      ATA.last$holdout.forecast <- ts(houldout.ata * SeasonalActual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
    }
  }
  ATA.last$fitted <- ts(ATA.fitted, frequency = firstTspX[3], start =firstTspX[1])
  if (negative.Forecast==TRUE){
    ATA.last$forecast <- ts(ATA.forecast, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    ATA.last$forecast <- ts(ATA.forecast, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }
  ATA.last$residuals <- ATA.last$actual - ATA.last$fitted
  if (Holdout == TRUE){
    ATA.last$holdout.training <- ts(ATA.last$actual[1:HoldOutLen], frequency = firstTspX[3], start = firstTspX[1])
    ATA.last$holdout.validation <- ts(ATA.last$actual[(HoldOutLen+1):InsampleLen], frequency = firstTspX[3], start = tsp(ATA.last$holdout.training)[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  }
  my_list <- ATA.last
  my_list$out.sample <- ts(OutSample, frequency = firstTspX[3], start = firstTspX[2] + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1))
  if (level.Fix==TRUE){
    method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
  }else if (trend.Fix==TRUE){
    method <- paste("ATA(", my_list$p, ",1," ,my_list$phi, ")", sep="")
  }else if (trend.Search==TRUE){
    method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
  }else {
    method <- paste("ATA(", my_list$p, "," ,my_list$q, ",", my_list$phi, ")", sep="")
  }
  my_list$initial.level <- initialLevel
  my_list$initial.trend <- initialTrend
  my_list$level.fixed <- level.Fix
  my_list$trend.fixed <- trend.Fix
  my_list$trend.search <- trend.Search
  my_list$transform.method <- transform.Method
  my_list$lambda <- Lambda
  my_list$shift <- Shift
  my_list$bcLower <- boxcox_attr_set$bcLower
  my_list$bcUpper <- boxcox_attr_set$bcUpper
  my_list$bcBiasAdj <- boxcox_attr_set$bcBiasAdj
  my_list$accuracy.type <- accuracy.Type
  my_list$nmse <- nmse
  my_list$is.season <- is.season
  my_list$seasonal.model <- ifelse(is.season==TRUE, switch(output[5]+1, "none", "decomp", "stl", "stlplus", "stR", "tbats", "x13", "x11"), "none")
  if  (!is.null(seasonal.Type)){
    my_list$seasonal.type <- seasonal.Type
  }else {
    my_list$seasonal.type <- crit_a
  }
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
  my_list$seasonal.period <- seasonal.Frequency
  my_list$seasonal.index <- SeasonalIndex
  my_list$seasonal <- ts(SeasonalActual, frequency = firstTspX[3], start=firstTspX[1])
  my_list$seasonal.adjusted <- ts(ATA.BackTransform(X=AdjInput, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))), frequency = firstTspX[3], start=firstTspX[1])
  ci.output <- ATA.CI(object = my_list, ci.level = ci.Level)
  my_list$ci.level <- ci.Level
  if (negative.Forecast==TRUE){
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
