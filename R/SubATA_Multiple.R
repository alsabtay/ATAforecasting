#' @importFrom forecast msts
#' @importFrom stats end start ts tsp tsp<- var
SubATA_Multi_Before <- function(train_set, pb, qb, model.type, seasonal.Test, seasonal.Model, seasonal.Type, seasonal.Frequency, h, accuracy.Type,
                                level.Fix, trend.Fix, trend.Search, phiStart, phiEnd, phiSize, initialLevel, initialTrend, transform.Method, Lambda, Shift, main_set,
                                test_set, seas_attr_set, freqYh, ci.Level, negative.Forecast, boxcox_attr_set, Holdout, hold_set_size, Adjusted_P, Holdin, nmse)
{
  tspX <- tsp(train_set)
  firstTspX <- tsp(main_set)
  if (is.null(seasonal.Test)){
    is.season <- ATA.Seasonality(train_set, seasonal.Frequency, seas_attr_set)
  }else if (seasonal.Test==TRUE){
    is.season <- ATA.Seasonality(train_set, seasonal.Frequency, seas_attr_set)
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
      if (is.season==FALSE | max(seasonal.Frequency)==1){
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
    if (is.season==FALSE | max(seasonal.Frequency)==1){
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
    train_set_mat <- rep(NA,length(train_set))
    DeSI <- rep(NA,max(seasonal.Frequency))
    DeSA <- rep(NA,length(train_set))
    TA_0 <- rep(NA,length(train_set))
    TM_0 <- rep(NA,length(train_set))
    typeName <- as.data.frame("omit")
    for (smo in 1:max_smo){
      for (st in 1:max_st){
        if (seas.model[smo]!="none"){
          org.seas.Type <- seas.type[st]
          if (seas.model[smo]!="decomp" & seas.type[st]=="M"){
            out.transform <- ATA.Transform(train_set, tMethod = "Box_Cox", tLambda = 0, tShift = 0)  # lambda = 0 for multiplicative model
			      seas_train_set <- forecast::msts(out.transform$trfmX, start = start(train_set), seasonal.periods = seasonal.Frequency)
            seas.Type <- "A"
            seas.Model <- seas.model[smo]
            seas.Lambda <- out.transform$tLambda
			      seas.Shift <- out.transform$tShift
            seas.Transform <- "Box_Cox"
          }else {
            seas_train_set <- train_set
            seas.Type <- seas.type[st]
            seas.Model <- seas.model[smo]
            seas.Lambda <- NULL
			      seas.Shift <- 0
            seas.Transform <- NULL
          }
        }else {
          seas_train_set <- train_set
          seas.Type <- "A"
          seas.Model <- "none"
          seas.Lambda <- NULL
    		  seas.Shift <- 0
          seas.Transform <- NULL
        }
        ata.seasonal.component <- ATA.Decomposition(seas_train_set, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
        seasadj_train_set <- ATA.BackTransform(X=ata.seasonal.component$AdjustedX, tMethod=seas.Transform, tLambda=seas.Lambda, tShift=seas.Shift)
        AdjSI <- ATA.BackTransform(X=ata.seasonal.component$SeasIndex, tMethod=seas.Transform, tLambda=seas.Lambda, tShift=seas.Shift)
        AdjSA <- ATA.BackTransform(X=ata.seasonal.component$SeasActual, tMethod=seas.Transform, tLambda=seas.Lambda, tShift=seas.Shift)
        if (seas.Model=="x13" | seas.Model=="x11"){
          seas.Type <- ata.seasonal.component$SeasType
        }
        train_set_mat <- as.matrix.data.frame(cbind(train_set_mat, as.numeric(seasadj_train_set)))
        DeSI <- as.matrix.data.frame(cbind(DeSI, as.numeric(AdjSI)))
        DeSA <- as.matrix.data.frame(cbind(DeSA, as.numeric(AdjSA)))
        typeName <- cbind(typeName, seas.Type)
        TA_0 <- cbind(TA_0, as.double(seasadj_train_set - ATA.Shift(seasadj_train_set,1)))
        TM_0 <- cbind(TM_0, as.double(seasadj_train_set / ATA.Shift(seasadj_train_set,1)))
      }
    }
    main_train_set_mat <- train_set_mat <- train_set_mat[,-1]
    DeSI <- DeSI[,-1]
    DeSA <- DeSA[,-1]
    TA_0 <- TA_0[,-1]
    TM_0 <- TM_0[,-1]
    if (Holdout == TRUE){
      holdout_part <- ifelse(hold_set_size > 0 & hold_set_size < 1, floor(length(train_set) * hold_set_size), hold_set_size)
      valid_len <- length(train_set) - holdout_part
      train_len <- length(train_set)
      train_set_mat <- forecast::msts(main_train_set_mat[1:valid_len,], start = start(main_set), seasonal.periods = seasonal.Frequency)
      validation_set <- forecast::msts(main_train_set_mat[(valid_len+1):train_len,], start = end(train_set_mat) - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1), seasonal.periods = seasonal.Frequency)
      output <- SubATAHoldout(as.matrix.data.frame(train_set_mat)
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
                               , as.double(seasonal.Frequency)
                               , as.matrix.data.frame(validation_set))

    }else if (Holdin == TRUE){
      validation_set <- NA
      output <- SubATAHoldhin(as.matrix.data.frame(train_set_mat)
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
                               , as.double(seasonal.Frequency)
                               , as.integer(h)
                               , as.integer(nmse))
    }else {
      validation_set <- NA
      output <- SubATA(as.matrix.data.frame(train_set_mat)
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
                        , as.double(seasonal.Frequency)
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

    AdjInput <- forecast::msts(as.numeric(main_train_set_mat[,output[7]]), start = start(main_set), seasonal.periods = seasonal.Frequency)
    SeasonalActual <- forecast::msts(as.numeric(DeSA[,output[7]]), start = start(main_set), seasonal.periods = seasonal.Frequency)
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
    ifelse(Holdout==TRUE & Adjusted_P==TRUE, new_pk <- round((output[1] * length(train_set))/ length(train_set_mat[,output[7]])), new_pk <- output[1])
    ATA.last <- ATA.Core(AdjInput, pk = new_pk, qk = output[2], phik = output[3], mdlType = ifelse(output[4]==1,"A","M"), initialLevel = initialLevel, initialTrend = initialTrend, nmse = nmse)
    ATA.last$holdout <- Holdout
    ATA.last$holdin <- Holdin
    if(Holdout==TRUE){
      ATA.last$holdout.accuracy <- output[8]
      ATA.last$holdout.forecast <- ATAHoldoutForecast(as.double(train_set_mat[,output[7]])
                                                      , as.integer(output[1])
                                                      , as.integer(output[2])
                                                      , as.double(output[3])
                                                      , as.integer(output[4])
                                                      , as.integer(ifelse(initialLevel, 1, 0))
                                                      , as.integer(ifelse(initialTrend, 1, 0))
                                                      , as.double(TA_0)
                                                      , as.double(TM_0)
                                                      , as.integer(frequency(train_set))
                                                      , as.integer(length(validation_set)))
    }
  }else {
    seas.Type <- "A"
    OS_SIValue <- rep(0,times=h)
    seas.Model <- "none"
    seas.Lambda <- NULL
	  seas.Shift <- 0
    seas.Transform <- NULL
    ata.seasonal.component <- ATA.Decomposition(train_set, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
    SeasonalActual <- ata.seasonal.component$SeasActual
    SeasonalIndex <- ata.seasonal.component$SeasIndex
	  AdjInput <- seasadj_train_set <- ata.seasonal.component$AdjustedX
    if (Holdout == TRUE){
      holdout_part <- ifelse(hold_set_size > 0 & hold_set_size < 1, floor(length(train_set) * hold_set_size), hold_set_size)
      valid_len <- length(train_set) - holdout_part
      train_len <- length(train_set)
      train_set_mat <- forecast::msts(seasadj_train_set[1:valid_len], start = start(train_set), seasonal.periods = seasonal.Frequency)
      validation_set <- forecast::msts(seasadj_train_set[(valid_len+1):train_len], start = end(train_set_mat) - ifelse(tspX[3]>1, (holdout_part - 1) * (1/tspX[3]), (holdout_part - 1) * 1), seasonal.periods = seasonal.Frequency)
    }else {
      train_set_mat <- seasadj_train_set
      validation_set <- NA
    }
    ATA.last <- SubATA.Damped(train_set_mat, pb = pb, qb = qb, model.Type = model.Type, accuracy.Type = accuracy.Type, level.fix = level.Fix, trend.fix = trend.Fix, trend.Search = trend.Search, phiStart = phiStart, phiEnd = phiEnd, phiSize = phiSize,
                              initialLevel = initialLevel, initialTrend = initialTrend, main_set = seasadj_train_set, Holdout = Holdout, HoldoutSet = validation_set, Adjusted_P = Adjusted_P, h = h, Holdin = Holdin, nmse = nmse, seas_periods = seasonal.Frequency)
  }
  ATA.last$h <- h
  ATA.last <- SubATA.Forecast(ATA.last, hh=h, initialLevel = initialLevel)
  ATA.last$actual <- main_set
  fit_ata <- ATA.last$fitted
  forecast_ata <- ATA.last$forecast
  fc.amse <- ATA.last$amse.fc
  ATA.last$level <- forecast::msts(ATA.BackTransform(X=ATA.last$level, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                   start = start(main_set), seasonal.periods = seasonal.Frequency)
  ATA.last$trend <- forecast::msts(ATA.BackTransform(X=ATA.last$trend, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                   start = start(main_set), seasonal.periods = seasonal.Frequency)
  crit_a <- ifelse(is.season==TRUE, ifelse(output[6]==0,"A","M"), seas.Type)
  crit_a <- ifelse(is.season==FALSE, "A", crit_a)
  if(crit_a=="A"){
    ATA.fitted <- ATA.BackTransform(X = fit_ata + SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
    ATA.forecast <- ATA.BackTransform(X = forecast_ata + OS_SIValue, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
    if (Holdout == TRUE){
      ATA.last$holdout.forecast <- forecast::msts(ATA.BackTransform(X = ATA.last$holdout.forecast + SeasonalActual[(valid_len+1):train_len], tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                      start = end(train_set_mat) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
    }
  }else {
    ATA.fitted <- ATA.BackTransform(X = fit_ata * SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
    ATA.forecast <- ATA.BackTransform(X = forecast_ata * OS_SIValue, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
    if (Holdout == TRUE){
      ATA.last$holdout.forecast <- forecast::msts(ATA.BackTransform(X = ATA.last$holdout.forecast * SeasonalActual[(valid_len+1):train_len], tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                      start = end(train_set_mat) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
    }
  }
  ATA.last$fitted <- forecast::msts(ATA.fitted, start = start(main_set), seasonal.periods = seasonal.Frequency)
  if (negative.Forecast==TRUE){
    ATA.last$forecast <- forecast::msts(ATA.forecast, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    ATA.last$forecast <- forecast::msts(ATA.forecast, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
  }
  ATA.last$residuals <- ATA.last$actual - ATA.last$fitted
  if (Holdout == TRUE){
    ATA.last$holdout.training <- forecast::msts(ATA.last$actual[1:valid_len], start = start(main_set), seasonal.periods = seasonal.Frequency)
    ATA.last$holdout.validation <- forecast::msts(ATA.last$actual[(valid_len+1):train_len], start = end(ATA.last$holdout.training) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
  }
  my_list <- ATA.last
  my_list$out.sample <- forecast::msts(test_set, start = end(main_set) + ifelse(firstTspX[3]>1, 1/firstTspX[3], 1), seasonal.periods = seasonal.Frequency)
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
  my_list$seasonal.period <- seasonal.Frequency
  my_list$seasonal.index <- ATA.BackTransform(X=SeasonalIndex, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
  my_list$seasonal <- forecast::msts(ATA.BackTransform(X=SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                     start = start(main_set), seasonal.periods = seasonal.Frequency)
  my_list$seasonal.adjusted <- forecast::msts(ATA.BackTransform(X=AdjInput, tMethod=transform.Method, tLambda=Lambda, tShift=Shift, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals))),
                                              start = start(main_set), seasonal.periods = seasonal.Frequency)
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
  my_list$par.specs <- list("p" = my_list$p, "q" = my_list$q, "phi" = my_list$phi,
                              "trend" = trend_mthd,
                              "seasonal" = seas_mthd,
                              "period" = seasonal.Frequency,
                              "decomp_model" = ifelse(seas_mthd == "N", NA, my_list$seasonal.model),
                              "initial_level" = ifelse(my_list$initial.level==FALSE, NA, TRUE),
                              "initial_trend" = ifelse(my_list$initial.trend==FALSE, NA, TRUE))
  accuracy_ata <- ATA.Accuracy(my_list, test_set, print.out = FALSE)
  my_list$accuracy <- accuracy_ata
  return(my_list)
}
