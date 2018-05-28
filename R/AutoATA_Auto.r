#' @export

AutoATA.Auto <- function(ts_input, pb, qb, model.Type, seasonal.Test, seasonal.Model, seasonal.Type, seasonal.Frequency, h, accuracy.Type, 
							level.Fix, trend.Fix, phiStart, phiEnd, phiSize, initialLevel, initialTrend, transform.Method, Lambda, orig_X, 
							OutSample, seas_attr_set, freqYh, ci.Level, negative.Forecast)
{
	tspX <- tsp(ts_input)
	optAccryStart <- 9999999999999.9
	if (is.null(seasonal.Model)){
		if (length(seasonal.Frequency)>1){
			seas.model <- c("stR","tbats")
		}else {
			if (seasonal.Frequency==1){
				seasonal.Model <- seas.model <- "none"
				seasonal.Test <- TRUE
			}else {
				if (seasonal.Frequency!=12){
					seas.model <- c("decomp","stl", "stlplus", "stR", "tbats")
				}else {
					seas.model <- c("decomp","stl", "stlplus", "stR", "tbats", "x13", "x11")
				}
			}
		}
	}else {
		if (length(seasonal.Frequency)>1){
				seas.model <- c("stR","tbats")
		}else {
			if (seasonal.Frequency==1){
				seasonal.Model <- seas.model <- "none"
				seasonal.Test <- TRUE
			}else {
				seas.model <- seasonal.Model
			}
		}
	}
	ifelse(is.null(seasonal.Type), ifelse(seas.model=="none", seas.type <- "A", seas.type <- c("A","M")), seas.type <- seasonal.Type)
	ifelse(is.null(model.Type),model.Type <- "B",model.Type <- model.Type)
	if (is.null(seasonal.Test)){
		max_st <- 2L
		s_Test <- c(TRUE, FALSE)
	}else if (seasonal.Frequency==1){
		max_st <- 1L
		s_Test <-  TRUE
	}else {
		max_st <- 1L
		s_Test <- seasonal.Test
	}
	max_s <- length(seas.model)
	for (s in 1:max_s){
		if (seas.model[s]=="x13" | seas.model[s]=="x11" | seas.model[s]=="none"){
			max_n <- 1L
		}else {
			if (length(seas.type)==1){
				max_n <- 1L
			}else {
				max_n <- 2L
			}
		}
		for (n in 1:max_n) {
			for (st in 1:max_st) {
				if (seas.model[s]=="none"){
					is.season <- FALSE
				}else {
					if (s_Test[st]==FALSE){
						is.season <- TRUE
					}else {
						is.season <- SeasonalityTest(ts_input, seasonal.Frequency, seas_attr_set)
					}
				}
				if (is.season==TRUE){
					org.seas.Type <- seas.type[n]
					if (seas.model[s]!="decomp" & seas.type[n]=="M"){
						X <- ATA.Transform(ts_input,tMethod="BoxCox",tLambda=0)$trfmX  # lambda = 0 for multiplicative model
						seas.Type <- "A"
						seas.Model <- seas.model[s]
						seas.Lambda <- 0
						seas.Transform <- "BoxCox"
					}else {
						X <- ts_input
						seas.Type <- seas.type[n]
						seas.Model <- seas.model[s]
						seas.Lambda <- NULL
						seas.Transform <- NULL
					}
				}else {
					X <- ts_input
					seas.Type <- "A"
					seas.Model <- "none"
					seas.Lambda <- NULL
					seas.Transform <- NULL
				}
				ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
				AdjInSample <- ATA.Inv.Transform(X=ata.seasonal.component$AdjustedX, tMethod=seas.Transform, tLambda=seas.Lambda)
				SeasonalIndex <- ATA.Inv.Transform(X=ata.seasonal.component$SeasIndex, tMethod=seas.Transform, tLambda=seas.Lambda)
				SeasonalActual <- ATA.Inv.Transform(X=ata.seasonal.component$SeasActual, tMethod=seas.Transform, tLambda=seas.Lambda)
				if (seas.model[s]=="x13" | seas.model[s]=="x11"){
					if (min(abs(SeasonalIndex))>0 & max(abs(SeasonalIndex))<3){
						seas.Type <- "M"
					}else {
						seas.Type <- "A"
					}
					org.seas.Type <- seas.Type
				}
				if (is.season==FALSE & seas.Type=="A"){
					OS_SIValue <- rep(0,times=h)
				}else if (is.season==FALSE & seas.Type=="M"){
					OS_SIValue <- rep(1,times=h)
				}else if (is.season==TRUE){
					OS_SIValue <- rep(NA,times=h)
					for (k in 1:h){
						OS_SIValue[k] <- SeasonalIndex[freqYh[k]]
					}
				}else{
				}
				Xdata <- as.numeric(AdjInSample)
				TA_0 <- Xdata-ATA.Shift(Xdata,1)
				TM_0 <- Xdata/ATA.Shift(Xdata,1)
				output <- AutoATADamped(as.double(Xdata)
								, as.integer(ifelse(pb=="opt", -1, pb))
								, as.integer(ifelse(qb=="opt", -1, qb))
								, as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
								, as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12))
								, as.integer(ifelse(level.Fix, 1, 0))
								, as.integer(ifelse(trend.Fix, 1, 0))
								, as.double(phiStart)
								, as.double(phiEnd)
								, as.double(phiSize)
								, as.integer(ifelse(initialLevel, 1, 0))
								, as.integer(ifelse(initialTrend, 1, 0))
								, as.double(TA_0)
								, as.double(TM_0))
				ATA.opt <- AutoATA.Core(AdjInSample, pk = output[1], qk = output[2], phik = output[3], mdlType = ifelse(output[4]==1,"A","M"), initialLevel = initialLevel, initialTrend = initialTrend)
				optAccryEnd <- AutoATA.Accuracy(ATA.opt, accryType = accuracy.Type)
				if (optAccryEnd <= optAccryStart){
					opt_p <- output[1] 
					opt_q <- output[2] 
					opt_phi <- output[3]
					ifelse(output[4]==1,model.Type <- "A", model.Type <- "M")
					seasonal.Model <- seas.Model
					org.seas.type <- org.seas.Type
					seasonal.Type <- seas.Type
					seasonal.Adj <- AdjInSample
					seasonal.Index <- SeasonalIndex
					seasonal.Actual <- SeasonalActual
					seasonal.forecast <- OS_SIValue
					opt_lambda <- Lambda
					opt_transform <- transform.Method
					optAccryStart <- optAccryEnd
				}
			}
		}
	}
	ATA.last <- ATA.Core(seasonal.Adj, pk = opt_p, qk = opt_q, phik = opt_phi, mdlType = model.Type, initialLevel = initialLevel, initialTrend = initialTrend)
	ATA.last$h <- h
	ATA.last <- AutoATA.Forecast(ATA.last, hh=h, initialLevel = initialLevel)
	ATA.last$actual <- msts(orig_X, start=tsp(orig_X)[1], seasonal.periods = seasonal.Frequency)
	fit.ata <- ATA.last$fitted
	forecast.ata <- ATA.last$forecast
	ATA.last$level <- ATA.Inv.Transform(X=ATA.last$level, tMethod=opt_transform, tLambda=opt_lambda)
	ATA.last$trend <- ATA.Inv.Transform(X=ATA.last$trend, tMethod=opt_transform, tLambda=opt_lambda)
	if(seasonal.Type=="A"){
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata + seasonal.Actual, tMethod=opt_transform, tLambda=opt_lambda)
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata + seasonal.forecast, tMethod=opt_transform, tLambda=opt_lambda)
	}else {
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata * seasonal.Actual, tMethod=opt_transform, tLambda=opt_lambda)
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata * seasonal.forecast, tMethod=opt_transform, tLambda=opt_lambda)				
	}
	seasonal.Actual <- ATA.Inv.Transform(X=seasonal.Actual, tMethod=opt_transform, tLambda=opt_lambda)
	ATA.last$fitted <- ATA.fitted
	if (negative.Forecast==TRUE){
		ATA.last$forecast <- ATA.forecast
	}else {
		ATA.forecast[ATA.forecast<0] <- 0
		ATA.last$forecast <- ATA.forecast
	}	
	accuracy.ata <- ATA.Accuracy(ATA.last, OutSample)	
	my_list <- ATA.last
	my_list$out.sample <- OutSample
	if (level.Fix==TRUE){
		method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
	}else if (trend.Fix==TRUE){
			method <- paste("ATA(", my_list$p, ",1," ,my_list$phi, ")", sep="")
	}else {
		method <- paste("ATA(", my_list$p, "," ,my_list$q, ",", my_list$phi, ")", sep="")
	}
	my_list$method <- method
	my_list$initial.value <- initialLevel
	my_list$level.fixed <- level.Fix
	my_list$trend.fixed <- trend.Fix
	my_list$transform.method <- opt_transform
	my_list$lambda <- opt_lambda
	my_list$accuracy.type <- accuracy.Type
	my_list$accuracy <- accuracy.ata
	my_list$is.season <- is.season
	my_list$seasonal.model <- seasonal.Model
	my_list$seasonal.type <- org.seas.type
	my_list$seasonal.period <- seasonal.Frequency
	my_list$seasonal.index <- seasonal.Index
	my_list$seasonal <- seasonal.Actual
	my_list$seasonal.adjusted <- ATA.Inv.Transform(X=seasonal.Adj, tMethod=opt_transform, tLambda=opt_lambda)
	ci.output <- ATA.CI(my_list, ci.Level)
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