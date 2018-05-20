#' @export

AutoATA.Auto <- function(ts_input, pb, qb, model.Type, seasonal.Test, seasonal.Model, seasonal.Type, seasonal.Frequency, h, accuracy.Type, 
							level.Fix, trend.Fix, phiStart, phiEnd, phiSize, initialLevel, initialTrend, transform.Method, Lambda, orig_X, 
							OutSample, seas_attr_set, freqYh, ci.Level)
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
				seas.model <- c("decomp","stl", "stlplus", "stR", "tbats", "x13", "x11")
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
			max_n <- 2L
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
					if (seas.model[s]=="x13" | seas.model[s]=="x11"){
						X <- ts_input
						seas.Type <- seas.type[n]
						seas.Model <- seas.model[s]
						Lambda <- NULL
						transform.Method <- NULL
					}else if (seas.model[s]!="decomp" & seas.type[n]=="M"){
						X <- ATA.Transform(ts_input,tMethod="BoxCox",tLambda=0)$trfmX  # lambda = 0 for multiplicative model
						seas.Type <- "A"
						seas.Model <- seas.model[s]
						model.Type <- "A"
						Lambda <- 0
						transform.Method <- "BoxCox"
					}else {
						tX <- ATA.Transform(ts_input,tMethod=transform.Method,tLambda=Lambda)
						X <- tX$trfmX
						Lambda <- tX$tLambda
						seas.Type <- seas.type[n]
						seas.Model <- seas.model[s]
					}
				}else {
					tX <- ATA.Transform(ts_input,tMethod=transform.Method,tLambda=Lambda)
					X <- tX$trfmX
					Lambda <- tX$tLambda
					seas.Type <- "A"
					seas.Model <- "none"
				}
				ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
				AdjInSample <- ata.seasonal.component$AdjustedX
				SeasonalIndex <- ata.seasonal.component$SeasIndex
				SeasonalActual <- ata.seasonal.component$SeasActual
				if (seas.model[s]=="x13" | seas.model[s]=="x11"){
					if (min(SeasonalIndex)>0 & max(SeasonalIndex)<3){
						seas.Type <- "M"
					}else {
						seas.Type <- "A"
					}
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
					seasonal.Type <- seas.Type
					seasonal.Adj <- AdjInSample
					seasonal.Index <- SeasonalIndex
					seasonal.Actual <- SeasonalActual
					seasonal.forecast <- OS_SIValue
					optAccryStart <- optAccryEnd
				}
			}
		}
	}
	ATA.last <- ATA.Core(seasonal.Adj, pk = opt_p, qk = opt_q, phik = opt_phi, mdlType = model.Type, initialLevel = initialLevel, initialTrend = initialTrend)
	ATA.last$h <- h
	ATA.last <- AutoATA.Forecast(ATA.last, hh=h, initialLevel = initialLevel)
	ATA.last$actual <- orig_X
	fit.ata <- ATA.last$fitted
	forecast.ata <- ATA.last$forecast
	ATA.last$level <- ATA.Inv.Transform(X=ATA.last$level, tMethod=transform.Method, tLambda=Lambda)
	ATA.last$trend <- ATA.Inv.Transform(X=ATA.last$trend, tMethod=transform.Method, tLambda=Lambda)
	if(seasonal.Type=="A"){
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata + seasonal.Actual, tMethod=transform.Method, tLambda=Lambda)
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata + seasonal.forecast, tMethod=transform.Method, tLambda=Lambda)
	}else {
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata * seasonal.Actual, tMethod=transform.Method, tLambda=Lambda)
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata * seasonal.forecast, tMethod=transform.Method, tLambda=Lambda)				
	}
	seasonal.Actual <- ATA.Inv.Transform(X=seasonal.Actual, tMethod=transform.Method, tLambda=Lambda)
	ATA.last$fitted <- ATA.fitted
	ATA.last$forecast <- ATA.forecast
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
	my_list$transform.method <- transform.Method
	my_list$lambda <- Lambda
	my_list$accuracy.type <- accuracy.Type
	my_list$accuracy <- accuracy.ata
	my_list$is.season <- is.season
	my_list$seasonal.model <- seasonal.Model
	my_list$seasonal.type <- seasonal.Type
	my_list$seasonal.period <- seasonal.Frequency
	my_list$seasonal.index <- seasonal.Index
	my_list$seasonal <- seasonal.Actual
	my_list$seasonal.adjusted <- ATA.Inv.Transform(X=seasonal.Adj, tMethod=transform.Method, tLambda=Lambda)
	ci.output <- ATA.CI(my_list, ci.Level)
	my_list$ci.level <- ci.Level
	my_list$forecast.lower <- ci.output$forecast.lower
	my_list$forecast.upper <- ci.output$forecast.upper
	return(my_list)
}