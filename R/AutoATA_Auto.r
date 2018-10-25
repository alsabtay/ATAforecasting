#' @export AutoATA.Auto

AutoATA.Auto <- function(ts_input, pb, qb, model.Type, seasonal.Test, seasonal.Model, seasonal.Type, seasonal.Frequency, h, accuracy.Type, 
							level.Fix, trend.Fix, phiStart, phiEnd, phiSize, initialLevel, initialTrend, transform.Method, Lambda, orig_X, 
							OutSample, seas_attr_set, freqYh, ci.Level, negative.Forecast, boxcox_attr_set)
{
	tspX <- tsp(ts_input)
	if (is.null(seasonal.Test)){
		is.season <- ATA.Seasonality(ts_input, seasonal.Frequency, seas_attr_set)
	}else if (seasonal.Test==TRUE){
		is.season <- ATA.Seasonality(ts_input, seasonal.Frequency, seas_attr_set)
	}else {
		if (max(seasonal.Frequency)==1){
			is.season <- FALSE
		}else {
			is.season <- TRUE
		}
	}
	if (is.null(seasonal.Model)){
		if (is.season==TRUE & length(seasonal.Frequency)>1){
			seas.model <- c("stl","stR","tbats")
		}else {
			if (is.season==FALSE | seasonal.Frequency==1){
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
		if (is.season==TRUE & length(seasonal.Frequency)>1){
				seas.model <- c("stl","stR","tbats")
		}else {
			if (is.season==FALSE | seasonal.Frequency==1){
				seas.model <- "none"
				seasonal.Type <- "A"
			}else {
				seas.model <- seasonal.Model
			}
		}
	}
	ifelse(is.null(seasonal.Type), seas.type <- c("A","M"), seas.type <- seasonal.Type)
	ifelse(is.null(model.Type),model.Type <- "B",model.Type <- model.Type)
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
						X <- ATA.Transform(ts_input,tMethod="BoxCox",tLambda=0)$trfmX  # lambda = 0 for multiplicative model
						seas.Type <- "A"
						seas.Model <- seas.model[smo]
						seas.Lambda <- 0
						seas.Transform <- "BoxCox"
					}else {
						X <- ts_input
						seas.Type <- seas.type[st]
						seas.Model <- seas.model[smo]
						seas.Lambda <- NULL
						seas.Transform <- NULL
					}
				}else {
					X <- ts_input
					org.seas.Type <- seas.Type <- "A"
					seas.Model <- "none"
					seas.Lambda <- NULL
					seas.Transform <- NULL
				}
				ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
				AdjX <- ATA.Inv.Transform(X=ata.seasonal.component$AdjustedX, tMethod=seas.Transform, tLambda=seas.Lambda)
				AdjSI <- ATA.Inv.Transform(X=ata.seasonal.component$SeasIndex, tMethod=seas.Transform, tLambda=seas.Lambda)
				AdjSA <- ATA.Inv.Transform(X=ata.seasonal.component$SeasActual, tMethod=seas.Transform, tLambda=seas.Lambda)
				if (seas.Model=="x13" | seas.Model=="x11"){
					org.seas.Type <- seas.Type <- ata.seasonal.component$SeasType
				}
				#colSMOName <- colnames(DeSeas)
				DeSeas <- as.matrix.data.frame(cbind(DeSeas, as.numeric(AdjX)))
				DeSI <- as.matrix.data.frame(cbind(DeSI, as.numeric(AdjSI)))
				DeSA <- as.matrix.data.frame(cbind(DeSA, as.numeric(AdjSA)))
				#colnames(DeSeas) <- c(colSMOName,seas.Model)
				#colnames(DeSI) <- colnames(DeSeas)
				#colnames(DeSA) <- colnames(DeSeas)
				#colSTypeName <- colnames(typeName)
				typeName <- cbind(typeName, org.seas.Type)
				#colnames(typeName) <- c(colSTypeName, org.seas.Type)
				TA_0 <- cbind(TA_0, as.double(AdjX-ATA.Shift(AdjX,1)))
				TM_0 <- cbind(TM_0, as.double(AdjX/ATA.Shift(AdjX,1)))
				#colnames(TA_0) <- colnames(DeSeas)
				#colnames(TM_0) <- colnames(DeSeas)
			}
		}
		DeSeas <- DeSeas[,-1]
		DeSI <- DeSI[,-1]
		DeSA <- DeSA[,-1]
		TA_0 <- TA_0[,-1]
		TM_0 <- TM_0[,-1]
		output <- AutoATA(as.matrix.data.frame(DeSeas)
						, as.integer(ifelse(pb=="opt", -1, pb))
						, as.integer(ifelse(qb=="opt", -1, qb))
						, as.integer(switch(model.Type,"B"=0,"A"=1,"M"=2))
						, as.integer(switch(accuracy.Type,"MAE"=1,"MdAE"=2,"MSE"=3,"MdSE"=4,"MPE"=5,"MdPE"=6,"MAPE"=7,"MdAPE"=8,"sMAPE"=9,"sMdAPE"=10,"RMSE"=11,"MASE"=12,"OWA"=13))
						, as.integer(ifelse(level.Fix, 1, 0))
						, as.integer(ifelse(trend.Fix, 1, 0))
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
						, as.integer(length(ts_input))
						, as.integer(frequency(ts_input)))
		
		#output[1] = d_opt_p;
		#output[2] = d_opt_q;
		#output[3] = d_opt_phi;
		#output[4] = d_opt_mo;
		#output[5] = LastIXSMO;
		#output[6] = LastIXST;
		#output[7] = mod_clmn;
		
		AdjInput <- msts(as.numeric(DeSeas[,output[7]]), start=tsp(orig_X)[1], seasonal.periods = seasonal.Frequency)
		SeasonalActual <- msts(as.numeric(DeSA[,output[7]]), start=tsp(orig_X)[1], seasonal.periods = seasonal.Frequency)
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
		ATA.last <- ATA.Core(AdjInput, pk = output[1], qk = output[2], phik = output[3], mdlType = ifelse(output[4]==1,"A","M"), initialLevel = initialLevel, initialTrend = initialTrend)
	}else {
		X <- ts_input
		org.seas.Type <- seas.Type <- "A"
		seas.Model <- "none"
		seas.Lambda <- NULL
		seas.Transform <- NULL
		ata.seasonal.component <- ATA.Decomposition(X, s.model=seas.Model, s.type=seas.Type, s.frequency=seasonal.Frequency, seas_attr_set=seas_attr_set)
		AdjInput <- ATA.Inv.Transform(X=ata.seasonal.component$AdjustedX, tMethod=seas.Transform, tLambda=seas.Lambda)
		SeasonalActual <- ATA.Inv.Transform(X=ata.seasonal.component$SeasActual, tMethod=seas.Transform, tLambda=seas.Lambda)
		SeasonalIndex <- ATA.Inv.Transform(X=ata.seasonal.component$SeasIndex, tMethod=seas.Transform, tLambda=seas.Lambda)
		OS_SIValue <- rep(0,times=h)
		ATA.last <- AutoATA.Damped(AdjInput, pb = pb, qb = qb, model.Type = model.Type, accuracy.Type = accuracy.Type, level.fix = level.Fix, trend.fix = trend.Fix, phiStart = phiStart, phiEnd = phiEnd, phiSize = phiSize, initialLevel = initialLevel, initialTrend = initialTrend)
	}
	ATA.last$h <- h
	ATA.last <- AutoATA.Forecast(ATA.last, hh=h, initialLevel = initialLevel)
	ATA.last$actual <- msts(orig_X, start=tsp(orig_X)[1], seasonal.periods = seasonal.Frequency)
	fit.ata <- ATA.last$fitted
	forecast.ata <- ATA.last$forecast
	ATA.last$level <- ATA.Inv.Transform(X=ATA.last$level, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
	ATA.last$trend <- ATA.Inv.Transform(X=ATA.last$trend, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
	crit_a <- ifelse(is.season==TRUE, ifelse(output[6]==0,"A","M"), org.seas.Type)
	if(crit_a=="A"){
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata + SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata + OS_SIValue, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
	}else {
		ATA.fitted <- ATA.Inv.Transform(X = fit.ata * SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
		ATA.forecast <- ATA.Inv.Transform(X = forecast.ata * OS_SIValue, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))				
	}
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
	my_list$initial.level <- initialLevel
	my_list$initial.trend <- initialTrend
	my_list$level.fixed <- level.Fix
	my_list$trend.fixed <- trend.Fix
	my_list$transform.method <- transform.Method
	my_list$lambda <- Lambda
	my_list$accuracy.type <- accuracy.Type
	my_list$accuracy <- accuracy.ata
	my_list$is.season <- is.season
	my_list$seasonal.model <- ifelse(is.season==TRUE, switch(output[5]+1, "none", "decomp", "stl", "stlplus", "stR", "tbats", "x13", "x11"), "none")
	my_list$seasonal.type <- ifelse(is.season==TRUE, ifelse(output[6]==0,"A","M"), org.seas.Type)
	my_list$seasonal.period <- seasonal.Frequency
	my_list$seasonal.index <- ATA.Inv.Transform(X=SeasonalIndex, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
	my_list$seasonal <- ATA.Inv.Transform(X=SeasonalActual, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
	my_list$seasonal.adjusted <- ATA.Inv.Transform(X=AdjInput, tMethod=transform.Method, tLambda=Lambda, tbiasadj=boxcox_attr_set$bcBiasAdj, tfvar=ifelse(boxcox_attr_set$bcBiasAdj==FALSE, NULL, var(ATA.last$residuals)))
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