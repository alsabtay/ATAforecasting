#' @title Automatic Seasonal Decomposition for ATA Method
#' @description Returns seasonally adjusted data constructed by removing the seasonal component. The methodology is fully automatic.
#' @aliases print.ata plot.ata 
#' @param input It must be \code{ts} or \code{msts} or \code{numeric} object. if it is \code{numeric} object, \code{findPeriod} must be 1 or 2 or 3 or 4. if it is \code{msts} object, \code{findPeriod} must be 3 or 4.
#' @param s.model A string identifying method for seasonal decomposition. If NULL, "decomp" method is default. c("none", "decomp", "stl", "stlplus", "tbats", "stR") phrases of methods denote.
#' \itemize{
#'		 \item{none}	: seasonal decomposition is not required.
#'		 \item{decomp} 	: classical seasonal decomposition. If \code{decomp}, the \code{stats} package will be used. 
#'		 \item{stl}		: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stl}, the \code{stats} package will be used. 
#'		 \item{stlplus}	: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stlplus}, the \code{stlplus} package will be used. 
#'		 \item{tbats}   : exponential smoothing state space model with box-cox transformation, ARMA errors, trend and seasonal components 
#' 					  	  as described in De Livera, Hyndman & Snyder (2011). Parallel processing is used by default to speed up the computations. If \code{tbats}, the \code{forecast} package will be used. 
#'		 \item{stR}    	: seasonal-trend decomposition procedure based on regression developed by Dokumentov and Hyndman (2015). If \code{stR}, the \code{stR} package will be used. 
#'		 \item{x13}    	: seasonal-trend decomposition procedure based on X13ARIMA/SEATS. If \code{x13}, the \code{seasonal} package will be used. 
#'		 \item{x11}    	: seasonal-trend decomposition procedure based on X11. If \code{x11}, the \code{seasonal} package will be used. 
#' }
#' @param s.type A one-character string identifying method for the seasonal component framework. If NULL, "M" is default. The letter "A" for additive model, the letter "M" for multiplicative model.
#' @param s.frequency Value(s) of seasonal periodicity. If \code{s.frequency} is not integer, \code{X} must be \code{msts} time series object. c(s1,s2,s3,...) for multiple period. If \code{X} has multiple periodicity, "tbats" or "stR" seasonal model have to be selected.
#' For example: period of the input data which have one seasonal pattern --> 12 for monthly / 4 for quarterly / 7 for daily / 5 for business days. periods of the input data which have complex/multiple seasonal patterns --> c(7,354.37,365.25).
#' @return Univariate time series.
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}, \code{\link{stlplus}}, \code{stR}, \code{\link{seasonal}}.
#' @keywords ata seasonal decomposition forecast accuracy ts msts
#' @export ATA.Decomposition

ATA.Decomposition <- function(input, s.model, s.type, s.frequency, seas_attr_set)
{
tsp_input <- tsp(input)	
if (s.model == "none" | max(s.frequency)==1){
	if (s.type=="A"){	
		adjX <- input
		SeasActual <- rep(0,times=length(input))
		SeasActual <- ts(SeasActual, f=tsp_input[3], s=tsp_input[1])
		s.frequency <- frequency(input)
		SeasIndex <- rep(0,times=s.frequency) 
	}else {
		adjX <- input
		SeasActual <- rep(1,times=length(input))
		SeasActual <- ts(SeasActual, f=tsp_input[3], s=tsp_input[1])
		s.frequency <- frequency(input)
		SeasIndex <- rep(1,times=s.frequency) 
	}	
}else { 
	if (class(input)!="ts" & class(input)!="msts"){
		return("The data set must be time series object (ts or msts) ATA Method was terminated!")
	}
	input <- msts(input, start=tsp_input[1], seasonal.periods = s.frequency)
	tsp_input <- tsp(input)	
	if (s.model=="decomp"){									  	# Do classical decomposition
		if (s.type=="A"){
			desX <- decompose(input, type = c("additive"))
			adjX <- forecast::seasadj(desX)
			SeasActual <- desX$seasonal
			SeasIndex <- rep(NA,times=s.frequency)
			for (s in 1:s.frequency){
				SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
			}
		}else {
			desX <- decompose(input, type = c("multiplicative"))
			adjX <- forecast::seasadj(desX)
			SeasActual <- desX$seasonal
			SeasIndex <- rep(NA,times=s.frequency)
			for (s in 1:s.frequency){
				SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
			}
		}
	}else if (s.model=="stl"){									# Do STL decomposition
		stldesX <- stl(input, s.window = "per", robust=TRUE)
		adjX <- forecast::seasadj(stldesX)  
		SeasActual <- forecast::seasonal(stldesX)
		SeasIndex <- rep(NA,times=s.frequency)
		for (s in 1:s.frequency){
			SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
		}
	}else if (s.model=="stlplus"){								# Do STLPlus decomposition
		stlplusdesX <- stlplus(input, s.window = "per", robust=TRUE)
		adjX <- input - stlplusdesX$data$seasonal
		SeasActual <- stlplus::seasonal(stlplusdesX)
		SeasActual <- msts(SeasActual, start=tsp_input[1], seasonal.periods = s.frequency)
		SeasIndex <- rep(NA,times=s.frequency)
		for (s in 1:s.frequency){
			SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
		}
	}else if (s.model=="stR"){									# Do stR decomposition
		stRdesX <- AutoSTR(input, robust=TRUE)
		stRcomp <- stR::components(stRdesX)
		nameCol <- colnames(stRcomp)
		nameCol <- grep('Seasonal', nameCol, value=TRUE)
		if (length(nameCol)==0){
			if (s.type=="A"){	
				adjX <- input 
				SeasActual <- msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(0,times=max(s.frequency))
			}else {
				adjX <- input
				SeasActual <- msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(1,times=max(s.frequency))
			}	
		}else {
			adjX <- stR::seasadj(stRdesX)
			if (length(s.frequency)==1){
				SeasActual <- stRcomp[,nameCol]
				SeasIndex <- rep(NA,times=s.frequency)
				for (s in 1:s.frequency){
					SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
				}
			}else {
				SeasActual <- rowSums(stRcomp[,nameCol],na.rm=TRUE)
				SeasActual <- msts(SeasActual, start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(NA,times=max(s.frequency))
				for (s in 1:max(s.frequency)){
					SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
				}
			}
		}
	}else if (s.model=="tbats"){								# Do tbats decomposition
		tbatsdesX <- tbats(input)
		tbatscomp <- tbats.components(tbatsdesX)
		nameCol <- colnames(tbatscomp)
		nameCol <- grep('season', nameCol, value=TRUE)
		if (length(nameCol)==0){
			if (s.type=="A"){	
				adjX <- input 
				SeasActual <- msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(0,times=max(s.frequency))
			}else {
				adjX <- input
				SeasActual <- msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(1,times=max(s.frequency))
			}	
		}else {
			adjX <- forecast::seasadj(tbatsdesX)
			if (length(s.frequency)==1){
				SeasActual <- tbatscomp[,nameCol]
				SeasIndex <- rep(NA,times=s.frequency)
				for (s in 1:s.frequency){
					SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
				}
			}else {
				SeasActual <- rowSums(tbatscomp[,nameCol],na.rm=TRUE)
				SeasActual <- msts(SeasActual, start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(NA,times=max(s.frequency))
				for (s in 1:max(s.frequency)){
					SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
				}
			}
		}
	}else if (s.model=="x13"){									# Do X13ARIMA/SEATS decomposition
		x13desX <- seas(input, estimate.maxiter=seas_attr_set$x13.estimate.maxiter, estimate.tol=seas_attr_set$x13.estimate.tol)
		SeasActual <- seasonal::series(x13desX,"seats.adjustfac")
		if (is.null(SeasActual)) {
			if (s.type=="A"){	
				adjX <- input 
				SeasActual <- msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(0,times=max(s.frequency))
			}else {
				adjX <- input
				SeasActual <- msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(1,times=max(s.frequency))
			}	
		}else {
			adjX <- seasonal::series(x13desX,"seats.seasonaladj")
			SeasIndex <- rep(NA,times=s.frequency)
			for (s in 1:s.frequency){
					SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
			}
		}
	}else if (s.model=="x11"){									# Do X13ARIMA/SEATS X11 decomposition
		x11desX <- seas(input, x11 = "", estimate.maxiter=seas_attr_set$x11.estimate.maxiter, estimate.tol=seas_attr_set$x11.estimate.tol)
		SeasActual <- seasonal::series(x11desX,"x11.adjustfac")
		if (is.null(SeasActual)) {
			if (s.type=="A"){	
				adjX <- input 
				SeasActual <- msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(0,times=max(s.frequency))
			}else {
				adjX <- input
				SeasActual <- msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
				SeasIndex <- rep(1,times=max(s.frequency))
			}	
		}else {
			adjX <- seasonal::series(x11desX,"x11.seasadj")
			SeasIndex <- rep(NA,times=s.frequency)
			for (s in 1:s.frequency){
					SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
			}
		}
	}else {
	}
}
my_list <- list("AdjustedX" = adjX, "SeasIndex" = SeasIndex, "SeasActual" = SeasActual)
return(my_list)
gc()
}	
