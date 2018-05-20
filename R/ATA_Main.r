#' @title Forecasting Time Series by ATA Method with Box-Cox transformation
#' @description Returns ATA(p,q,phi) applied to \code{X}.
#' Based on the modified simple exponential smoothing as described in Yapar, G. (2016). 
#' ATA method is a new univariate time series forecasting method which provides innovative 
#' solutions to issues faced during the initialization and optimization stages of existing methods. 
#' ATA's forecasting performance is superior to existing methods both in terms of easy implementation 
#' and accurate forecasting. It can be applied to non-seasonal or deseasonalized time series, 
#' where the deseasonalization can be performed via any preferred decomposition method.
#' This methodology performed extremely well on the M3 and M4-competition data.
#'
#' @aliases print.ata plot.ata 
#' @param X a numeric vector or time series of class \code{ts} or \code{msts} for in-sample.
#' @param Y a numeric vector or time series of class \code{ts} or \code{msts} for out-sample. If you do not have out-sample data, you can split in-sample data into training and test dataset with \code{partition.h} argument. 
#' @param parP Value of Level parameter \code{p}. If NULL or "opt", it is estimated. \code{p} has all integer values from 1 to \code{length(X)}.
#' @param parQ Value of Trend parameter \code{q}. If NULL or "opt", it is estimated. \code{q} has all integer values from 0 to \code{p}.
#' @param parPHI Value of Damping Trend parameter \code{phi}. If NULL or "opt", it is estimated. phi has all values from 0 to 1.
#' @param model.type A one-character string identifying method using the framework terminology. The letter "A" for additive model, the letter "M" for multiplicative model. 
#' If NULL, both letters will be tried and the best model (according to the accuracy measure \code{accuracy.type}) returned.
#' @param seasonal.test Testing for stationary and seasonality. If TRUE, the method firstly uses \code{test="adf"}, Augmented Dickey-Fuller, unit-root test then the test returns the least number of differences required to pass the test at level \code{alpha}.
#' After the unit-root test, seasonal test applies on the stationary \code{X}.
#' @param seasonal.model A string identifying method for seasonal decomposition. If NULL, "decomp" method is default. c("none", "decomp", "stl", "stlplus", "tbats", "stR") phrases of methods denote
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
#' @param seasonal.period Value(s) of seasonal periodicity. If NULL, \code{frequency} of X is default  If \code{seasonal.period} is not integer, \code{X} must be \code{msts} time series object. c(s1,s2,s3,...) for multiple period. If \code{X} has multiple periodicity, "tbats" or "stR" seasonal model have to be selected.
#' @param seasonal.type	A one-character string identifying method for the seasonal component framework. If NULL, "M" is default. The letter "A" for additive model, the letter "M" for multiplicative model.
#' If other seasonal decomposition method except \code{decomp} with "M", Box-Cox transformation with \code{lambda}=0 is selected.
#' @param seasonal.test.attr Attributes set for unit root, seasonality tests, X13ARIMA/SEATS and X11. If NULL, s.tcrit=1.645, uroot.test="adf", uroot.alpha=0.05, uroot.maxd=2, x13.estimate.maxiter=1500, x13.estimate.tol=1.0e-5, x11.estimate.maxiter=1500, x11.estimate.tol=1.0e-5. If you want to change, please use \code{ata.seasonal.attr} function and its output.
#' For example, you can use \code{seasonal.test.attr = ata.seasonal.attr(s.tcrit=1.96)} equation in \code{ATA} function. 
#' @param find.period Find seasonal period(s) automatically. If NULL, 0 is default. When \code{find.period},
#' \itemize{
#'		 \item{0} : none
#'		 \item{1} : single period with find.freq
#'		 \item{2} : single period with \code{forecast::findfrequency}
#'		 \item{3} : multiple period with find.freq & stR
#'		 \item{4} : multiple period with find.freq & tbats
#' }
#' @param accuracy.type	Accuracy measure for selection of the best model. IF NULL, \code{sMAPE} is default. 
#' \itemize{
#'		 \item{MAE}		: mean absolute error.
#'		 \item{MSE}		: mean square error.
#'		 \item{RMSE}	: root mean squared error.
#'		 \item{MPE}		: mean percentage error.
#'		 \item{MAPE}	: mean absolute percentage error.
#'		 \item{sMAPE}	: symmetric mean absolute percentage error.
#'		 \item{MASE}	: mean absolute scaled error.
#'		 \item{MdAE}	: median absolute error.
#'		 \item{MdSE}	: median square error.
#'		 \item{MdPE}	: median percentage error.
#'		 \item{MdAPE}	: median absolute percentage error.
#'		 \item{sMdAPE}	: symmetric median absolute percentage error.
#' }
#' @param level.fixed "pStarQ"  --> First, fits ATA(p,0) where p = p* is optimized for q=0. Then, fits ATA(p*,q) where q is optimized for p = p*.
#' @param trend.fixed "pBullet" --> Fits ATA(p,1) where p = p* is optimized for q = 1.
#' @param partition.h If \code{Y} is NULL, this parameter divides \code{X} into two parts: training set (in-sample) and test set (out-sample). \code{partition.h} is number of periods for forecasting and size of test set. 
#' When the parameter is NULL; if the frequency of \code{X} is 4 the parameter is set to 8; if the frequency of \code{X} is 12 the parameter is set to 18; the parameter is set to 6 for other cases.		
#' @param transform.method Transformation method  --> BoxCox, sqrt, inverse, log, log10. 
#' When \code{BoxCox} or \code{log} or \code{log10}  is specified,
#' \code{model.type} and \code{seasonal.type} is set to "A".
#' @param lambda Box-Cox transformation parameter. If NULL, data transformed before model is estimated. 
#' When \code{lambda} is specified, \code{model.type} and \code{seasonal.type} is set to "A".
#' @param initial.value If NULL, FALSE is default. If FALSE, ATA Method calculates the pth observation in \code{X} for level and qth observation in \code{X(T)-X(T-1)} for trend. 
#' If TRUE, ATA Method calculates average of first p value in \code{X}for level and average of first q value in \code{X(T)-X(T-1)} for trend.
#' @param ci.level Confidence Interval levels for forecasting.
#' @param start.phi Lower boundary for searching \code{parPHI}.If NULL, 0 is default.
#' @param end.phi Upper boundary for searching \code{parPHI}. If NULL, 1 is is default.
#' @param size.phi Increment step for searching \code{parPHI}. If NULL, 0.05 is default.
#' @param print.out Default is TRUE. If FALSE, summary of ATA Method is not shown. 
#' @param plot.out Default is TRUE. If FALSE, graphics of ATA Method are not shown. 
#' @param ... Other undocumented arguments.
#' @return Returns an object of class "\code{ata}", containing.

#' @param actual The original time series.
#' @param fitted Fitted values (one-step forecasts). The mean is of the fitted values is calculated over the ensemble.
#' @param level Estimated level values.
#' @param trend Estimated trend values.
#' @param residuals Original values minus fitted values.
#' @param coefp The weights attached to level observations.
#' @param coefq	The weights attached to trend observations.
#' @param p	Optimum level parameter.
#' @param q	Optimum trend parameter.
#' @param phi Optimum damped trend parameter.
#' @param model.type Form of trend.
#' @param h The number of steps to forecast ahead.
#' @param forecast Point forecasts as a time series. 
#' @param out.sample Test values as a time series.
#' @param method The name of the optimum forecasting method as a character string.
#' @param initial.value Selected initial values for the time series forecasting method.
#' @param level.fixed A choice of optional level fixed trended methods.
#' @param trend.fixed A choice of optional trend fixed trended methods.
#' @param transform.method Transformation method  --> BoxCox, sqrt, inverse, log, log10.
#' @param lambda Box-Cox transformation parameter.
#' @param accuracy.type Accuracy measure that is chosen for model selection.
#' @param accuracy In and out sample accuracy measures and its descriptives that are calculated for optimum model are given.
#' @param is.season Indicates whether it contains seasonal pattern.
#' @param seasonal.model The name of the selected decomposition method.
#' @param seasonal.type Form of seasonality.
#' @param seasonal.period The number of seasonality periods (which defaults to \code{frequency(X)}).
#' @param seasonal.index Weights of seasonality.
#' @param seasonal Estimated seasonal values.
#' @param seasonal.adjusted Deseasonalized time series values.
#' @param execution.time The real and CPU time (in seconds) spent by the system executing that task, including the time spent executing run-time or system services on its behalf.
#' @param calculation.time How much real time (in seconds) the currently running R process has already taken.
#'
#' The generic accessor functions \code{ATA.Forecast} and \code{ata.accuracy}
#' extract useful features of the value returned by \code{ata} and associated
#' functions.
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}, \code{\link{seasonal}}.
#' @references Yapar, G., (2016)
#' "Modified simple exponential smoothing"
#' \emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi:10.15672/HJMS.201614320580
#'
#' Yapar, G., Capar, S., Selamlar, H. T., Yavuz, I., (2016) 
#' "Modified holt's linear trend method"
#' \emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi: 10.15672/HJMS.2017.493 
#'
#' @keywords ata forecast accuracy ts msts
#' @examples
#' fit <- ATA(M3[[1899]]$x,M3[[1899]]$xx)
#' plot(ATA.Forecast(fit,h=36))
#'
#' @export ATA

ATA <- function(X, Y=NULL, 
					parP=NULL, 
					parQ=NULL, 
					parPHI=NULL, 
					model.type=NULL,
					seasonal.test=NULL,
					seasonal.model=NULL,
					seasonal.period=NULL,
					seasonal.type=NULL,
					seasonal.test.attr=NULL,
					find.period=NULL,
					accuracy.type=NULL, 
					level.fixed=FALSE,
					trend.fixed=FALSE,
					partition.h=NULL,
					transform.method=NULL,
					lambda=NULL,
					initial.value=NULL,
					ci.level=95,
					start.phi=NULL,
					end.phi=NULL,
					size.phi=NULL,
					# print.out = TRUE,
					plot.out = TRUE)
{
	if (class(X)[1]!="ts" & class(X)[1]!="msts"){	
		return("Class of X must be ts/msts object with single/multiple period seasonality. ATA Method was terminated!")
	}
	if (is.null(parQ)){
		parQ <- "opt"
	}
	if (is.null(parP)){
		parP <- "opt"
	}
	if (is.null(parPHI)){
		parPHI <- "opt"
		if (is.null(size.phi)){
			size.phi <- 0.05
		}
		if (is.null(start.phi)){
			start.phi <- 0
		}
		if (is.null(end.phi)){
			end.phi <- 1
		}
	}else {
		if (parPHI == "opt"){
			if (is.null(size.phi)){
				size.phi <- 0.05
			}
			if (is.null(start.phi)){
				start.phi <- 0
			}
			if (is.null(end.phi)){
				end.phi <- 1
			}
		}else {
			start.phi <- parPHI
			end.phi <- parPHI
			size.phi <- 1
		}
	}
	if (!is.null(seasonal.period)){
		find.period <- 0
		s.frequency <- seasonal.period
	}else{
		if (is.null(find.period)){
			find.period <- 0
		}
		s.frequency <- frequency(X)
		seasonal.period <- frequency(X)
	}
	if (find.period!=0) {
		if(find.period==1){			
			seasonal.period <- find.freq(input) 
		}else if(find.period==2){			
			seasonal.period <- findfrequency(input)
		}else if (find.period==3){
			seasonal.period <- find.freq.all(input)
			seasonal.model=="stR" 
		}else if (find.period==4){
			seasonal.period <- find.freq.all(input)
			seasonal.model=="tbats" 
		}else {
			return("find.period must be integer and between 0 and 4. ATA Method was terminated!")
		}
		s.frequency <- seasonal.period
	}
	if (is.null(accuracy.type)){ 
		accuracy.type <- "sMAPE"
	}
	if (level.fixed==TRUE & trend.fixed==TRUE){
		level.fixed <- FALSE 
		trend.fixed <- TRUE
	}
	if (is.null(initial.value)){
		initial.value = FALSE
	}
	if (is.null(seasonal.test.attr)) {
		seas_attr_set <- ata.seasonal.attr()
	}else {
		seas_attr_set <- seasonal.test.attr
	}
	if (!is.null(seasonal.type)){
		if (is.null(seasonal.test)){
			seasonal.test <- TRUE
		}
	}
	Qlen <- length(parQ)
	Plen <- length(parP)
	if (class(parP) =="character" & parP!="opt"){
			return("p value must be integer and between 1 and length of input. ATA Method was terminated!")
	}else if ((class(parP)=="numeric" |class(parP)=="integer") & Plen>1 & (max(parP)>length(X))){
			return("p value must be integer and between 1 and length of input. ATA Method was terminated!")
	}else{
	}
	if (class(parQ)=="character" & parQ!="opt"){
			return("p value must be integer and between 0 and p. ATA Method was terminated!")
	}else if ((class(parP)=="numeric" |class(parP)=="integer") & Qlen>1 & (max(parQ)>=max(parP))){
			return("q value must be integer and between 0 and p. ATA Method was terminated!")
	}else{
	}
	if (class(parPHI)=="character" & parPHI!="opt"){
		return("phi value must be numeric and between 0 and 1. ATA Method was terminated!")
	}else if ((class(parPHI)=="numeric" |class(parPHI)=="integer") & parPHI<0 & parPHI>1 & length(parPHI)>1){	
		return("phi value must be numeric and between 0 and 1. ATA Method was terminated!")
	}	
	if (!is.null(seasonal.type)){
		if ((seasonal.type != "A" & seasonal.type != "M") | !is.character(seasonal.type) | length(seasonal.type) > 1){	
			return("Seasonal Type value must be string. A for additive or M for multiplicative. ATA Method was terminated!")
		}
	}
	if (!is.null(seasonal.model)){
		if (length(seasonal.model) == 1){
			if ((seasonal.model != "none" & seasonal.model != "decomp" & seasonal.model != "stl" & seasonal.model != "stlplus" & seasonal.model != "tbats" & seasonal.model != "stR" & seasonal.model != "x13" & seasonal.model != "x11") | !is.character(seasonal.model)){	
				return("Seasonal Decomposition Model value must be string: decomp, stl, stlplus, tbats, stR. ATA Method was terminated!")
			}
		}else {
			if(any(seasonal.model %in% c("decomp","stl", "stlplus", "stR", "tbats", "x13", "x11"))){

			}else {
				return("Seasonal Decomposition Model value must be string: decomp, stl, stlplus, tbats, stR. ATA Method was terminated!")
			}
		}
	}
	if ((accuracy.type != "MAE" & accuracy.type != "MSE" & accuracy.type != "MPE" & accuracy.type != "MAPE" & accuracy.type != "sMAPE" & accuracy.type != "MASE" & accuracy.type != "MdAE" & accuracy.type != "MdSE" & accuracy.type != "MdPE" & accuracy.type != "MdAPE" & accuracy.type != "sMdAPE") | !is.character(accuracy.type) | length(accuracy.type) > 1){		
		return("Accuracy Type value must be string and it must get one value: MAE or MSE or MPE or MAPE or sMAPE or MASE or MdAE or MdSE or MdPE or MdAPE or sMdAPE. ATA Method was terminated!")
	}
	if (!is.null(model.type)){
		if ((model.type != "A" & model.type != "M") | !is.character(model.type) | length(model.type) > 1){	
			return("Model Type value must be string. A for additive or M for multiplicative or NULL for both of them. ATA Method was terminated!")
		}
	}
	if (!is.null(initial.value)){
		if (initial.value != FALSE & initial.value != TRUE) {	
			return("Initial value must be boolean and it must get one value: TRUE or FALSE. ATA Method was terminated!")
		}
	}
	if (!is.null(transform.method)){
		if ((transform.method != "sqrt" & transform.method != "BoxCox" & transform.method != "inverse" & transform.method != "log" & transform.method != "log10") | !is.character(transform.method) | length(transform.method) > 1){	
			return("Model Type value must be string. A for additive or M for multiplicative or NULL for both of them. ATA Method was terminated!")
		}
	}
	if (class(seas_attr_set)!="ataattrset"){
		return("Attributes Set for unit root and seasonality tests are not suitable set. ATA Method was terminated!")
	}
	
	WD <- getwd()
	start.time <- Sys.time()
	ptm <- proc.time()
	class_X <- class(X)
	tspX <- tsp(X)
	if (initial.value==FALSE){
		initial.level <- FALSE
		initial.trend <- FALSE
	}else {
		initial.level <- TRUE
		initial.trend <- TRUE
	}
	if (!is.null(Y[1])){
		OutSample <- Y
		h <- length(Y)
		partition.h <- NULL
	}else {
		if (!is.null(partition.h)){
			OSLen <- length(X)- partition.h
			ISLen <- length(X)
			OutSample <- X[(OSLen+1):ISLen]
			OutSample <- ts(OutSample, f = tspX[3], s = tspX[2] - ifelse(tspX[3]>1, (partition.h - 1) * (1/tspX[3]), (partition.h - 1) * 1))
			X <- X[1:OSLen]
			X <- ts(X, f = tspX[3], s = tspX[1])
			h <- length(OutSample)
		}else {
			if (s.frequency==4){
				h <- 8
			}else if (s.frequency==12){
				h <- 18
			}else {
				h <- 6
			}
			OutSample <- rep(NA,times=h)
			OutSample <- ts(OutSample, f = tspX[3], s = tspX[2] + ifelse(tspX[3]>1, 1/tspX[3], 1))
		}
	}
	orig.X <- X
	freqYh <- cycle(OutSample)
	if (!is.null(transform.method)){
		if (transform.method=="log"){ 
			model.type <- "A"
			seasonal.type <- "A"
			lambda <- 0
			transform.method <- "BoxCox"
		}
		if (transform.method=="log10"){
			model.type <- "A"
			seasonal.type <- "A"
		}
	}
	if (length(seasonal.type)==1 & length(seasonal.model)==1){
		orig.seastype <- seasonal.type
		if (seasonal.model=="none"){
			is.season <- FALSE
		}else {
			if (seasonal.test==FALSE){
				if (s.frequency==1){
					is.season <- FALSE
				}else {
					is.season <- TRUE
				}
			}else {
				is.season <- SeasonalityTest(X, s.frequency, seas_attr_set)
			}
		}
		if (is.season==TRUE){
			if (seasonal.model=="x13" | seasonal.model=="x11"){
				lambda <- NULL
				transform.method <- NULL
			}else if (seasonal.model!="decomp" & seasonal.type=="M"){
				X <- ATA.Transform(X,tMethod="BoxCox",tLambda=0)$trfmX  # lambda = 0 for multiplicative model
				seasonal.type <- "A"
				model.type <- "A"
				lambda <- 0
				transform.method <- "BoxCox"
			}else {
				tX <- ATA.Transform(X,tMethod=transform.method,tLambda=lambda)
				X <- tX$trfmX
				lambda <- tX$tLambda
			}
		}else {
			tX <- ATA.Transform(X,tMethod=transform.method,tLambda=lambda)
			X <- tX$trfmX
			lambda <- tX$tLambda
			seasonal.model <- "none"
			seasonal.type <- "A"
		}
		ata.seasonal.component <- ATA.Decomposition(X, s.model=seasonal.model, s.type=seasonal.type, s.frequency=s.frequency, seas_attr_set=seas_attr_set)
		AdjInSample <- ata.seasonal.component$AdjustedX
		SeasonalIndex <- ata.seasonal.component$SeasIndex
		SeasonalActual <- ata.seasonal.component$SeasActual
		if (seasonal.model=="x13" | seasonal.model=="x11"){
			if (abs(min(SeasonalIndex))<=1){
				seasonal.type <- "M"
			}else {
				seasonal.type <- "A"
			}
			orig.seastype <- seasonal.type
		}
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
		if (is.numeric(parP) & is.numeric(parQ) & is.numeric(parPHI) & !is.null(model.type)){
			ata.output <- ATA.Core(AdjInSample, pk = parP, qk = parQ, phik = parPHI, mdlType = model.type, initialLevel = initial.level, initialTrend = initial.trend)
			ata.output$h <- h
			ata.output <- AutoATA.Forecast(ata.output, hh=h, initialLevel = initial.level)
			ata.output$actual <- orig.X
			ata.output$accuracy.type <- accuracy.type
			fit.ata <- ata.output$fitted
			forecast.ata <- ata.output$forecast
			ata.output$level <- ATA.Inv.Transform(X=ata.output$level, tMethod=transform.method, tLambda=lambda)
			ata.output$trend <- ATA.Inv.Transform(X=ata.output$trend, tMethod=transform.method, tLambda=lambda)
			if(seasonal.type=="A"){
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata + SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata + OS_SIValue, tMethod=transform.method, tLambda=lambda)
			}else {
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata * SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata * OS_SIValue, tMethod=transform.method, tLambda=lambda)				
			}
			SeasonalActual <- ATA.Inv.Transform(X=SeasonalActual, tMethod=transform.method, tLambda=lambda)
			ata.output$fitted <- ATA.fitted
			ata.output$forecast <- ATA.forecast
			accuracy.ata <- ATA.Accuracy(ata.output, OutSample)
		}else if (is.numeric(parP) & is.numeric(parQ) & is.numeric(parPHI) & is.null(model.type)){
			mdl.type <- c("A","M")
			optAccryStart <- 9999999999999.9
			for (m in 1:2){	
				ATA.opt <- AutoATA.Core(AdjInSample, pk = parP, qk = parQ, phik = parPHI, mdlType = mdl.type[m], initialLevel = initial.level, initialTrend = initial.trend )
				optAccryEnd <- AutoATA.Accuracy(ATA.opt, accryType = accuracy.type)
				if (optAccryEnd <= optAccryStart){
					model.type <- mdl.type[m]
					optAccryStart <- optAccryEnd
				}
			}
			ata.output <- ATA.Core(AdjInSample, pk = parP, qk = parQ, phik = parPHI, mdlType = model.type, initialLevel = initial.level, initialTrend = initial.trend)
			ata.output$h <- h
			ata.output <- AutoATA.Forecast(ata.output, hh=h, initialLevel = initial.level)
			ata.output$actual<- orig.X
			fit.ata <- ata.output$fitted
			forecast.ata <- ata.output$forecast
			ata.output$level <- ATA.Inv.Transform(X=ata.output$level, tMethod=transform.method, tLambda=lambda)
			ata.output$trend<- ATA.Inv.Transform(X=ata.output$trend, tMethod=transform.method, tLambda=lambda)
			if(seasonal.type=="A"){
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata + SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata + OS_SIValue, tMethod=transform.method, tLambda=lambda)
			}else {
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata * SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata * OS_SIValue, tMethod=transform.method, tLambda=lambda)				
			}
			SeasonalActual <- ATA.Inv.Transform(X=SeasonalActual, tMethod=transform.method, tLambda=lambda)
			ata.output$fitted <- ATA.fitted
			ata.output$forecast <- ATA.forecast
			accuracy.ata <- ATA.Accuracy(ata.output, OutSample)	
		}else {		
			ata.output <- AutoATA.Damped(AdjInSample, pb = parP, qb = parQ, model.Type = model.type, accuracy.Type = accuracy.type, level.fix = level.fixed, trend.fix = trend.fixed, phiStart = start.phi, phiEnd = end.phi, phiSize = size.phi, initialLevel = initial.level, initialTrend = initial.trend)	
			ata.output$h <- h
			ata.output <- AutoATA.Forecast(ata.output, hh=h, initialLevel = initial.level)
			ata.output$actual <- orig.X
			fit.ata <- ata.output$fitted
			forecast.ata <- ata.output$forecast
			ata.output$level <- ATA.Inv.Transform(X=ata.output$level, tMethod=transform.method, tLambda=lambda)
			ata.output$trend <- ATA.Inv.Transform(X=ata.output$trend, tMethod=transform.method, tLambda=lambda)
			if(seasonal.type=="A"){
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata + SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata + OS_SIValue, tMethod=transform.method, tLambda=lambda)
			}else {
				ATA.fitted <- ATA.Inv.Transform(X=fit.ata * SeasonalActual, tMethod=transform.method, tLambda=lambda)
				ATA.forecast <- ATA.Inv.Transform(X=forecast.ata * OS_SIValue, tMethod=transform.method, tLambda=lambda)				
			}
			SeasonalActual <- ATA.Inv.Transform(X=SeasonalActual, tMethod=transform.method, tLambda=lambda)
			ata.output$fitted <- ATA.fitted
			ata.output$forecast <- ATA.forecast
			accuracy.ata <- ATA.Accuracy(ata.output, OutSample)
		}
		my_list <- ata.output
		my_list$out.sample <- OutSample
		if (level.fixed==TRUE){
			method <- paste("ATA(",my_list$p, ",", my_list$q,",", my_list$phi, ")", sep="")
		}else if (trend.fixed==TRUE){
				method <- paste("ATA(", my_list$p, ",1," ,my_list$phi, ")", sep="")
		}else {
			method <- paste("ATA(", my_list$p, "," ,my_list$q, ",", my_list$phi, ")", sep="")
		}
		my_list$method <- method
		my_list$initial.value <- initial.value
		my_list$level.fixed <- level.fixed
		my_list$trend.fixed <- trend.fixed
		my_list$transform.method <- transform.method
		my_list$lambda <- lambda
		my_list$accuracy.type <- accuracy.type
		my_list$accuracy <- accuracy.ata
		my_list$is.season <- is.season
		my_list$seasonal.model <- seasonal.model
		my_list$seasonal.type <- orig.seastype
		my_list$seasonal.period <- s.frequency
		my_list$seasonal.index <- SeasonalIndex
		my_list$seasonal <- SeasonalActual
		my_list$seasonal.adjusted <- ATA.Inv.Transform(X=AdjInSample, tMethod=transform.method, tLambda=lambda)
		ci.output <- ATA.CI(my_list, ci.level)
		my_list$ci.level <- ci.level
		my_list$forecast.lower <- ci.output$forecast.lower
		my_list$forecast.upper <- ci.output$forecast.upper
	}else {
		my_list <- AutoATA.Auto(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type, 
									level.fixed, trend.fixed, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method, 
									lambda, orig.X, OutSample, seas_attr_set, freqYh, ci.level)
	}
	executionTime <- proc.time() - ptm
	end.time <- Sys.time()
	my_list$execution.time <- executionTime
	my_list$calculation.time <- round(as.double(difftime(end.time, start.time,units="sec")),4)
	attr(my_list, "class") <- "ata"
	if (plot.out==TRUE) {
		plot.ata(my_list)
	}
	# if (print.out==TRUE) {
		# print.ata(my_list)
	# }
	gc()
	return(my_list)
}