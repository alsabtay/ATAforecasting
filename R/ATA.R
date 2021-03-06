#' ATAforecasting: Automatic Time Series Analysis and Forecasting using Ata Method with Box-Cox Power Transformations Family and Seasonal Decomposition Techniques
#'
#' @description Returns ATA(p,q,phi) applied to \code{X}.
#' The Ata method based on the modified simple exponential smoothing as described in Yapar, G. (2016) <doi:10.15672/HJMS.201614320580> , 
#' Yapar G., Capar, S., Selamlar, H. T., Yavuz, I. (2017) <doi:10.15672/HJMS.2017.493> and Yapar G., Selamlar, H. T., Capar, S., Yavuz, I. (2019) 
#' <doi:10.15672/hujms.461032> is a new univariate time series forecasting method which provides innovative solutions to issues faced during 
#' the initialization and optimization stages of existing methods.
#' Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting. 
#' It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
#' This methodology performed well on the M3 and M4-competition data.
#' Returns ATA(p,q,phi) applied to \code{X}.
#'
#' @docType package
#'
#' @name ATAforecasting-package
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' Maintainer: alisabritaylan@gmail.com
#'
#' @keywords package
NULL # Instead of "_PACKAGE" to remove inclusion of \alias{ATAforecasting}
# "_PACKAGE"


## Generic Ata Methods functions
## Part of ATAforecasting package


#' Automatic Time Series Analysis and Forecasting using Ata Method with Box-Cox Power Transformations Family and Seasonal Decomposition Techniques
#'
#' \code{ATA} is a generic function for Ata Method forecasting. 
#' The Ata method based on the modified simple exponential smoothing as described in Yapar, G. (2016) <doi:10.15672/HJMS.201614320580> , 
#' Yapar G., Capar, S., Selamlar, H. T., Yavuz, I. (2017) <doi:10.15672/HJMS.2017.493> and Yapar G., Selamlar, H. T., Capar, S., Yavuz, I. (2019) 
#' <doi:10.15672/hujms.461032> is a new univariate time series forecasting method which provides innovative solutions to issues faced during 
#' the initialization and optimization stages of existing methods.
#' Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting. 
#' It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
#' This methodology performed well on the M3 and M4-competition data.
#'
#' Returns ATA(p,q,phi) applied to \code{X}.
#'
#' @param X A numeric vector or time series of class \code{ts} or \code{msts} for in-sample.
#' @param Y A numeric vector or time series of class \code{ts} or \code{msts} for out-sample. If you do not have out-sample data, you can split in-sample data into training and test dataset with \code{partition.h} argument.
#' @param parP Value of Level parameter \code{p}. If NULL or "opt", it is estimated. \code{p} has all integer values from 1 to \code{length(X)}.
#' @param parQ Value of Trend parameter \code{q}. If NULL or "opt", it is estimated. \code{q} has all integer values from 0 to \code{p}.
#' @param parPHI Value of Damping Trend parameter \code{phi}. If NULL or "opt", it is estimated. phi has all values from 0 to 1.
#' @param model.type An one-character string identifying method using the framework terminology. The letter "A" for additive model, the letter "M" for multiplicative model.
#' If NULL, both letters will be tried and the best model (according to the accuracy measure \code{accuracy.type}) returned.
#' @param seasonal.test Testing for stationary and seasonality. If TRUE, the method firstly uses \code{test="adf"}, Augmented Dickey-Fuller, unit-root test then the test returns the least number of differences required to pass the test at level \code{alpha}.
#' After the unit-root test, seasonal test applies on the stationary \code{X}.
#' @param seasonal.model A string identifying method for seasonal decomposition. If NULL, "decomp" method is default. c("none", "decomp", "stl", "stlplus", "tbats", "stR") phrases of methods denote
#' \itemize{
#'		 \item{none}	: seasonal decomposition is not required.
#'		 \item{decomp} 	: classical seasonal decomposition. If \code{decomp}, the \code{stats} package will be used.
#'		 \item{stl}		: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stl}, the \code{stats} and \code{forecast} packages will be used. Multiple seasonal periods are allowed.
#'		 \item{stlplus}	: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stlplus}, the \code{stlplus} package will be used.
#'		 \item{tbats}   : exponential smoothing state space model with Box-Cox transformation, ARMA errors, trend and seasonal components.
#' 					  	  as described in De Livera, Hyndman & Snyder (2011). Parallel processing is used by default to speed up the computations. If \code{tbats}, the \code{forecast} package will be used. Multiple seasonal periods are allowed.
#'		 \item{stR}    	: seasonal-trend decomposition procedure based on regression developed by Dokumentov and Hyndman (2015). If \code{stR}, the \code{stR} package will be used. Multiple seasonal periods are allowed.
#'		 \item{x13}    	: seasonal-trend decomposition procedure based on X13ARIMA/SEATS. If \code{x13}, the \code{seasonal} package will be used.
#'		 \item{x11}    	: seasonal-trend decomposition procedure based on X11. If \code{x11}, the \code{seasonal} package will be used.
#' }
#' @param seasonal.period Value(s) of seasonal periodicity. If NULL, \code{frequency} of X is default  If \code{seasonal.period} is not integer, \code{X} must be \code{msts} time series object. c(s1,s2,s3,...) for multiple period. If \code{X} has multiple periodicity, "tbats" or "stR" seasonal model have to be selected.
#' @param seasonal.type	An one-character string identifying method for the seasonal component framework. The letter "A" for additive model, the letter "M" for multiplicative model.
#' If NULL, both letters will be tried and the best model (according to the accuracy measure \code{accuracy.type}) returned.
#' If seasonal decomposition methods except \code{decomp} with "M", Box-Cox transformation with \code{lambda}=0 is selected.
#' @param seasonal.test.attr Attributes set for unit root, seasonality tests, X13ARIMA/SEATS and X11. If NULL, corrgram.tcrit=1.28, uroot.test="adf", suroot.test="correlogram", suroot.uroot=TRUE, uroot.type="trend", uroot.alpha=0.05, suroot.alpha=0.05, uroot.maxd=2, suroot.maxD=1, suroot.m=frequency(X), uroot.pkg="urca", multi.period="min", x13.estimate.maxiter=1500, x13.estimate.tol=1.0e-5, x11.estimate.maxiter=1500, x11.estimate.tol=1.0e-5. If you want to change, please use \code{ATA.SeasAttr} function and its output.
#' For example, you can use \code{seasonal.test.attr = ATA.SeasAttr(corrgram.tcrit=1.65)} equation in \code{ATA} function.
#' @param find.period Find seasonal period(s) automatically. If NULL, 0 is default. When \code{find.period},
#' \itemize{
#'		 \item{0} : none
#'		 \item{1} : single period with find.freq
#'		 \item{2} : single period with \code{forecast::findfrequency}
#'		 \item{3} : multiple period with find.freq & stR
#'		 \item{4} : multiple period with find.freq & tbats
#'		 \item{5} : multiple period with find.freq & stl
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
#'		 \item{OWA}		: overall weighted average of MASE and sMAPE.
#'		 \item{MdAE}	: median absolute error.
#'		 \item{MdSE}	: median square error.
#'		 \item{RMdSE}	: root median squared error.
#'		 \item{MdPE}	: median percentage error.
#'		 \item{MdAPE}	: median absolute percentage error.
#'		 \item{sMdAPE}	: symmetric median absolute percentage error.
#' }
#' @param level.fixed If TRUE, "pStarQ"  --> First, fits ATA(p,0) where p = p* is optimized for q=0. Then, fits ATA(p*,q) where q is optimized for p = p*.
#' @param trend.fixed If TRUE, "pBullet" --> Fits ATA(p,1) where p = p* is optimized for q = 1.
#' @param trend.search If TRUE, "qBullet" --> Fits ATA(p,q) where p = p* is optimized for q = q* (q > 0). Then, fits ATA(p*,q) where q is optimized for p = p*.
#' @param h The number of steps to forecast ahead.
#' When the parameter is NULL; if the frequency of \code{X} is 4 the parameter is set to 8; if the frequency of \code{X} is 12 the parameter is set to 18; the parameter is set to 6 for other cases.
#' @param partition.h If \code{Y} is NULL, this parameter divides \code{X} into two parts: training set (in-sample) and test set (out-sample). \code{partition.h} is number of periods for forecasting and size of test set.
#' If the value is between 0 and 1, percentage of length is active. If \code{holdout} is TRUE, this parameter will be same as \code{h} for defining holdout set.
#' @param holdout Default is FALSE. If TRUE, ATA Method uses the holdout forecasting for accuracy measure to select the best model. In holdout forecasting, the last few data points are removed from the data series.
#' The remaining historical data series is called in-sample data (training set), and the holdout data is called out-of-sample data (holdout set).
#' If TRUE, partition.h will used for holdout data.
#' @param holdout.adjustedP Default is TRUE. If TRUE, parP will be adjusted by length of training - validation sets and in-sample set when the holdout forecasting is active.
#' @param holdin Default is FALSE. If TRUE, ATA Method uses the hold-in forecasting for accuracy measure to select the best model. In hold-in forecasting, the last h-length data points are used for accuracy measure.
#' @param transform.order If "before", Box-Cox transformation family will be applied and then seasonal decomposition techniques will be applied. If "after", seasonal decomposition techniques will be applied and then Box-Cox transformation family will be applied.
#' @param transform.method Transformation method  --> "Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog". If the transformation process needs shift parameter, 
#' \code{ATA.Transform} will calculate required shift parameter automatically.
#' @param transform.attr Attributes set for Box-Cox transformation. If NULL, bcMethod = "loglik", bcLower = 0, bcUpper = 1, bcBiasAdj = FALSE. If you want to change, please use \code{ATA.BoxCoxAttr} function and its output.
#' @param lambda Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.
#' @param shift Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated.
#' @param initial.level If NULL, FALSE is default. If FALSE, ATA Method calculates the pth observation in \code{X} for level.
#' If TRUE, ATA Method calculates average of first p value in \code{X}for level.
#' @param initial.trend If NULL, FALSE is default. If FALSE, ATA Method calculates the qth observation in \code{X(T)-X(T-1)} for trend.
#' If TRUE, ATA Method calculates average of first q value in \code{X(T)-X(T-1)} for trend.
#' @param ci.level Confidence Interval levels for forecasting.
#' @param start.phi Lower boundary for searching \code{parPHI}.If NULL, 0 is default.
#' @param end.phi Upper boundary for searching \code{parPHI}. If NULL, 1 is is default.
#' @param size.phi Increment step for searching \code{parPHI}. If NULL, 0.05 is default.
#' @param negative.forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#' @param print.out Default is TRUE. If FALSE, summary of ATA Method is not shown.
#' @param plot.out Default is TRUE. If FALSE, graphics of ATA Method are not shown.
#'
#' @return Returns an object of class \code{ATA}. The generic accessor functions \code{ATA.Forecast} and \code{ATA.Accuracy} extract useful features of the value returned by \code{ATA} and associated functions.
#' \code{ATA} object is a list containing at least the following elements
#' \itemize{
#' 		 \item{actual}		: The original time series.
#' 		 \item{fitted}		: Fitted values (one-step forecasts). The mean is of the fitted values is calculated over the ensemble.
#' 		 \item{level}		  : Estimated level values.
#' 		 \item{trend}		  : Estimated trend values.
#' 		 \item{residuals}	: Original values minus fitted values.
#' 		 \item{coefp}		  : The weights attached to level observations.
#' 		 \item{coefq}		  : The weights attached to trend observations.
#' 		 \item{p}		      : Optimum level parameter.
#' 		 \item{q}		      : Optimum trend parameter.
#' 		 \item{phi}		    : Optimum damped trend parameter.
#' 		 \item{model.type}: Form of trend.
#' 		 \item{h}		      : The number of steps to forecast ahead.
#' 		 \item{forecast}	: Point forecasts as a time series.
#' 		 \item{out.sample}: Test values as a time series.
#' 		 \item{method}		: The name of the optimum forecasting method as a character string.
#' 		 \item{initial.level}     : Selected initial level values for the time series forecasting method.
#' 		 \item{initial.trend}     : Selected initial trend values for the time series forecasting method.
#' 		 \item{level.fixed}       : A choice of optional level-fixed trended methods.
#' 		 \item{trend.fixed}       : A choice of optional trend-fixed trended methods.
#' 		 \item{trend.search}      : A choice of optional trend and level optimized trended methods if q > 1.
#' 		 \item{transform.method}  : Box-Cox power transformation family method  --> Box_Cox, Sqrt, Reciprocal, Log, NegLog, Modulus, BickelDoksum, Manly, Dual, YeoJohnson, GPower, GLog.
#' 		 \item{transform.order}   : Define how to apply Box-Cox power transformation techniques, before or after seasonal decomposition.
#' 		 \item{lambda}  	: Box-Cox power transformation family parameter.
#' 		 \item{shift}		  : Box-Cox power transformation family shifting parameter.
#' 		 \item{accuracy.type}		  : Accuracy measure that is chosen for model selection.
#' 		 \item{accuracy} 	: In and out sample accuracy measures and its descriptives that are calculated for optimum model are given.
#' 		 \item{holdout}		: Holdout forecasting is TRUE or FALSE.
#' 		 \item{holdout.training} 	: Training set in holdout forecasting.
#' 		 \item{holdout.validation}: Validation set in holdout forecasting.
#' 		 \item{holdout.forecast}	: Holdout forecast.
#' 		 \item{holdout.accuracy}	: Accuracy measure chosen for model selection in holdout forecasting.
#' 		 \item{holdin}		: Hold-in forecasting is TRUE or FALSE.
#' 		 \item{is.season}	: Indicates whether it contains seasonal pattern.
#' 		 \item{seasonal.model}		: The name of the selected decomposition method.
#' 		 \item{seasonal.type}	  	: Form of seasonality.
#' 		 \item{seasonal.period}		: The number of seasonality periods.
#' 		 \item{seasonal.index}		: Weights of seasonality.
#' 		 \item{seasonal}	: Estimated seasonal values.
#' 		 \item{seasonal.adjusted}	: Deseasonalized time series values.
#' 		 \item{execution.time}		: The real and CPU time 'in seconds' spent by the system executing that task, including the time spent executing run-time or system services on its behalf.
#' 		 \item{calculation.time}	: How much real time 'in seconds' the currently running R process has already taken.
#' }
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{forecast}, \code{stlplus}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{tbats}, \code{seasadj}, \code{seasonal}.
#'
#' @references
#' 
#' #'\insertRef{yapar2017mses}{ATAforecasting}
#'
#' #'\insertRef{yapar2018mhes}{ATAforecasting}
#'
#' #'\insertRef{yapar2018mses}{ATAforecasting}
#'
#' #'\insertRef{yapar2019ata}{ATAforecasting}
#'
#' @keywords Ata forecast accuracy ts msts
#' 
#' @importFrom forecast findfrequency
#' @importFrom stats cycle frequency ts tsp tsp<-
#' @importFrom Rdpack reprompt
#'
#' @examples
#' demoATA <-  window(touristTR, start = 2008, end = 2014.917)
#' ata.fit <- ATA(demoATA, h=18, seasonal.test = TRUE, seasonal.model = "stl")
#' ATA.plot(ATA.Forecast(ata.fit,h=18, out.sample = tail(touristTR,18)))
#'
#' @export
ATA <- function(X, Y = NULL,
                parP = NULL,
                parQ = NULL,
                parPHI = NULL,
                model.type = NULL,
                seasonal.test = NULL,
                seasonal.model = NULL,
                seasonal.period = NULL,
                seasonal.type = NULL,
                seasonal.test.attr = NULL,
                find.period = NULL,
                accuracy.type = NULL,
                level.fixed = FALSE,
                trend.fixed = FALSE,
                trend.search = FALSE,
                h = NULL,
                partition.h = NULL,
                holdout = FALSE,
                holdout.adjustedP = TRUE,
                holdin = FALSE,
                transform.order = "before",
                transform.method = NULL,
                transform.attr = NULL,
                lambda = NULL,
                shift = 0,
                initial.level = NULL,
                initial.trend = NULL,
                ci.level = 95,
                start.phi = NULL,
                end.phi = NULL,
                size.phi = NULL,
                negative.forecast = TRUE,
                print.out = TRUE,
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
      seasonal.period <- find.freq(X)
	  seasonal.test <- TRUE
    }else if(find.period==2){
      seasonal.period <- forecast::findfrequency(X)
	  seasonal.test <- TRUE
    }else if (find.period==3){
      seasonal.period <- find.multi.freq(X)
      seasonal.model=="stR"
	  seasonal.test <- TRUE
    }else if (find.period==4){
      seasonal.period <- find.multi.freq(X)
      seasonal.model=="tbats"
	  seasonal.test <- TRUE
    }else if (find.period==5){
      seasonal.period <- find.multi.freq(X)
      seasonal.model=="stl"
	  seasonal.test <- TRUE
    }else {
      return("find.period must be integer and between 0 and 5. ATA Method was terminated!")
    }
    s.frequency <- seasonal.period
  }
  if (length(s.frequency)>1){
    if (seasonal.model != "tbats" & seasonal.model != "stR" & seasonal.model != "stl"){
      seasonal.model <- "stl"
	  seasonal.test <- TRUE
    }
  }else {
	if (s.frequency > 1 & is.null(seasonal.model)){
		seasonal.model <- "decomp"
		seasonal.test <- TRUE
	}
  }
  if (is.null(accuracy.type)){
    accuracy.type <- "sMAPE"
  }
  if (trend.search==TRUE){
    level.fixed <- FALSE
    trend.fixed <- FALSE
  }else if (level.fixed==TRUE & trend.fixed==TRUE){
    level.fixed <- FALSE
    trend.fixed <- TRUE
  }else {
  }
  if (is.null(initial.level)){
    initial.level = FALSE
  }
  if (is.null(initial.trend)){
    initial.trend = FALSE
  }
  if (is.null(seasonal.test.attr)) {
    seas_attr_set <- ATA.SeasAttr()
  }else {
    seas_attr_set <- seasonal.test.attr
  }
  if (!is.null(seasonal.type)){
    if (is.null(seasonal.test)){
      seasonal.test <- TRUE
    }
  }else {
	seasonal.test <- FALSE
  }
  if (is.null(transform.attr)) {
    boxcox_attr_set <- ATA.BoxCoxAttr()
  }else {
    boxcox_attr_set <- transform.attr
  }
  if (holdout == TRUE & holdin == TRUE){
    return("Only one parameter of the two parameters (holdout or holdin) must be selected. Please choose one one of them. ATA Method was terminated!")
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
  if ((accuracy.type != "MAE" & accuracy.type != "MSE" & accuracy.type != "RMSE" & accuracy.type != "MPE" & accuracy.type != "MAPE" & accuracy.type != "sMAPE" & accuracy.type != "MASE" & accuracy.type != "OWA" & accuracy.type != "MdAE" & accuracy.type != "MdSE" & accuracy.type != "MdPE" & accuracy.type != "MdAPE" & accuracy.type != "sMdAPE") | !is.character(accuracy.type) | length(accuracy.type) > 1){
    return("Accuracy Type value must be string and it must get one value: MAE or MSE or MPE or MAPE or sMAPE or MASE or MdAE or MdSE or MdPE or MdAPE or sMdAPE. ATA Method was terminated!")
  }
  if (!is.null(model.type)){
    if ((model.type != "A" & model.type != "M") | !is.character(model.type) | length(model.type) > 1){
      return("Model Type value must be string. A for additive or M for multiplicative or NULL for both of them. ATA Method was terminated!")
    }
  }
  if (!is.null(initial.level)){
    if (initial.level != FALSE & initial.level != TRUE) {
      return("Initial value for Level must be boolean and it must get one value: TRUE or FALSE. ATA Method was terminated!")
    }
  }
  if (!is.null(initial.trend)){
    if (initial.trend != FALSE & initial.trend != TRUE) {
      return("Initial value for Trend must be boolean and it must get one value: TRUE or FALSE. ATA Method was terminated!")
    }
  }
  if (!is.null(transform.order)){
    if ((transform.order != "before" & transform.order != "after") | !is.character(transform.order) | length(transform.order) > 1){
      return("Transformation Order value must be string. 'before' for Transformation --> Decompostion or 'after' for Decomposition --> Transformation. ATA Method was terminated!")
    }
  }
  if (!is.null(transform.method)){
    if ((transform.method != "Box_Cox" & transform.method != "Box_Cox_Shift" & transform.method != "Modulus" & transform.method != "BickelDoksum" & transform.method != "Dual"  & transform.method != "Manly"  & transform.method != "Sqrt" & transform.method != "SqrtShift" & 
         transform.method != "YeoJohnson" & transform.method != "GPower" & transform.method != "GLog" & transform.method != "Log" & transform.method != "Reciprocal" & transform.method !="ReciprocalShift" & transform.method != "NegLog" & 
		 transform.method != "LogShift") | !is.character(transform.method) | length(transform.method) > 1){
      return("Transform Method value must be string. Please select a valid Box-Cox transformation technique. ATA Method was terminated!")
    }
  }
  if (class(seas_attr_set)!="ataattrset"){
    return("Attributes set for unit root and seasonality tests are not suitable set. ATA Method was terminated!")
  }
  if (class(boxcox_attr_set)!="ataattrset"){
    return("Attributes set for Box-Cox transformation are not suitable set. ATA Method was terminated!")
  }

  WD <- getwd()
  start.time <- Sys.time()
  ptm <- proc.time()
  class_X <- class(X)
  X <- ts(X, frequency = min(s.frequency), start = tsp(X)[1])
  firstTspX <- tsp(X)
  X <- forecast::msts(X, start=firstTspX[1], seasonal.periods = s.frequency)
  tspX <- tsp(X)
  if (!is.null(Y[1])){
    OutSample <- Y
    h <- length(Y)
    if (holdout == TRUE){
      if (is.null(partition.h)){
        partition.h	<- h
      }
    }else {
      partition.h <- NULL
    }
  }else {
    if (!is.null(partition.h)){
      part_h <- ifelse(partition.h > 0 & partition.h < 1, floor(length(X) * partition.h), partition.h)
      OSLen <- length(X)- part_h
      ISLen <- length(X)
      OutSample <- X[(OSLen+1):ISLen]
      OutSample <- ts(OutSample, frequency = tspX[3], start = tspX[2] - ifelse(tspX[3]>1, (part_h - 1) * (1/tspX[3]), (part_h - 1) * 1))
      X <- X[1:OSLen]
      X <- ts(X, frequency = tspX[3], start = tspX[1])
      h <- length(OutSample)
    }else {
      if (is.null(h)){
        if (max(s.frequency)==4){
          h <- 8
        }else if (max(s.frequency)==12){
          h <- 18
        }else {
          h <- 6
        }
      }
      OutSample <- rep(NA,times=h)
      OutSample <- ts(OutSample, frequency = tspX[3], start = tspX[2] + ifelse(tspX[3]>1, 1/tspX[3], 1))
      if (holdout == TRUE){
        partition.h	<- h
      }
    }
  }
  orig.X <- ts(X, frequency = firstTspX[3], start = firstTspX[1])
  tsp_ny <- tsp(forecast::msts(orig.X, start=tsp(orig.X)[1], seasonal.periods = s.frequency))
  fsample <- ts(OutSample, frequency = max(s.frequency), start = tsp_ny[2] + ifelse(tsp_ny[3]>1, 1/tsp_ny[3], 1))
  freqYh <- cycle(fsample)
  if (transform.order == "before"){
    ChgX <- ATA.Transform(X,tMethod=transform.method, tLambda=lambda, tShift=shift, bcMethod = boxcox_attr_set$bcMethod, bcLower = boxcox_attr_set$bcLower, bcUpper = boxcox_attr_set$bcUpper)
    X <- forecast::msts(ChgX$trfmX, start=firstTspX[1], seasonal.periods = s.frequency)
    lambda <- ChgX$tLambda
    shift <- ChgX$tShift
    if (length(seasonal.type)==1 & length(seasonal.model)==1){
      my_list <- AutoATA.Single(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                                level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method,
                                lambda, shift, orig.X, OutSample, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, partition.h, holdout.adjustedP, holdin)
    }else {
      my_list <- AutoATA.Multiple(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                                  level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method,
                                  lambda, shift, orig.X, OutSample, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, partition.h, holdout.adjustedP, holdin)
    }
  }else {
    if (!is.null(transform.method)){
      model.type <- "M"
      warning("model.type parameter has been set as 'M' because of a transformation techniques from Box-Cox power transformation family selected.")
    }
    if (length(seasonal.type)==1 & length(seasonal.model)==1){
      my_list <- AutoATA.SingleO(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                                 level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method,
                                 lambda, shift, orig.X, OutSample, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, partition.h, holdout.adjustedP, holdin)
    }else {
      my_list <- AutoATA.MultipleO(X, parP, parQ, model.type, seasonal.test, seasonal.model, seasonal.type, s.frequency, h, accuracy.type,
                                   level.fixed, trend.fixed, trend.search, start.phi, end.phi, size.phi, initial.level, initial.trend, transform.method,
                                   lambda, shift, orig.X, OutSample, seas_attr_set, freqYh, ci.level, negative.forecast, boxcox_attr_set, holdout, partition.h, holdout.adjustedP, holdin)
    }
  }
  executionTime <- proc.time() - ptm
  end.time <- Sys.time()
  my_list$transform.order <- transform.order
  my_list$execution.time <- executionTime
  my_list$calculation.time <- round(as.double(difftime(end.time, start.time,units="sec")),4)
  attr(my_list, "class") <- "ATA"
  if (plot.out==TRUE) {
    ATA.plot(my_list)
  }
  gc()
  ATA.print(my_list)
  return(my_list)
}



#' Specialized Screen Print Function of The ATAforecasting
#'
#' @param object an object of \code{ATA}
#' @param ... other inputs
#'
#' @return a summary for the results of the ATAforecasting
#'
#' @export
ATA.print <- function(object,...)
{
	opscipen <- options("scipen" = 100, "digits"=7)
    on.exit(options(opscipen))
	x <- object
    cat(x$method,"\n\n")
    if (x$level.fixed==TRUE){
      cat("   level.fixed: TRUE","\n\n")
    }
    if (x$trend.fixed==TRUE){
      cat("   trend.fixed: TRUE","\n\n")
    }
    if(!is.null(x$transform.method)){
      cat(paste("   '",x$transform.method, "' transformation method was selected.","\n\n", sep=""))
    }
    if(!is.null(x$lambda)){
      cat("   Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
    }
    cat(paste("   model.type:",x$model.type, "\n\n"))
    if (x$is.season==FALSE){
      cat("   seasonal.model: no seasonality","\n\n")
    }else {
      cat(paste("   seasonal.model:",x$seasonal.model, "\n\n"))
    }
    if (x$is.season==TRUE){
      cat(paste("   seasonal.type:",x$seasonal.type, "\n\n"))
    }
    cat(paste("   forecast horizon:",x$h, "\n\n"))
    cat(paste("   accuracy.type:",x$accuracy.type, "\n\n"))

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$inSample$MAE, x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MPE$inSample$MPE, x$accuracy$MAPE$inSample$MAPE, x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$MASE$inSample$MASE, x$accuracy$OWA$inSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$inSample$MdAE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$RMdSE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$sMAPE$inSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")


    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MPE$outSample$MPE, x$accuracy$MAPE$outSample$MAPE, x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$MASE$outSample$MASE, x$accuracy$OWA$outSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MdAE, x$accuracy$MSE$outSample$MdSE, x$accuracy$MSE$outSample$RMdSE, x$accuracy$MPE$outSample$MdPE, x$accuracy$MAPE$outSample$MdAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    stats <- c(x$execution.time[1], x$execution.time[2], x$execution.time[3])
    names(stats) <- c("user","system","elapsed")
    cat("\n")
    print(stats)
    cat("\n")
    cat(paste("calculation.time:",x$calculation.time, "\n\n"))
    cat("\n")

    cat("Forecasts:","\n")
    print(x$forecast)
    cat("\n\n")
}



#' Specialized Plot Function of The ATAforecasting
#'
#' @param object an object of \code{ATA}
#' @param fcol line color
#' @param flty line type
#' @param flwd line width
#' @param ... other inputs
#'
#' @return a graphic output for the components of the ATAforecasting
#'
#' @importFrom stats cycle frequency ts tsp tsp<- 
#' @importFrom graphics axis legend layout lines mtext par plot polygon 
#'
#' @export
ATA.plot <- function(object, fcol=4, flty = 2, flwd = 2, ...)
{
  x <- object
  oldpar <- par(no.readonly = TRUE)# save default, for resetting...
  on.exit(par(oldpar))
  caption <- paste(ifelse(x$model.type=="A"," Additive "," Multiplicative "), x$method, sep="")
  xx <- x$actual
  hpred <- length(x$forecast)
  freq <- frequency(xx)
  xxx <- ts(c(x$actual, rep(NA,hpred)), end=tsp(xx)[2] + hpred/freq, frequency=freq)
  xxy <- ts(c(x$fitted, rep(NA,hpred)), end=tsp(xx)[2] + hpred/freq, frequency=freq)
  min_y <- min(x$actual, x$fitted, x$out.sample, x$forecast, x$forecast.lower, na.rm=TRUE)
  max_y <- max(x$actual, x$fitted, x$out.sample, x$forecast, x$forecast.upper, na.rm=TRUE)
  range_y <- abs(max_y - min_y)
  min_last <- floor(min_y - range_y * 0.20)
  max_last <- ceiling(max_y + range_y * 0.20)
  range_last <- abs(max_last - min_last)
  rnd_par <- ifelse(range_last<10, 1, 0)
  dataset <- cbind(xxx,xxy)
  colnames(dataset, do.NULL = FALSE)
  colnames(dataset) <- c("actual","fitted")
  legend_names <- c("actual","fitted","out-sample","forecast")
  tmp <- seq(from = tsp(x$forecast)[1], by = 1/freq, length = hpred)
  if (x$is.season==FALSE){
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(dataset,plot.type="s", ylim=c(min_last, max_last), col=1:ncol(dataset), xlab=NULL, ylab="fitted", yaxt="n")
    axis(side=2,at=seq(min_last, max_last,round(range_last/10, rnd_par)), labels=seq(min_last, max_last,round(range_last/10, rnd_par)), las=1, lwd=1)
    polygon(x=c(tmp, rev(tmp)), y=c(x$forecast.lower, rev(x$forecast.upper)), col="lightgray", border=NA)
    lines(x$forecast, lty = flty, lwd = flwd, col = fcol)
    lines(x$out.sample, lty = 1, lwd = flwd, col = fcol+2)
    legend("topleft", legend_names, col=c(1,2,fcol+2,fcol), lty=1, cex=.80, box.lty=0, text.font=2, ncol=2,  bg="transparent")
    mtext(caption, side = 3, line = -1.5, outer = TRUE)
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$trend, ylab="trend")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$level, ylab="level")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$residuals, ylab="residuals")
  }else { 
    layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow=TRUE))
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(dataset,plot.type="s", ylim=c(min_last, max_last), col=1:ncol(dataset), xlab=NULL, ylab="fitted", yaxt="n")
    axis(side=2,at=seq(min_last, max_last,round(range_last/10, rnd_par)), labels=seq(min_last, max_last,round(range_last/10, rnd_par)), las=1, lwd=1)
    polygon(x=c(tmp, rev(tmp)), y=c(x$forecast.lower, rev(x$forecast.upper)), col="lightgray", border=NA)
    lines(x$forecast, lty = flty, lwd = flwd, col = fcol)
    lines(x$out.sample, lty = 1, lwd = flwd, col = fcol+2)
    legend("topleft", legend_names, col=c(1,2,fcol+2,fcol), lty=1, cex=.80, box.lty=0, text.font=2, ncol=2, bg="transparent")
    mtext(paste(caption,"with ",ifelse(x$seasonal.type=="A","Additive","Multiplicative"), " Decomposition by '",ifelse(x$seasonal.model=="decomp","classical",x$seasonal.model),"' Method"), side = 3, line = -1.5, outer = TRUE)
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$seasonal.adjusted,ylab="deseasonalized")
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$level, ylab="level")
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$trend,ylab="trend")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$seasonal,ylab="seasonality")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$residuals, ylab="residuals")
  }
  #par(oldpar)
}

find.precision <- function(x)
{
	results <- nchar(sub(".", "", x, fixed=TRUE))
	return(results)
}

find.scale <- function(x)
{
	results <- nchar(sub("\\d+\\.?(.*)$", "\\1", x))
	return(results)
}
