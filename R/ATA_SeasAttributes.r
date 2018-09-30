#' @title Attributes set for unit root and seasonality tests
#' @description This function is a class of seasonality tests using \code{pgram.test}, \code{\link{ndiffs}} and \code{\link{nsdiffs}} functions. 
#' 
#' If \code{uroot.test="kpss"}, the KPSS test is used with the null hypothesis that
#' \code{x} has a stationary root against a unit-root alternative. Then the
#' test returns the least number of differences required to pass the test at
#' the level \code{uroot.alpha}. If \code{uroot.test="adf"}, the Augmented Dickey-Fuller
#' test is used and if \code{uroot.test="pp"} the Phillips-Perron test is used. In
#' both of these cases, the null hypothesis is that \code{x} has a unit root
#' against a stationary root alternative. Then the test returns the least
#' number of differences required to fail the test at the level \code{alpha}.
#'
#' \code{nsdiffs} uses seasonal unit root tests to determine the number of
#' seasonal differences required for time series to be made stationary.
#' 
#' Several different tests are available:
#' * If \code{suroot.test="seas"} (default), a measure of seasonal strength is used, where differencing is
#' selected if the seasonal strength (Wang, Smith & Hyndman, 2006) exceeds 0.64 
#' (based on minimizing MASE when forecasting using auto.arima on M3 and M4 data).
#' * If \code{suroot.test="ch"}, the Canova-Hansen (1995) test is used 
#' (with null hypothesis of deterministic seasonality) 
#' * If \code{suroot.test="hegy"}, the Hylleberg, Engle, Granger & Yoo (1990) test is used.
#' * If \code{suroot.test="ocsb"}, the Osborn-Chui-Smith-Birchenhall
#' (1988) test is used (with null hypothesis that a seasonal unit root exists).
#' 
#' @param pgram.tcrit t-value for periodogram seasonality test.
#' @param uroot.test Type of unit root test before all type seasonality test. Possible values are "adf", "pp" and "kpss". 
#' @param suroot.test Type of seasonal unit root test to use. Possible values are "periodogram", "seas", "hegy", "ch" and "ocsb". 
#' @param suroot.uroot If TRUE, unit root test for stationary before seasonal unit root test is allowed.
#' @param uroot.type Specification of the deterministic component in the regression for unit root test. Possible values are "level" and "trend".
#' @param uroot.alpha Significant level of the unit root test, possible values range from 0.01 to 0.1
#' @param suroot.alpha Significant level of the seasonal unit root test, possible values range from 0.01 to 0.1
#' @param uroot.maxd Maximum number of non-seasonal differences allowed.
#' @param suroot.maxD Maximum number of seasonal differences allowed.
#' @param suroot.m Deprecated. Length of seasonal period: frequency of data for nsdiff.
#' @param uroot.pkg Using \code{ucra} or \code{tseries} packages for unit root test. The default value is "ucra".
#' @param multi.period Selection type of multi seasonal period. \code{min} or \code{max} function for selection
#' @param x13.estimate.maxiter Maximum iteration for X13ARIMA/SEATS estimation
#' @param x13.estimate.tol Convergence tolerence for X13ARIMA/SEATS estimation
#' @param x11.estimate.maxiter Maximum iteration for X11 estimation
#' @param x11.estimate.tol Convergence tolerence for X11 estimation
#' @return An object of class \code{ataattrset}.
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}.
#'

#' @export ATA.SeasAttributes

ATA.SeasAttributes <- function(pgram.tcrit=1.28, uroot.test="adf", suroot.test="periodogram", suroot.uroot=TRUE, uroot.type="trend", uroot.alpha=0.05, suroot.alpha=0.05, uroot.maxd=2, suroot.maxD=1, suroot.m=NULL, uroot.pkg="ucra", multi.period="min", x13.estimate.maxiter=1500, x13.estimate.tol=1.0e-5, x11.estimate.maxiter=1500, x11.estimate.tol=1.0e-5) 
{
	if ((uroot.test != "adf" & uroot.test != "pp" & uroot.test != "kpss") | !is.character(uroot.test)){	
		return("Selection method of unit root test must be string. adf, pp or kpss test for searching unit root. ATA Method was terminated!")
	}
	if ((suroot.test != "periodogram" & suroot.test != "seas" & suroot.test != "hegy" & suroot.test != "ch" & suroot.test != "ocsb") | !is.character(suroot.test)){	
		return("Selection method of seasonal unit root test must be string. periodogram, seas, hegy, ch or ocsb test for searching seasonal unit root. ATA Method was terminated!")
	}
	if ((uroot.type != "level" & uroot.type != "trend") | !is.character(uroot.type)){	
		return("Selection type of deterministic component in the regression for unit root test must be string. level or trend for searching unit root. ATA Method was terminated!")
	}	
	if ((multi.period != "min" & multi.period != "max" ) | !is.character(multi.period)){	
		return("Selection type of multi seasonal period must be string. min or max function for selection. ATA Method was terminated!")
	}
	if(uroot.alpha < 0.01){
		warning("Specified alpha value is less than the minimum, setting uroot.alpha=0.01")
		uroot.alpha <- 0.01
	}else if(uroot.alpha > 0.1){
		warning("Specified alpha value is larger than the maximum, setting uroot.alpha=0.1")
		uroot.alpha <- 0.1
	}
	if(suroot.alpha < 0.01){
		warning("Specified alpha value is less than the minimum, setting suroot.alpha=0.01")
		suroot.alpha <- 0.01
	}else if(suroot.alpha > 0.1){
		warning("Specified alpha value is larger than the maximum, setting suroot.alpha=0.1")
		suroot.alpha <- 0.1
	}
	mylist <- list("pgram.tcrit"=pgram.tcrit, "uroot.test"=uroot.test, "suroot.test"=suroot.test, "suroot.uroot"=suroot.uroot, "uroot.type"=uroot.type, "uroot.alpha"=uroot.alpha, "suroot.alpha"=suroot.alpha, "uroot.maxd"=uroot.maxd, "suroot.maxD"=suroot.maxD, "suroot.m"=suroot.m, "uroot.pkg"=uroot.pkg, "multi.period"=multi.period, "x13.estimate.maxiter"=x13.estimate.maxiter, "x13.estimate.tol"=x13.estimate.tol, "x11.estimate.maxiter"=x11.estimate.maxiter, "x11.estimate.tol"=x11.estimate.tol)
	attr(mylist, "class") <- "ataattrset"
	return(mylist)	
}
