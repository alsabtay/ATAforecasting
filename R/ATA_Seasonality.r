#' @title Seasonality Tests for ATA Method 
#' @description This function is upgraded version of \code{\link{ndiffs}} and \code{\link{nsdiffs}} written by Hyndman et al. \code{forecast} package. Please review manual and vignette documents of latest \code{forecast} package.
#'
#' \strong{According to \code{forecast} package,}
#'
#' \code{ndiffs} and \code{nsdiffs} functions to estimate the number of differences required to make a given time series stationary. 
#' 
#' \code{ndiffs} estimates the number of first differences necessary.
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
#' \code{nsdiffs} estimates the number of seasonal differences necessary. \code{nsdiffs} uses seasonal unit root tests to determine the number of
#' seasonal differences required for time series to be made stationary.
#' 
#' Several different seasonal tests are available:
#' * If \code{suroot.test="seas"} (default), a measure of seasonal strength is used, where differencing is
#' selected if the seasonal strength (Wang, Smith & Hyndman, 2006) exceeds 0.64 
#' (based on minimizing MASE when forecasting using auto.arima on M3 and M4 data).
#' * If \code{suroot.test="ch"}, the Canova-Hansen (1995) test is used 
#' (with null hypothesis of deterministic seasonality) 
#' * If \code{suroot.test="hegy"}, the Hylleberg, Engle, Granger & Yoo (1990) test is used.
#' * If \code{suroot.test="ocsb"}, the Osborn-Chui-Smith-Birchenhall
#' (1988) test is used (with null hypothesis that a seasonal unit root exists).
#' * If \code{suroot.test="periodogram"}, This function is written based on M4 Competition Seasonality Test.
#' 
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link{forecast}}, \code{\link{ucra}}, \code{\link{tseries}}, \code{\link{uroot}}, \code{\link{stlplus}}, \code{stR}, 
#' \code{\link[stats]{stl}}, \code{\link[stats]{decompose}}, \code{\link{tbats}}, \code{\link{seasadj}}.
#' @references Dickey DA and Fuller WA (1979), "Distribution of the Estimators for
#' Autoregressive Time Series with a Unit Root", \emph{Journal of the American
#' Statistical Association} \bold{74}:427-431.
#'
#' Kwiatkowski D, Phillips PCB, Schmidt P and Shin Y (1992) "Testing the Null
#' Hypothesis of Stationarity against the Alternative of a Unit Root",
#' \emph{Journal of Econometrics} \bold{54}:159-178.
#'
#' Osborn DR, Chui APL, Smith J, and Birchenhall CR (1988) "Seasonality and the
#' order of integration for consumption", \emph{Oxford Bulletin of Economics
#' and Statistics} \bold{50}(4):361-377.
#'
#' Osborn, D.R. (1990) "A survey of seasonality in UK macroeconomic variables",
#' \emph{International Journal of Forecasting}, \bold{6}:327-336.
#'
#' Said E and Dickey DA (1984), "Testing for Unit Roots in Autoregressive
#' Moving Average Models of Unknown Order", \emph{Biometrika}
#' \bold{71}:599-607.

#' Wang, X, Smith, KA, Hyndman, RJ (2006) "Characteristic-based clustering
#' for time series data", \emph{Data Mining and Knowledge Discovery},
#' \bold{13}(3), 335-364.
#' 
#' Canova F and Hansen BE (1995) "Are Seasonal Patterns Constant
#' over Time? A Test for Seasonal Stability", \emph{Journal of Business and
#' Economic Statistics} \bold{13}(3):237-252.
#'
#' Hylleberg S, Engle R, Granger C and Yoo B (1990) "Seasonal integration
#' and cointegration.", \emph{Journal of Econometrics} \bold{44}(1), pp. 215-238.

#' @export ATA.Seasonality

ATA.Seasonality <- function(input, ppy, attr_set)
{
	if (ppy==1){
		test_seasonal <- FALSE
	}else {
		test <- attr_set$suroot.test
		if (test=="periodogram"){
			test_seasonal <- pgram.test(input, ppy, attr_set)
		}else {
			if (length(ppy)>1){
				if (attr_set$multi.period=="max"){
					ppy <- max(ppy)
				}else {
					ppy <- min(ppy)
				}
			}
			if (attr_set$suroot.uroot==TRUE){
				uroot.test <- attr_set$uroot.test
				uroot.type <- attr_set$uroot.type
				uroot.alpha <- attr_set$uroot.alpha
				uroot.pkg <- attr_set$uroot.pkg
				uroot.maxd <- attr_set$uroot.maxd
				#Used to determine whether a time series is stationary (trend)
				if (uroot.pkg=="ucra") {
					d <- forecast::ndiffs(input, alpha=uroot.alpha, test=uroot.test, type=uroot.type, max.d=uroot.maxd)
				}else {
					d <- ndiffs.tseries(input, alpha=uroot.alpha, test=uroot.test, max.d=uroot.maxd)
				}
				if (d > 0){
					input <- diff(input, differences=d, lag=1)
				}
			}
			suroot.alpha <- attr_set$suroot.alpha
			suroot.maxD <- attr_set$suroot.maxD
			suroot.m <- attr_set$suroot.m
			if (!is.null(attr_set$suroot.m)){
				D <- forecast::nsdiffs(input, alpha=suroot.alpha, m=suroot.m, test=test, max.D=suroot.maxD)
			}else {
				D <- forecast::nsdiffs(input, alpha=suroot.alpha, test=test, max.D=suroot.maxD)
			}
			test_seasonal <- ifelse(D==0, FALSE, TRUE) 
		}
	}
}

#' @export pgram.test

pgram.test <- function(input, ppy, attr_set)
{
	if (ppy==1){
		test_seasonal <- FALSE
	}else {
		pgram.tcrit <- attr_set$pgram.tcrit
		uroot.test <- attr_set$uroot.test
		uroot.type <- attr_set$uroot.type
		uroot.alpha <- attr_set$uroot.alpha
		uroot.pkg <- attr_set$uroot.pkg
		uroot.maxd <- attr_set$uroot.maxd
		if (length(ppy)>1){
			if (attr_set$multi.period=="max"){
				ppy <- max(ppy)
			}else {
				ppy <- min(ppy)
			}
		}
	  #Used to determine whether a time series is stationary (trend)
		if (uroot.pkg=="ucra") {
			d <- forecast::ndiffs(input, alpha=uroot.alpha, test=uroot.test, type=uroot.type, max.d=uroot.maxd)
		}else {
			d <- ndiffs.tseries(input, alpha=uroot.alpha, test=uroot.test, max.d=uroot.maxd)
		}
		if (d > 0){
			input <- diff(input, differences=d, lag=1)
		}
	  #Used to determine whether a time series is seasonal
		if (length(input)<3*ppy){
			test_seasonal <- FALSE
		}else {
			xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
			clim <- pgram.tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
			test_seasonal <- (abs(xacf[ppy]) > clim[ppy])
			if (is.na(test_seasonal)==TRUE){
				test_seasonal <- FALSE
			}
		}
	}
	return(test_seasonal)
}

#' @export ndiffs.tseries

ndiffs.tseries <- function(x, alpha=0.05, test=c("kpss","adf","pp"), max.d=2)
{
  #ndiffs function using tseries package
  test <- match.arg(test)
  x <- c(na.omit(c(x)))
  d <- 0

  if(is.constant(x))
    return(d)

  if(test=="kpss")
    suppressWarnings(dodiff <- tseries::kpss.test(x)$p.value < alpha)
  else if(test=="adf")
    suppressWarnings(dodiff <- tseries::adf.test(x)$p.value > alpha)
  else if(test=="pp")
    suppressWarnings(dodiff <- tseries::pp.test(x)$p.value > alpha)
  else
    stop("This shouldn't happen")
  if(is.na(dodiff))
  {
    return(d)
  }
  while(dodiff & d < max.d)
  {
    d <- d+1
    x <- diff(x)
    if(is.constant(x))
      return(d)
    if(test=="kpss")
      suppressWarnings(dodiff <- tseries::kpss.test(x)$p.value < alpha)
    else if(test=="adf")
      suppressWarnings(dodiff <- tseries::adf.test(x)$p.value > alpha)
    else if(test=="pp")
      suppressWarnings(dodiff <- tseries::pp.test(x)$p.value > alpha)
    else
      stop("This shouldn't happen")
    if(is.na(dodiff))
      return(d-1)
  }
  return(d)
}
