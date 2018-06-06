#' Attributes set for unit root and seasonality tests
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
#' @param s.tcrit t-value for seasonality test
#' @param uroot.test Type of unit root test to use
#' @param uroot.type Specification of the deterministic component in the regression. Possible values are "level" and "trend".
#' @param uroot.alpha Level of the test, possible values range from 0.01 to 0.1
#' @param uroot.maxd Maximum number of non-seasonal differences allowed
#' @param x13.estimate.maxiter Maximum iteration for X13ARIMA/SEATS estimation
#' @param x13.estimate.tol Convergence tolerence for X13ARIMA/SEATS estimation
#' @param x11.estimate.maxiter Maximum iteration for X11 estimation
#' @param x11.estimate.tol Convergence tolerence for X11 estimation
#' @return An object of class \code{ataattrset}.
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}.
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
#'
#' @export
 type = c("level","trend")
ata.seasonal.attr <- function(s.tcrit=1.645, uroot.test="adf", uroot.type="level", uroot.alpha=0.05, uroot.maxd=2, x13.estimate.maxiter=1500, x13.estimate.tol=1.0e-5, x11.estimate.maxiter=1500, x11.estimate.tol=1.0e-5) 
{
	mylist <- list("s.tcrit"=s.tcrit,"uroot.test"=uroot.test,"uroot.type"=uroot.type, "uroot.alpha"=uroot.alpha, "uroot.maxd"=uroot.maxd, "x13.estimate.maxiter"=x13.estimate.maxiter, "x13.estimate.tol"=x13.estimate.tol, "x11.estimate.maxiter"=x11.estimate.maxiter, "x11.estimate.tol"=x11.estimate.tol)
	attr(mylist, "class") <- "ataattrset"
	return(mylist)
}
