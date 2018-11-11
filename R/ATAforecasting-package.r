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
#' @docType package
#' @name ATAforecasting-package
#' @aliases ATAforecasting
#' @aliases ATAforecasting-package
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
#' \dontrun{
#' fit <- ATA(M3[[1899]]$x,M3[[1899]]$xx)
#' plot(ATA.Forecast(fit,h=36))
#'}
#'
#' @useDynLib ATAforecasting, .registration = TRUE
#' @exportPattern("^[[:alpha:]]+")
#' @importFrom(Rcpp, evalCpp)
#' @export(ATA)
#' @export(ATA.Accuracy)
#' @export(ATA.CI)
#' @export(ATA.Core)
#' @export(ATA.Decomposition)
#' @export(ATA.Forecast)
#' @export(ATA.Inv.Transform)
#' @export(ATA.Shift)
#' @export(ATA.Transform)
#' @export(ATA.Seasonality)
#' @export(ATA.SeasAttr)
#' @export(ATA.BoxCoxAttr)
#' @export(AutoATA.Accuracy)
#' @export(AutoATA.Accuracy.Holdout)
#' @export(AutoATA.Auto)
#' @export(AutoATA.Core)
#' @export(AutoATA.Core.Holdout)
#' @export(AutoATA.Damped)
#' @export(AutoATA.Forecast)
#' @export(find.freq.fourier)
#' @export(find.freq)
#' @export(find.freq.all)
#' @export(ndiffs.tseries)
#' @export(pgram.test)
#' @export(plot.ata)
#' @export(print.ata)
#' @import(Rcpp)
#' @import(RcppArmadillo)
#' @import(tseries)
#' @import(forecast)
#' @import(urca)
#' @import(uroot)
#' @import(seasonal)
#' @import(stR)
#' @import(stlplus)
#' @import(xts)
#' @import(timeSeries)
#' @import(TSA)
#' @import(Mcomp)
NULL

#' @export