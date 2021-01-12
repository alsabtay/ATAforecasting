#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @importFrom graphics axis legend layout lines mtext par plot polygon 
#' @importFrom forecast BoxCox BoxCox.lambda findfrequency InvBoxCox is.constant mstl msts ndiffs nsdiffs seasadj seasonal tbats tbats.components
#' @importFrom seasonal seas series udg
#' @importFrom stats acf as.ts cycle decompose frequency median na.omit qnorm qt sd ts tsp tsp<- spec.ar stl var
#' @importFrom stlplus stlplus
#' @importFrom stR AutoSTR components
#' @importFrom timeSeries colKurtosis colSkewness 
#' @importFrom TSA periodogram
#' @importFrom tseries adf.test kpss.test pp.test
#' @importFrom utils head tail 
#' @importFrom xts period.apply
#'
#' @useDynLib ATAforecasting, .registration = TRUE
NULL

#' @exportPattern("^[[:alpha:]]+")
