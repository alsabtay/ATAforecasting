#' Forecasting Method for The ATAforecasting
#'
#' @description \code{ATA.Forecast} is a generic function for forecasting of the ATA Method.
#'
#' @param object An \code{ATA} object is required for forecast.
#' @param h Number of periods for forecasting.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#' @param ci.level Confidence Interval levels for forecasting. Default value is 95.
#' @param negative.forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#' @param print.out Default is TRUE. If FALSE, forecast summary of ATA Method is not shown.
#'
#' @return An object of class "\code{ata}".
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{forecast}, \code{stlplus}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{tbats}, \code{seasadj}.
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
#'
#' @keywords Ata forecast accuracy ts msts
#'
#' @importFrom stats cycle end frequency start ts tsp tsp<- var
#' @importFrom Rdpack reprompt
#' @importFrom forecast msts
#'
#' @examples
#' demoATA <- window(fundingTR, start = tsp(fundingTR)[1], end = 2013)
#' ata.fit <- ATA(demoATA, parPHI = 1, seasonal.test = TRUE, seasonal.model = "decomp")
#' ATA.plot(ATA.Forecast(ata.fit, h=18))
#'
#'
#' @export
ATA.Forecast <- function(object, h=NULL, out.sample=NULL, ci.level=95, negative.forecast=TRUE, print.out = TRUE)
{
  y <- object
  if (class(object)!="ata"){
    return("The Input must be 'ATA' object. Please use ATA(x) function to produce 'ATA' object. ATA Forecast was terminated!")
  }
  m <- frequency(object$actual)
  if (is.null(h)){
    if (m==4){
      h <- 8
    }else if (m==5){
      h <- 10
    }else if (m==7){
      h <- 14
    }else if (m==12){
      h <- 24
    }else if (m==24){
      h <- 48
    }else {
      h <- 6
    }
    warning(paste("Input forecast horizon has been changed with ", h))
  }
  if(!is.null(out.sample)){
    if (length(out.sample)!=h){
      warning("Forecast horizon has been set to the length of out.sample.")
    }
  }
  tsp_y <- tsp(object$actual)
  fsample <- forecast::msts(rep(NA,h), start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  freqYh <- cycle(fsample)
  if (is.null(object$transform.method)){
    if (object$seasonal.model!="decomp" & object$seasonal.type=="M"){
      seasonal.type <- "A"
      lambda <- 0
	    shift <- 0
      transform.method <- "Box_Cox"
      bcBiasAdj <- FALSE
      out.transform <- ATA.Transform(object$seasonal.adjusted, tMethod = transform.method, tLambda = lambda, tShift = shift)       # lambda = 0 for multiplicative model
	    y$actual <- forecast::msts(out.transform$trfmX, start = start(object$actual), seasonal.periods = object$seasonal.period)
	    shift <- out.transform$tShift
	}else {
      seasonal.type <- y$seasonal.type
      lambda <- object$lambda
	    shift <- object$shift
      transform.method <- object$transform.method
      out.transform <- ATA.Transform(object$seasonal.adjusted, tMethod = transform.method, tLambda = lambda, tShift = shift)
	    y$actual <- forecast::msts(out.transform$trfmX, start = start(object$actual), seasonal.periods = object$seasonal.period)
	    shift <- out.transform$tShift
    }
  }else {
      seasonal.type <- object$seasonal.type
      lambda <- object$lambda
	    shift <- object$shift
      transform.method <- object$transform.method
      out.transform <- ATA.Transform(object$seasonal.adjusted, tMethod=transform.method, tLambda=lambda, tShift = shift)
	    y$actual <- forecast::msts(out.transform$trfmX, start = start (object$actual), seasonal.periods = object$seasonal.period)
	    shift <- out.transform$tShift
  }
  if (object$is.season==FALSE & seasonal.type=="A"){
    OS_SIValue <- rep(0,times=h)
  }else if (y$is.season==FALSE & seasonal.type=="M"){
    OS_SIValue <- rep(1,times=h)
  }else if (y$is.season==TRUE){
    OS_SIValue <- rep(NA,times=h)
    for (k in 1:h){
      OS_SIValue[k] <- object$seasonal.index[freqYh[k]]
    }
  }else{
  }
  OS_SIValue <- ATA.Transform(OS_SIValue, tMethod=transform.method, tLambda=lambda, tShift = shift)$trfmX
  y$level <- forecast::msts(ATA.Transform(y$level, tMethod=transform.method, tLambda=lambda, tShift = shift)$trfmX, start=tsp_y[1], seasonal.periods = y$seasonal.period)
  y$trend <- forecast::msts(ATA.Transform(y$trend, tMethod=transform.method, tLambda=lambda, tShift = shift)$trfmX, start=tsp_y[1], seasonal.periods = y$seasonal.period)
  y <- SubATA.Forecast(y, hh = h, initialLevel = object$initial.level)
  forecast.ata <- y$forecast
  if(y$seasonal.type=="A"){
    ATA.forecast <- ATA.BackTransform(X=forecast.ata + OS_SIValue, tMethod=transform.method, tLambda=lambda, tShift = shift, tbiasadj=y$bcBiasAdj, tfvar=ifelse(y$bcBiasAdj==FALSE, NULL, var(y$residuals)))
  }else {
    ATA.forecast <- ATA.BackTransform(X=forecast.ata * OS_SIValue, tMethod=transform.method, tLambda=lambda, tShift = shift, tbiasadj=y$bcBiasAdj, tfvar=ifelse(y$bcBiasAdj==FALSE, NULL, var(y$residuals)))
  }
  if (negative.forecast==TRUE){
    y$forecast <- forecast::msts(ATA.forecast, start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    y$forecast <- forecast::msts(ATA.forecast, start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  }
  y$h <- h
  accuracy.ata <- ATA.Accuracy(y, out.sample = out.sample)
  y$accuracy <- accuracy.ata
  y$out.sample <- ifelse(is.null(out.sample), fsample, out.sample)
  ci.output <- ATA.CI(object = y, ci.level = ci.level)
  y$ci.level <- ci.level
  if (negative.forecast==TRUE){
    y$forecast.lower <- ci.output$forecast.lower
    y$forecast.upper <- ci.output$forecast.upper
  }else {
    ci_low <- ci.output$forecast.lower
    ci_up <- ci.output$forecast.upper
    ci_low[ci_low<0] <- 0
    ci_up[ci_up<0] <- 0
    y$forecast.lower <- ci_low
    y$forecast.upper <- ci_up
  }
  attr(y, "class") <- "ata"
  print_out <- cbind(y$forecast.lower, y$forecast, y$forecast.upper)
  colnames(print_out) <- c("lower", "forecast", "upper")
  if (print.out==TRUE) {
    print(print_out)
  }
  gc()
  return(y)
}
