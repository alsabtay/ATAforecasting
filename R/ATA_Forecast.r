#' Forecasting Time Series Using the ATA Method
#'
#' @description \code{ATA.Forecast} is a generic function for forecasting of the ATA Method.
#' The function invokes particular \emph{methods} which depend on the class of the first argument.
#'
#' @param y An \code{ata} object is required for forecast.
#' @param h Number of periods for forecasting.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#' @param ci.level Confidence Interval levels for forecasting. Default value is 95.
#' @param negative.forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#'
#' @return An object of class "\code{ata}".
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}.
#'
#' @references Yapar, G., (2016)
#' "Modified simple exponential smoothing"
#' \emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi:10.15672/HJMS.201614320580
#'
#' Yapar, G., Capar, S., Selamlar, H. T., Yavuz, I., (2016)
#' "Modified holt's linear trend method"
#' \emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi: 10.15672/HJMS.2017.493
#'
#' @keywords ata forecast accuracy ts msts
#'
#' @export
#'
#' @examples
#' ata.fit <- ATA(insample, outsample)
#' fc <- ATA.Forecast(ata.fit, h=18)
ATA.Forecast <- function(y, h=NULL, out.sample=NULL, ci.level=95, negative.forecast=TRUE)
{
  if (class(y)!="ata"){
    return("The Input must be 'ata' object. Please use ATA(x) function to produce 'ata' object. ATA Forecast was terminated!")
  }
  if (is.null(h)){
    if (frequency(y$actual)==4){
      h <- 8
    }else if (frequency(y$actual)==12){
      h <- 24
    }else {
      h <- 7
    }
  }
  if(!is.null(out.sample)){
    if (length(out.sample)==h){
      return("The length of out.sample must be equal h. ATA Forecast was terminated!")
    }
  }
  ata.output <- y
  tsp_y <- tsp(y$actual)
  fsample <- ts(rep(NA,h), f = tsp_y[3], s = tsp_y[2] + ifelse(tsp_y[3]>1, 1/tsp_y[3], 0))
  freqYh <- cycle(fsample)
  if (is.null(y$transform.method)){
    if (y$seasonal.model!="decomp" & y$seasonal.type=="M"){
      seasonal.type <- "A"
      lambda <- 0
      transform.method <- "BoxCox"
      bcBiasAdj <- FALSE
      ata.output$actual <- ATA.Transform(y$seasonal.adjusted, tMethod=transform.method, tLambda=lambda)$trfmX       # lambda = 0 for multiplicative model
    }else {
      seasonal.type <- y$seasonal.type
      lambda <- y$lambda
      transform.method <- y$transform.method
      ata.output$actual <- ATA.Transform(y$seasonal.adjusted, tMethod=transform.method, tLambda=lambda)$trfmX
    }
  }else {
    if (y$seasonal.model=="x13" | y$seasonal.model=="x11"){
      seasonal.type <- y$seasonal.type
      lambda <- y$lambda
      transform.method <- y$transform.method
      ata.output$actual <- ATA.Transform(y$seasonal.adjusted, tMethod=transform.method, tLambda=lambda)$trfmX
    }else if (y$seasonal.model!="decomp" & y$seasonal.type=="M" & (y$transform.method=="BoxCox" & ty$ransform.method=="log")){
      seasonal.type <- "A"
      lambda <- 0
      transform.method <- "BoxCox"
      bcBiasAdj <- FALSE
      ata.output$actual <- ATA.Transform(y$seasonal.adjusted, tMethod=transform.method, tLambda=lambda)$trfmX        # lambda = 0 for multiplicative model
    }else {
      seasonal.type <- y$seasonal.type
      lambda <- y$lambda
      transform.method <- y$transform.method
      ata.output$actual <- ATA.Transform(y$seasonal.adjusted, tMethod=transform.method, tLambda=lambda)$trfmX
    }
  }
  if (y$is.season==FALSE & seasonal.type=="A"){
    OS_SIValue <- rep(0,times=h)
  }else if (y$is.season==FALSE & seasonal.type=="M"){
    OS_SIValue <- rep(1,times=h)
  }else if (y$is.season==TRUE){
    OS_SIValue <- rep(NA,times=h)
    for (k in 1:h){
      OS_SIValue[k] <- y$seasonal.index[freqYh[k]]
    }
  }else{
  }
  OS_SIValue <- ATA.Transform(OS_SIValue, tMethod=transform.method, tLambda=lambda)$trfmX
  ata.output$level <- ATA.Transform(y$level, tMethod=transform.method, tLambda=lambda)$trfmX
  ata.output$trend <- ATA.Transform(y$trend, tMethod=transform.method, tLambda=lambda)$trfmX
  ata.output <- AutoATA.Forecast(ata.output, hh=h, initialLevel=ata.output$initial.value)
  forecast.ata <- ata.output$forecast
  if(y$seasonal.type=="A"){
    ATA.forecast <- ATA.Inv.Transform(X=forecast.ata + OS_SIValue, tMethod=transform.method, tLambda=lambda, tbiasadj=y$bcBiasAdj, tfvar=ifelse(y$bcBiasAdj==FALSE, NULL, var(y$residuals)))
  }else {
    ATA.forecast <- ATA.Inv.Transform(X=forecast.ata * OS_SIValue, tMethod=transform.method, tLambda=lambda, tbiasadj=y$bcBiasAdj, tfvar=ifelse(y$bcBiasAdj==FALSE, NULL, var(y$residuals)))
  }
  if (negative.forecast==TRUE){
    y$forecast <- ATA.forecast
  }else {
    ATA.forecast[ATA.forecast<0] <- 0
    y$forecast <- ATA.forecast
  }
  y$h <- h
  accuracy.ata <- ATA.Accuracy(y,out.sample=out.sample)
  y$accuracy <- accuracy.ata
  y$out.sample <- ifelse(is.null(out.sample),fsample,out.sample)
  ci.output <- ATA.CI(y, ci.level)
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
  print.ata(y)
  return(y)
  print.ata(y)
  gc()
}
