#' Forecasting Method for The ATAforecasting
#'
#' @description \code{ATA.Forecast} is a generic function for forecasting of the ATA Method.
#'
#' @param object An \code{ata} object is required for forecast.
#' @param h Number of periods for forecasting.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#' @param ci.level Confidence Interval levels for forecasting. Default value is 95.
#' @param negative.forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#' @param print.out Default is TRUE. If FALSE, forecast summary of ATA Method is not shown.
#'
#' @return An object of class \code{ata} and forecast values.
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
#' trainATA <-  head(touristTR, 84)
#' ata_fit <- ATA(trainATA, parPHI = 1, seasonal.test = TRUE, seasonal.model = "decomp")
#' ata_fc <- ATA.Forecast(ata_fit, h=12)
#'
#' @export
ATA.Forecast <- function(object, h=NULL, out.sample=NULL, ci.level=95, negative.forecast=TRUE, print.out = TRUE)
{
  if (class(object)!="ata"){
    stop("The Input must be 'ata' object. Please use ATA() function to produce 'ata' object.")
  }
  m <- frequency(object$actual)
  if(!is.null(out.sample)){
    if(!is.na(out.sample[1])){
      h <- length(out.sample)
      message("Forecast horizon has been set to the length of out.sample.")
    }
  }else{
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
      message(paste("Input forecast horizon has been changed with ", h))
    }
  }
  tsp_y <- tsp(object$actual)
  fsample <- forecast::msts(rep(NA,h), start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  freqYh <- cycle(fsample)
  if (object$is.season==FALSE & object$seasonal.type=="A"){
    OS_SIValue <- rep(0,times=h)
  }else if (object$is.season==FALSE & object$seasonal.type=="M"){
    OS_SIValue <- rep(1,times=h)
  }else if (object$is.season==TRUE){
    OS_SIValue <- rep(NA,times=h)
    for (k in 1:h){
      OS_SIValue[k] <- object$seasonal.index[freqYh[k]]
    }
  }else{
  }
  y <- SubATA.Forecast(object, hh = h, initialLevel = object$initial.level)
  if(object$seasonal.type=="A"){
    ATA_forecast <- y$forecast + OS_SIValue
  }else {
    ATA_forecast <- y$forecast * OS_SIValue
  }
  if (negative.forecast==TRUE){
    object$forecast <- forecast::msts(ATA_forecast, start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  }else {
    ATA_forecast[ATA_forecast<0] <- 0
    object$forecast <- forecast::msts(ATA_forecast, start = end(object$actual) + ifelse(tsp_y[3]>1, 1/tsp_y[3], 1), seasonal.periods = object$seasonal.period)
  }
  object$h <- h
  accuracy.ata <- ATA.Accuracy(object, out.sample = out.sample, print.out = FALSE)
  object$accuracy <- accuracy.ata
  object$out.sample <- ifelse(is.null(out.sample), fsample, out.sample)
  ci.output <- ATA.CI(object = object, ci.level = ci.level)
  object$ci.level <- ci.level
  if (negative.forecast==TRUE){
    object$forecast.lower <- ci.output$forecast.lower
    object$forecast.upper <- ci.output$forecast.upper
  }else {
    ci_low <- ci.output$forecast.lower
    ci_up <- ci.output$forecast.upper
    ci_low[ci_low<0] <- 0
    ci_up[ci_up<0] <- 0
    object$forecast.lower <- ci_low
    object$forecast.upper <- ci_up
  }
  attr(object, "class") <- "ata"
  if (print.out==TRUE) {
    opscipen <- options("scipen" = 100, "digits"=4)
      on.exit(options(opscipen))
    print_out <- cbind(object$forecast.lower, object$forecast, object$forecast.upper)
    colnames(print_out) <- c("lower", "forecast", "upper")
    print(print_out)
  }
  gc()
  return(object)
}
