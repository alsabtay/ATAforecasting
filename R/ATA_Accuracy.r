#' Accuracy Measures for The ATAforecasting
#'
#' @description Returns ATA(p,q,phi) applied to ATA \code{object}.
#' Accuracy measures for a forecast model
#' Returns range of summary measures of the forecast accuracy. If \code{out.sample} is
#' provided, the function measures test set forecast accuracy.
#' If \code{out.sample} is not provided, the function only produces
#' training set accuracy measures.
#' The measures calculated are:
#' \itemize{
#'		 \item{lik}		: maximum likelihood functions
#'		 \item{sigma}	: residual variance.
#'		 \item{MAE}		: mean absolute error.
#'		 \item{MSE}		: mean square error.
#'		 \item{AMSE}	: Average MSE over first `nmse` forecast horizons.
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
#'
#' @param object An object of class \code{ATA} is required.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#'
#' @return Matrix giving forecast accuracy measures.
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{forecast}, \code{stlplus}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}}, \code{tbats}, \code{seasadj}.
#'
#' @references
#'
#' #'\insertRef{hyndmanandkoehler2006}{ATAforecasting}
#'
#' #'\insertRef{hyndman2019forecasting}{ATAforecasting}
#'
#'
#' @keywords Ata forecast accuracy ts msts
#'
#' @importFrom stats median var
#' @importFrom timeSeries colKurtosis colSkewness
#' @importFrom Rdpack reprompt
#'
#' @examples
#' demoATA <- window(touristTR, start = 2008, end = 2014.917)
#' ata.fit <- ATA(demoATA, h=18, seasonal.test = TRUE, seasonal.model = "stl")
#' ata.accuracy <- ATA.Accuracy(ata.fit, tail(touristTR,18))
#'
#' @export
ATA.Accuracy <- function(object, out.sample=NULL)
{
  ata.output <- object
  if (class(ata.output)!="ATA"){
    return("The Input must be 'ATA' object. Please use ATA function to produce 'ATA' object. ATA Accuracy will terminate!")
  }
  inSample <- ata.output$actual
  in_sample_fit <- ata.output$fitted
  out.sample_forecast <- ata.output$forecast
  in_sample <- as.numeric(inSample[-1])
  in_sample_fit <- as.numeric(in_sample_fit[-1])
  if (!is.null(out.sample)){
    out.sample <- as.numeric(out.sample)
  }else {
    out.sample <- NA
  }
  ata.error <- in_sample - in_sample_fit
  ata.pe <- ata.error / in_sample * 100
  pre_mae <- abs(ata.error)
  pre_mse <- ata.error^2
  pre_mpe <- ata.pe
  pre_mape <- abs(ata.pe)
  pre_smape <- abs(in_sample - in_sample_fit)/(abs(in_sample) + abs(in_sample_fit)) * 200

  if (!is.null(out.sample)){
    ata.error.os <- out.sample - out.sample_forecast
    ata.pe.os <- ata.error.os / out.sample * 100
    pre_mae_os <- abs(ata.error.os)
    pre_mse_os <- ata.error.os^2
    pre_mpe_os <- ata.pe.os
    pre_mape_os <- abs(ata.pe.os)
    pre_smape_os <- abs(out.sample - out.sample_forecast)/(abs(out.sample) + abs(out.sample_forecast)) * 200
    pre_mase_os <- outMASE(as.double(in_sample), as.double(out.sample), as.double(out.sample_forecast), as.integer(frequency(inSample)))
  }else {
    pre_mae_os <- NA
    pre_mse_os <- NA
    pre_mpe_os <- NA
    pre_mape_os <- NA
    pre_smape_os <- NA
  }
  np <- length(c(stats::na.omit(unlist(ata.output$par.specs))))
  ny <- length(ata.output$actual)

  mae <- round(mean(pre_mae, na.rm=TRUE),8)
  mse <- round(mean(pre_mse, na.rm=TRUE),8)
  lik <- ny * round(log(sum(pre_mse, na.rm=TRUE)),8)
  rmse <- round(sqrt(mse),8)
  mpe <- round(mean(pre_mpe, na.rm=TRUE),8)
  mape <- round(mean(pre_mape, na.rm=TRUE),8)
  smape <- round(mean(pre_smape, na.rm=TRUE),8)
  mdae <- round(median(pre_mae, na.rm=TRUE),8)
  mdse <- round(median(pre_mse, na.rm=TRUE),8)
  rmdse <- round(sqrt(mdse),8)
  mdpe <- round(median(pre_mpe, na.rm=TRUE),8)
  mdape <- round(median(pre_mape, na.rm=TRUE),8)
  smdape <- round(median(pre_smape, na.rm=TRUE),8)
  mase <- round(inMASE(as.double(in_sample), as.double(in_sample_fit), as.integer(frequency(inSample))),8)
  naiveAccry <- round(NaiveSD(as.double(in_sample), as.integer(frequency(inSample))),8)
  owa <- round(((mase/naiveAccry) + (smape/naiveAccry))/2, 8)

  aic <- lik + 2 * np
  bic <- lik + log(ny) * np
  aicc <- aic + 2 * np * (np + 1) / (ny - np - 1)

  stdDev_mae <- round(sqrt(var(pre_mae, na.rm=TRUE)),8)
  skew_mae <- round(timeSeries::colSkewness(pre_mae),8)
  kurt_mae <- round(timeSeries::colKurtosis(pre_mae),8)

  stdDev_mse <- round(sqrt(var(pre_mse, na.rm=TRUE)),8)
  skew_mse <- round(timeSeries::colSkewness(pre_mse),8)
  kurt_mse <- round(timeSeries::colKurtosis(pre_mse),8)

  stdDev_mpe <- round(sqrt(var(pre_mpe, na.rm=TRUE)),8)
  skew_mpe <- round(timeSeries::colSkewness(pre_mpe),8)
  kurt_mpe <- round(timeSeries::colKurtosis(pre_mpe),8)

  stdDev_mape <- round(sqrt(var(pre_mape, na.rm=TRUE)),8)
  skew_mape <- round(timeSeries::colSkewness(pre_mape),8)
  kurt_mape <- round(timeSeries::colKurtosis(pre_mape),8)

  stdDev_smape <- round(sqrt(var(pre_smape, na.rm=TRUE)),8)
  skew_smape <- round(timeSeries::colSkewness(pre_smape),8)
  kurt_smape <- round(timeSeries::colKurtosis(pre_smape),8)

  if (!is.na(out.sample[1])){
    mae_os <- round(mean(pre_mae_os, na.rm=TRUE),8)
    mse_os <- round(mean(pre_mse_os, na.rm=TRUE),8)
    rmse_os <- round(sqrt(mse_os),8)
    mpe_os <- round(mean(pre_mpe_os, na.rm=TRUE),8)
    mape_os <- round(mean(pre_mape_os, na.rm=TRUE),8)
    smape_os <- round(mean(pre_smape_os, na.rm=TRUE),8)
    mdae_os <- round(median(pre_mae_os, na.rm=TRUE),8)
    mdse_os <- round(median(pre_mse_os, na.rm=TRUE),8)
    rmdse_os <- round(sqrt(mdse_os),8)
    mdpe_os <- round(median(pre_mpe_os, na.rm=TRUE),8)
    mdape_os <- round(median(pre_mape_os, na.rm=TRUE),8)
    smdape_os <- round(median(pre_smape_os, na.rm=TRUE),8)
    mase_os <- round(pre_mase_os,8)
    owa_os <- round(((mase_os/naiveAccry) + (smape_os/naiveAccry))/2,8)
  }else {
    mae_os <- NA
    mse_os <- NA
    lik_os <- NA
    rmse_os <- NA
    mpe_os <- NA
    mape_os <- NA
    smape_os <- NA
    mdae_os <- NA
    mdse_os <- NA
    rmdse_os <- NA
    mdpe_os <- NA
    mdape_os <- NA
    smdape_os <- NA
    mase_os <- NA
    owa_os <- NA
  }

  RawAccuracy_is <- list("MAE"=pre_mae, "MSE"=pre_mse, "MPE"= pre_mpe, "MAPE"=pre_mape, "sMAPE"=pre_smape, "MASE" = (pre_mae/naiveAccry))
  RawAccuracy_os <- list("MAE"=pre_mae_os, "MSE"=pre_mse_os, "MPE"= pre_mpe_os, "MAPE"=pre_mape_os, "sMAPE"=pre_smape_os, "MASE" = (pre_mae_os/naiveAccry))
  RawAccuracy_all <- list("inSample"=RawAccuracy_is, "outSample"=RawAccuracy_os)
  MAE_is <- list("MAE"=mae, "MdAE"=mdae, "stdDev.MAE"=stdDev_mae, "skewness.MAE"=skew_mae, "kurtosis.MAE"=kurt_mae)
  MSE_is <- list("MSE"=mse, "MdSE"=mdse, "RMSE" = rmse, "RMdSE" = rmdse, "stdDev.MSE"=stdDev_mse, "skewness.MSE"=skew_mse, "kurtosis.MSE"=kurt_mse)
  MPE_is <- list("MPE"=mpe, "MdPE"=mdpe, "stdDev.MPE"=stdDev_mpe, "skewness.MPE"=skew_mpe, "kurtosis.MPE"=kurt_mpe)
  MAPE_is <- list("MAPE"=mape, "MdAPE"=mdape, "stdDev.MAPE"=stdDev_mape, "skewness.MAPE"=skew_mape, "kurtosis.MAPE"=kurt_mape)
  sMAPE_is <- list("sMAPE"=smape, "sMdAPE"=smdape, "stdDev.sMAPE"=stdDev_smape, "skewness.sMAPE"=skew_smape, "kurtosis.sMAPE"=kurt_smape)
  MASE_is <- list("MASE"=mase, "MdASE"=NA, "stdDev.MASE"=NA, "skewness.MASE"=NA, "kurtosis.MASE"=NA)
  OWA_is <- list("OWA"=owa, "stdDev.OWA"=NA, "skewness.OWA"=NA, "kurtosis.OWA"=NA)
  MAE_os <- list("MAE"=mae_os, "MdAE"=mdae_os)
  MSE_os <- list("MSE"=mse_os, "MdSE"=mdse_os, "RMSE" = rmse_os, "RMdSE" = rmdse_os)
  MPE_os <- list("MPE"=mpe_os, "MdPE"=mdpe_os)
  MAPE_os <- list("MAPE"=mape_os, "MdAPE"=mdape_os)
  sMAPE_os <- list("sMAPE"=smape_os, "sMdAPE"=smdape_os)
  MASE_os <- list("MASE"=mase_os, "MdASE"=NA)
  OWA_os <- list("OWA"=owa_os)
  MAE_all <- list("inSample"=MAE_is, "outSample"=MAE_os)
  MSE_all <- list("inSample"=MSE_is, "outSample"=MSE_os)
  MPE_all <- list("inSample"=MPE_is, "outSample"=MPE_os)
  MAPE_all <- list("inSample"=MAPE_is, "outSample"=MAPE_os)
  sMAPE_all <- list("inSample"=sMAPE_is, "outSample"=sMAPE_os)
  MASE_all <- list("inSample"=MASE_is, "outSample"=MASE_os)
  OWA_all <- list("inSample"=OWA_is, "outSample"=OWA_os)
  fits <- list("sigma2" = round(sum(pre_mse, na.rm=TRUE),8) / (ny - np),
            "loglik" = -0.5 * lik,
            "AIC" = aic,
            "AICc" = aicc,
            "BIC" = bic,
            "MSE" = mse,
            "MAE" = mae,
            "sMAPE" = smape,
            "MASE" = mase,
            "OWA" = owa)
  my_list <- list("MAE"=MAE_all, "MSE"=MSE_all, "MPE"= MPE_all, "MAPE"=MAPE_all, "sMAPE"=sMAPE_all, "MASE"=MASE_all, "OWA"=OWA_all, "RawAccuracy"=RawAccuracy_all, "fits"=fits)
  return(my_list)
}
