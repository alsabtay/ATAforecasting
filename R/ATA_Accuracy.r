#' Accuracy Measures for The ATAforecasting
#'
#' @description Returns ATA(p,q,phi)(E,T,S) applied to `ata` \code{object}.
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
#' @param object An object of class \code{ata} is required.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#' @param print.out Default is TRUE. If FALSE, summary of ATA Method's accuracy measures is not shown.
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
#' trainATA <-  head(touristTR, 84)
#' testATA <- window(touristTR, start = 2015, end = 2016.917)
#' ata_fit <- ATA(trainATA, h=24, seasonal.test = TRUE, seasonal.model = "decomp")
#' ata_accuracy <- ATA.Accuracy(ata_fit, testATA)
#'
#' @export
ATA.Accuracy <- function(object, out.sample=NULL, print.out = TRUE)
{
  if (class(object)!="ata"){
    stop("The Input must be 'ata' object. Please use ATA() function to produce 'ata' object.")
  }
  train_set <- object$actual
  ts_fit <- object$fitted
  ts_fc <- object$forecast
  ts_act <- as.numeric(train_set[-1])
  ts_fit <- as.numeric(ts_fit[-1])
  if (is.null(out.sample)) {
    test_set <- NA
  }else {
    test_set <- as.numeric(out.sample)
  }
  ata.error <- ts_act - ts_fit
  ata.pe <- ata.error / ts_act * 100
  pre_mae <- abs(ata.error)
  pre_mse <- ata.error^2
  pre_mpe <- ata.pe
  pre_mape <- abs(ata.pe)
  pre_smape <- abs(ts_act - ts_fit)/(abs(ts_act) + abs(ts_fit)) * 200

  if (!is.null(test_set)){
    ata.error.os <- test_set - ts_fc
    ata.pe.os <- ata.error.os / test_set * 100
    pre_mae_os <- abs(ata.error.os)
    pre_mse_os <- ata.error.os^2
    pre_mpe_os <- ata.pe.os
    pre_mape_os <- abs(ata.pe.os)
    pre_smape_os <- abs(test_set - ts_fc)/(abs(test_set) + abs(ts_fc)) * 200
  }else {
    pre_mae_os <- NA
    pre_mse_os <- NA
    pre_mpe_os <- NA
    pre_mape_os <- NA
    pre_smape_os <- NA
  }
  np <- length(c(stats::na.omit(unlist(object$par.specs))))
  ny <- length(train_set)

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
  if (object$is.season==TRUE){
    Dec <- stats::decompose(train_set,type="multiplicative")
    des_input <- train_set/Dec$seasonal
  }else{
    des_input <- train_set
  }
  if (length(object$seasonal.period) > 1){
    naive1Accry <- NaiveSV_Accry(as.double(ts_act), as.double(object$seasonal.period), 1)
    naive2Accry <- NaiveSV_Accry(as.double(des_input), as.double(object$seasonal.period), 1)
  }else {
    naive1Accry <- NaiveSD_Accry(as.double(ts_act), as.double(object$seasonal.period), 1)
    naive2Accry <- NaiveSD_Accry(as.double(des_input), as.double(object$seasonal.period), 1)
  }
  mase <- mae / naive1Accry
  owa <- ((smape / naive2Accry) + (mase / naive2Accry)) / 2
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

  if (!is.na(test_set[1])){
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
    mase_os <- mae_os / naive1Accry
    owa_os <- ((smape_os / naive2Accry) + (mase_os / naive2Accry)) / 2
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

  RawAccuracy_is <- list("MAE"=pre_mae, "MSE"=pre_mse, "MPE"= pre_mpe, "MAPE"=pre_mape, "sMAPE"=pre_smape)
  RawAccuracy_os <- list("MAE"=pre_mae_os, "MSE"=pre_mse_os, "MPE"= pre_mpe_os, "MAPE"=pre_mape_os, "sMAPE"=pre_smape_os)
  RawAccuracy_all <- list("inSample"=RawAccuracy_is, "outSample"=RawAccuracy_os)
  MAE_is <- list("MAE"=mae, "MdAE"=mdae, "stdDev.MAE"=stdDev_mae, "skewness.MAE"=skew_mae, "kurtosis.MAE"=kurt_mae)
  MSE_is <- list("MSE"=mse, "MdSE"=mdse, "RMSE" = rmse, "RMdSE" = rmdse, "stdDev.MSE"=stdDev_mse, "skewness.MSE"=skew_mse, "kurtosis.MSE"=kurt_mse)
  MPE_is <- list("MPE"=mpe, "MdPE"=mdpe, "stdDev.MPE"=stdDev_mpe, "skewness.MPE"=skew_mpe, "kurtosis.MPE"=kurt_mpe)
  MAPE_is <- list("MAPE"=mape, "MdAPE"=mdape, "stdDev.MAPE"=stdDev_mape, "skewness.MAPE"=skew_mape, "kurtosis.MAPE"=kurt_mape)
  sMAPE_is <- list("sMAPE"=smape, "sMdAPE"=smdape, "stdDev.sMAPE"=stdDev_smape, "skewness.sMAPE"=skew_smape, "kurtosis.sMAPE"=kurt_smape)
  MASE_is <- list("MASE"=round(mase, 8), "MdASE"=NA, "stdDev.MASE"=NA, "skewness.MASE"=NA, "kurtosis.MASE"=NA)
  OWA_is <- list("OWA"=round(owa, 8), "stdDev.OWA"=NA, "skewness.OWA"=NA, "kurtosis.OWA"=NA)
  MAE_os <- list("MAE"=mae_os, "MdAE"=mdae_os)
  MSE_os <- list("MSE"=mse_os, "MdSE"=mdse_os, "RMSE" = rmse_os, "RMdSE" = rmdse_os)
  MPE_os <- list("MPE"=mpe_os, "MdPE"=mdpe_os)
  MAPE_os <- list("MAPE"=mape_os, "MdAPE"=mdape_os)
  sMAPE_os <- list("sMAPE"=smape_os, "sMdAPE"=smdape_os)
  MASE_os <- list("MASE"=round(mase_os, 8), "MdASE"=NA)
  OWA_os <- list("OWA"=round(owa_os, 8))
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
  if (print.out==TRUE) {
    opscipen <- options("scipen" = 100, "digits"=4)
      on.exit(options(opscipen))
    cat("Model Fitting Measures:","\n")
    print_out <- c(fits$sigma2, fits$loglik, MAE_all$inSample$MAE, MSE_all$inSample$MSE, MSE_all$inSample$RMSE, MPE_all$inSample$MPE, MAPE_all$inSample$MAPE, sMAPE_all$inSample$sMAPE, MASE_all$inSample$MASE, OWA_all$inSample$OWA)
    names(print_out) <- c("sigma2", "loglik", "MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(print_out)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    print_out <- c(MAE_all$outSample$MAE, MSE_all$outSample$MSE, MSE_all$outSample$RMSE, MPE_all$outSample$MPE, MAPE_all$outSample$MAPE, sMAPE_all$outSample$sMAPE, MASE_all$outSample$MASE, OWA_all$outSample$OWA)
    names(print_out) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(print_out)
    cat("\n\n")
  }
  return(my_list)
}
