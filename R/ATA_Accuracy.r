#' @title Accuracy Measures for a ATA forecast model
#' @description Returns ATA(p,q,phi) applied to \code{ata.output}.
#' Accuracy measures for a forecast model
#'
#' Returns range of summary measures of the forecast accuracy. If \code{out.sample} is
#' provided, the function measures test set forecast accuracy.
#' If \code{out.sample} is not provided, the function only produces
#' training set accuracy measures.
#'
#' The measures calculated are:
#' \itemize{
#' 		\item{MAE}	 : mean absolute error.
#' 		\item{MSE}	 : mean square error.
#' 		\item{RMSE}	 : root mean squared error.
#' 		\item{MPE}	 : mean percentage error.
#' 		\item{MAPE}	 : mean absolute percentage error.
#' 		\item{sMAPE} : symmetric mean absolute percentage error.
#' 		\item{MASE}	 : mean absolute scaled error.
#' 		\item{OWA}	 : overall weighted average of MASE and sMAPE.
#' 		\item{MdAE}	 : median absolute error.
#' 		\item{MdSE}	 : median square error.
#' 		\item{RMdSE} : root median squared error.
#' 		\item{MdPE}	 : median percentage error.
#' 		\item{MdAPE} : median absolute percentage error.
#' 		\item{sMdAPE}: symmetric median absolute percentage error.
#' }
#' @param ata.output An object of class \dQuote{\code{ata}} is required.
#' @param out.sample A numeric vector or time series of class \code{ts} or \code{msts} for out-sample.
#' @return Matrix giving forecast accuracy measures.
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
#' \code{\link{tbats}}, \code{\link{seasadj}}.
#' @references Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures
#' of forecast accuracy". \emph{International Journal of Forecasting},
#' \bold{22}(4), 679-688. Hyndman, R.J. and Athanasopoulos, G. (2014)
#' "Forecasting: principles and practice", OTexts. Section 2.5 "Evaluating
#' forecast accuracy". \url{http://www.otexts.org/fpp/2/5}.
#' @keywords ata forecast accuracy ts msts
#' @examples
#'
#' ata.fit1 <- ATA(EuStockMarkets[1:200,1],h=100)
#' ATA.Accuracy(ata.fit1)
#' ATA.Accuracy(ata.fit1,EuStockMarkets[201:300,1])
#' @export ATA.Accuracy

ATA.Accuracy <- function(ata.output, out.sample=NULL)
{
	if (class(ata.output)!="ata"){
		return("The Input must be 'ata' object. Please use ATA function to produce 'ata' object. ATA Accuracy will terminate!")
	}
	inSample <- ata.output$actual
	in_sample_fit <- ata.output$fitted
	out.sample_forecast <- ata.output$forecast
	lenX <- length(inSample)
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
		ata.error.outsample <- out.sample - out.sample_forecast
		ata.pe <- ata.error.outsample / out.sample * 100
		pre_mae_os <- abs(ata.error.outsample)
		pre_mse_os <- ata.error.outsample^2
		pre_mpe_os <- ata.pe
		pre_mape_os <- abs(ata.pe)
		pre_smape_os <- abs(out.sample - out.sample_forecast)/(abs(out.sample) + abs(out.sample_forecast)) * 200
		pre_mase_os <- outMASE(as.double(in_sample), as.double(out.sample), as.double(out.sample_forecast), as.integer(frequency(inSample)))
	}else {
		pre_mae_os <- NA
		pre_mse_os <- NA
		pre_mpe_os <- NA
		pre_mape_os <- NA
		pre_smape_os <- NA
	}
	mae <- round(mean(pre_mae, na.rm=TRUE),6)
	mse <- round(mean(pre_mse, na.rm=TRUE),6)
	rmse <- round(sqrt(mse),6)
	mpe <- round(mean(pre_mpe, na.rm=TRUE),6)
	mape <- round(mean(pre_mape, na.rm=TRUE),6)
	smape <- round(mean(pre_smape, na.rm=TRUE),6)
	mdae <- round(median(pre_mae, na.rm=TRUE),6)
	mdse <- round(median(pre_mse, na.rm=TRUE),6)
	rmdse <- round(sqrt(mdse),6)
	mdpe <- round(median(pre_mpe, na.rm=TRUE),6)
	mdape <- round(median(pre_mape, na.rm=TRUE),6)
	smdape <- round(median(pre_smape, na.rm=TRUE),6)
	mase <- round(inMASE(as.double(in_sample), as.double(in_sample_fit), as.integer(frequency(inSample))),6)
	owa <- round((mase + smape)/2, 6)
	
	stdDev_mae <- round(sqrt(var(pre_mae, na.rm=TRUE)),6)
	skew_mae <- round(colSkewness(pre_mae),6)
	kurt_mae <- round(colKurtosis(pre_mae),6)	

	stdDev_mse <- round(sqrt(var(pre_mse, na.rm=TRUE)),6)
	skew_mse <- round(colSkewness(pre_mse),6)
	kurt_mse <- round(colKurtosis(pre_mse),6)

	stdDev_mpe <- round(sqrt(var(pre_mpe, na.rm=TRUE)),6)
	skew_mpe <- round(colSkewness(pre_mpe),6)
	kurt_mpe <- round(colKurtosis(pre_mpe),6)		

	stdDev_mape <- round(sqrt(var(pre_mape, na.rm=TRUE)),6)
	skew_mape <- round(colSkewness(pre_mape),6)
	kurt_mape <- round(colKurtosis(pre_mape),6)

	stdDev_smape <- round(sqrt(var(pre_smape, na.rm=TRUE)),6)
	skew_smape <- round(colSkewness(pre_smape),6)
	kurt_smape <- round(colKurtosis(pre_smape),6)	
	
	if (!is.na(out.sample[1])){
		mae_os <- round(mean(pre_mae_os, na.rm=TRUE),6)
		mse_os <- round(mean(pre_mse_os, na.rm=TRUE),6)
		rmse_os <- round(sqrt(mse_os),6)
		mpe_os <- round(mean(pre_mpe_os, na.rm=TRUE),6)
		mape_os <- round(mean(pre_mape_os, na.rm=TRUE),6)
		smape_os <- round(mean(pre_smape_os, na.rm=TRUE),6)
		mdae_os <- round(median(pre_mae_os, na.rm=TRUE),6)
		mdse_os <- round(median(pre_mse_os, na.rm=TRUE),6)
		rmdse_os <- round(sqrt(mdse_os),6)
		mdpe_os <- round(median(pre_mpe_os, na.rm=TRUE),6)
		mdape_os <- round(median(pre_mape_os, na.rm=TRUE),6)
		smdape_os <- round(median(pre_smape_os, na.rm=TRUE),6)
		mase_os <- round(pre_mase_os, 6)
		owa_os <- round((mase_os + smape_os)/2, 6)	
	}else {
		mae_os <- NA
		mse_os <- NA
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
	pre.ata.error <- rep(NA,times=lenX)
	pre.ata.error[2:lenX] <- ata.error
	ata.error <- pre.ata.error
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
	my_list <- list("MAE"=MAE_all, "MSE"=MSE_all, "MPE"= MPE_all, "MAPE"=MAPE_all, "sMAPE"=sMAPE_all, "MASE"=MASE_all, "OWA"=OWA_all) 
	return(my_list) 
}