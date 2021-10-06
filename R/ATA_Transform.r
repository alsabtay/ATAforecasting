#' Transformation Techniques for The ATAforecasting
#'
#' @description The function provides the applicability of different types of transformation techniques for the data to which the Ata method will be applied.
#' The \code{ATA.Transform} function works with many different types of inputs.
#'
#' @param X a numeric vector or time series of class \code{ts} or \code{msts} for in-sample.
#' @param tMethod Box-Cox power transformation family is consist of "Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog",
#' "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog" in ATAforecasting package. If the transformation process needs shift parameter,
#' \code{ATA.Transform} will calculate required shift parameter automatically.
#' @param tLambda Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.
#' @param tShift Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated.
#' @param bcMethod Choose method to be used in calculating lambda. "loglik" is default. Other method is "guerrero" (Guerrero, V.M. (1993)).
#' @param bcLower Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0.
#' @param bcUpper Upper limit for possible lambda values. The upper value is limited by 5. Default value is 1.
#'
#' @return A list object consists of transformation parameters and transformed data.
#' \code{ATA.Transform} is a list containing at least the following elements:
#' \itemize{
#'		 \item{trfmX}   : Transformed data
#'		 \item{tLambda} : Box-Cox power transformation family parameter
#'		 \item{tShift}  : Box-Cox power transformation family shifting parameter
#'}
#'
#' @references
#'
#' #'\insertRef{tukey1957}{ATAforecasting}
#'
#' #'\insertRef{boxcox1964}{ATAforecasting}
#'
#' #'\insertRef{manly1976}{ATAforecasting}
#'
#' #'\insertRef{johndraper1980}{ATAforecasting}
#'
#' #'\insertRef{bickeldoksum1982}{ATAforecasting}
#'
#' #'\insertRef{sakia1992}{ATAforecasting}
#'
#' #'\insertRef{guerrero1993}{ATAforecasting}
#'
#' #'\insertRef{yeojohn2000}{ATAforecasting}
#'
#' #'\insertRef{glog2002}{ATAforecasting}
#'
#' #'\insertRef{neglog2005}{ATAforecasting}
#'
#' #'\insertRef{yang2006}{ATAforecasting}
#'
#' #'\insertRef{gpower2013}{ATAforecasting}
#'
#' @keywords Ata Bickel--Doksum Box--Cox dual glog gpower Guerrero Manly neglog transformation Yeo--Johnson
#'
#' @importFrom forecast BoxCox.lambda
#' @importFrom Rdpack reprompt
#'
#' @export
ATA.Transform <- function(X
                          , tMethod = c("Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog")
                          , tLambda
                          , tShift = 0
                          , bcMethod = c("loglik", "guerrero")
                          , bcLower = 0
                          , bcUpper = 5)
{
	if (is.null(tMethod)){
		my_list <- list("trfmX" = X, "tLambda" = tLambda, "tShift" = tShift)
	}else {
		if (is.null(tLambda)){
		  if (bcMethod == "guerrero") {
			tLambda <- forecast::BoxCox.lambda(X, method = "guerrero", lower=bcLower, upper=bcUpper)
		  } else {
			tLambda <- forecast::BoxCox.lambda(X, method = "loglik", lower=bcLower, upper=bcUpper)
		  }
		}
		out_list <- SubATA.Transform(X, tMethod = tMethod, tType = "Vanilla", tLambda = tLambda, tShift = tShift)
		my_list <- list("trfmX" = out_list$tX, "tLambda" = out_list$tLambda, "tShift" = out_list$tShift)
	}
	return(my_list)
}

#' Back Transformation Techniques for The ATAforecasting
#'
#' @description The function provides the applicability of different types of back transformation techniques for the transformed data to which the Ata method will be applied.
#' The \code{ATA.BackTransform} function works with many different types of inputs.
#' @param X a numeric vector or time series of class \code{ts} or \code{msts} for in-sample.
#' @param tMethod Box-Cox power transformation family is consist of "Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog",
#' "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog" in ATAforecasting package.
#' @param tLambda Box-Cox power transformation family parameter. If NULL, data transformed before model is estimated.
#' @param tShift Box-Cox power transformation family shifting parameter. If NULL, data transformed before model is estimated.
#' @param tbiasadj Use adjusted back-transformed mean for Box-Cox transformations using \code{forecast::BoxCox}. If transformed data is used to produce forecasts and fitted values,
#' a regular back transformation will result in median forecasts. If tbiasadj is TRUE, an adjustment will be made to produce mean forecasts and fitted values.
#' @param tfvar Optional parameter required if tbiasadj=TRUE. Can either be the forecast variance, or a list containing the interval \code{level}, and the
#' corresponding \code{upper} and \code{lower} intervals.
#'
#' @return A list object consists of transformation parameters and transformed data.
#' \code{ATA.Transform} is a list containing at least the following elements:
#' \itemize{
#'		 \item{trfmX}   : Transformed data
#'		 \item{tLambda} : Box-Cox power transformation family parameter
#'		 \item{tShift}  : Box-Cox power transformation family shifting parameter
#'}
#'
#' @references
#'
#' #'\insertRef{tukey1957}{ATAforecasting}
#'
#' #'\insertRef{boxcox1964}{ATAforecasting}
#'
#' #'\insertRef{manly1976}{ATAforecasting}
#'
#' #'\insertRef{johndraper1980}{ATAforecasting}
#'
#' #'\insertRef{bickeldoksum1982}{ATAforecasting}
#'
#' #'\insertRef{sakia1992}{ATAforecasting}
#'
#' #'\insertRef{guerrero1993}{ATAforecasting}
#'
#' #'\insertRef{yeojohn2000}{ATAforecasting}
#'
#' #'\insertRef{glog2002}{ATAforecasting}
#'
#' #'\insertRef{neglog2005}{ATAforecasting}
#'
#' #'\insertRef{yang2006}{ATAforecasting}
#'
#' #'\insertRef{gpower2013}{ATAforecasting}
#'
#' @keywords Ata Bickel--Doksum Box--Cox dual glog gpower Guerrero Manly neglog transformation Yeo--Johnson
#'
#'
#' @export
ATA.BackTransform <- function(X, tMethod, tLambda, tShift, tbiasadj=FALSE, tfvar=NULL)
{
	if (is.null(tMethod)){
		trfmX <- X
	}else {
		out_list <- SubATA.Transform(tX = X, tMethod = tMethod, tType = "Back", tLambda = tLambda, tShift = tShift)
		trfmX <- out_list$tX
	}
	return(trfmX)
}
