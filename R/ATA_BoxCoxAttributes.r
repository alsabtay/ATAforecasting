#' The ATA.BoxCoxAttr function works with many different types of inputs.
#'
#' @param bcMethod Choose method to be used in calculating lambda. "loglik" is default. Other method is "guerrero" (Guerrero, V.M. (1993)).
#' @param bcLower Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0.
#' @param bcUpper Upper limit for possible lambda values. The upper value is limited by 5. Default value is 1.
#' @param bcBiasAdj Use adjusted back-transformed mean for Box-Cox transformations.
#' If transformed data is used to produce forecasts and fitted values, a regular back transformation will result in median forecasts.
#' If bcBiasAdj is TRUE, an adjustment will be made to produce mean forecasts and fitted values.
#' If bcBiasAdj=TRUE. Can either be the forecast variance, or a list containing the interval \code{level}, the corresponding \code{upper} and \code{lower} intervals.
#'
#' @return An object of class \code{ataattrset}.
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{\link{BoxCox}}, \code{\link{InvBoxCox}}, \code{\link{BoxCox.lambda}}
#'
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. \emph{JRSS B} \bold{26} 211--246.
#'
#' Guerrero, V.M. (1993) Time-series analysis supported by power transformations. \emph{Journal of Forecasting}, \bold{12}, 37--48.
#'
#' @export
ATA.BoxCoxAttr <- function(bcMethod = "loglik", bcLower = 0, bcUpper = 1, bcBiasAdj = FALSE)
{
  if ((bcMethod != "guerrero" & bcMethod != "loglik") | !is.character(bcMethod)){
    warning("Selected method for calculating lambda must be string. guerrero or loglik for calculating lambda.")
    bcMethod <- "loglik"
  }
  if(bcLower < -5){
    warning("Specified lower value is less than the minimum, setting bcLower=0")
    bcLower <- 0
  }else if(bcUpper > 5){
    warning("Specified upper value is larger than the maximum, setting bcUpper=1")
    bcUpper <- 1
  }
  if (!is.logical(bcBiasAdj)) {
    warning("bcBiasAdj information not found, defaulting to FALSE.")
    biasadj <- FALSE
  }
  mylist <- list("bcMethod"=bcMethod, "bcLower"=bcLower, "bcUpper"=bcUpper, "bcBiasAdj"=bcBiasAdj)
  attr(mylist, "class") <- "ataattrset"
  return(mylist)
}
