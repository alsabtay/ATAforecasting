#' The ATA.BoxCoxAttr function works with many different types of inputs.
#'
#' @param bcMethod Choose method to be used in calculating lambda. "guerrero" (Guerrero, V.M. (1993) is default. Other method is "loglik").
#' @param bcLower Lower limit for possible lambda values. The lower value is limited by -5. Default value is 0.
#' @param bcUpper Upper limit for possible lambda values. The upper value is limited by 5. Default value is 5.
#' @param bcBiasAdj Use adjusted back-transformed mean for Box-Cox transformations.
#' If transformed data is used to produce forecasts and fitted values, a regular back transformation will result in median forecasts.
#' If bcBiasAdj is TRUE, an adjustment will be made to produce mean forecasts and fitted values.
#' If bcBiasAdj=TRUE. Can either be the forecast variance, or a list containing the interval \code{level}, the corresponding \code{upper} and \code{lower} intervals.
#'
#' @return An object of class \code{ataoptim}.
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#'
#' @seealso \code{\link{BoxCox}}, \code{\link{InvBoxCox}}, \code{\link{BoxCox.lambda}}
#'
#' @references
#'
#' #'\insertRef{boxcox1964}{ATAforecasting}
#'
#' #'\insertRef{guerrero1993}{ATAforecasting}
#'
#'
#'
#' @export
ATA.BoxCoxAttr <- function(bcMethod = "guerrero", bcLower = 0, bcUpper = 5, bcBiasAdj = FALSE)
{
  if ((bcMethod != "guerrero" & bcMethod != "loglik") | !is.character(bcMethod)){
    warning("Selected method for calculating lambda must be string. guerrero or loglik for calculating lambda.")
    bcMethod <- "guerrero"
  }
  if(bcLower < 0){
    warning("Specified lower value is less than the minimum, setting bcLower=0")
    bcLower <- 0
  }else if(bcUpper > 5){
    warning("Specified upper value is larger than the maximum, setting bcUpper=5")
    bcUpper <- 5
  }
  if (!is.logical(bcBiasAdj)) {
    warning("bcBiasAdj information not found, defaulting to FALSE.")
    biasadj <- FALSE
  }
  mylist <- list("bcMethod"=bcMethod, "bcLower"=bcLower, "bcUpper"=bcUpper, "bcBiasAdj"=bcBiasAdj)
  attr(mylist, "class") <- "ataoptim"
  return(mylist)
}
