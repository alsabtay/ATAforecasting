#' Lag/Lead (Shift) Function for Univariate Series
#'
#' @param x given variable(s)
#' @param shift_by lag or lead parameter
#'
#' @return Generating a lag/lead variables
#' @importFrom utils head tail
#' @export
#'
ATA.Shift <- function(x, shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  if (length(shift_by) > 1)
    return(sapply(shift_by, ATA.Shift, x=x))
  out <- NULL
  abs_shift_by <- abs(shift_by)
  if (shift_by > 0 )
    out <- c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out <- c(rep(NA,abs_shift_by), utils::head(x,-abs_shift_by))
  else
    out <- x
  return(out)
}
