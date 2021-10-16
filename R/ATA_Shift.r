#' Lag/Lead (Shift) Function for Univariate Series
#'
#' @param x given vector
#' @param shift_by lag or lead parameter
#' @param fill a value to be used to fill the rows
#'
#' @return Generating a lag/lead variables
#' @importFrom utils head tail
#' @export
ATA.Shift <- function(x, shift_by, fill = NA){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  if (length(shift_by) > 1)
    return(sapply(shift_by, ATA.Shift, x=x))
  out <- NULL
  abs_shift_by <- abs(shift_by)
  if (shift_by > 0 )
    out <- c(tail(x,-abs_shift_by),rep(fill,abs_shift_by))
  else if (shift_by < 0 )
    out <- c(rep(fill,abs_shift_by), utils::head(x,-abs_shift_by))
  else
    out <- x
  return(out)
}


#' Lag/Lead (Shift) Function for Multivariate Series
#'
#' @param X given matrice
#' @param direction direction of shifting. Default is "down".
#' @param shift_by number of rows to be shifed upwards/downwards
#' @param fill a value to be used to fill the rows
#'
#' @return Generating a lag/lead matrice
#'
#' @export
ATA.Shift_Mat <- function(X, direction = "down", shift_by = 1, fill = NA)
{
  if (direction != "down" & direction != "up"){
    stop("'direction' is not a 'down' or 'up'")
  }
  if (!is.matrix(X)) {
    stop("X is not a matrix")
  }
  if (!is.numeric(X)) {
    stop("X is not a numeric matrix")
  }
  if (shift_by < 0)
    stop("'rows' is not positive")
  if (shift_by != trunc(shift_by))
    stop("'rows' is not an integer")
  if (direction == "down"){
    if (shift_by > 0)
      return(ATA.Shift_Mat(rbind(rep(fill, ncol(X)), X[1:nrow(X)-1,]), direction = "down", shift_by - 1, fill))
  }else {
    if (shift_by > 0)
      return(ATA.Shift_Mat(rbind(X[2:nrow(X),], rep(fill, ncol(X))), direction = "up", shift_by - 1, fill))
  }
  return(X)
}
