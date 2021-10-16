#' Seasonal Decomposition for The ATAforecasting
#'
#' @description Automatic seasonal decomposition for ATA Method is called \code{ATA.Decomposition} function in ATAforecasting package.
#' The function returns seasonally adjusted data constructed by removing the seasonal component. The methodology is fully automatic.
#' The \code{ATA.Decomposition} function works with many different types of inputs.
#' @param input It must be \code{ts} or \code{msts} or \code{numeric} object. if it is \code{numeric} object, \code{findPeriod} must be 1 or 2 or 3 or 4. if it is \code{msts} object, \code{findPeriod} must be 3 or 4.
#' @param s.model A string identifying method for seasonal decomposition. If NULL, "decomp" method is default. c("none", "decomp", "stl", "stlplus", "tbats", "stR") phrases of methods denote.
#' \itemize{
#'		 \item{none}	: seasonal decomposition is not required.
#'		 \item{decomp} 	: classical seasonal decomposition. If \code{decomp}, the \code{stats} package will be used.
#'		 \item{stl}		: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stl}, the \code{stats} and \code{forecast} packages will be used. Multiple seasonal periods are allowed.
#'		 \item{stlplus}	: seasonal-trend decomposition procedure based on loess developed by Cleveland et al. (1990). If \code{stlplus}, the \code{stlplus} package will be used.
#'		 \item{tbats}   : exponential smoothing state space model with Box--Cox transformation, ARMA errors, trend and seasonal components.
#' 					  	  as described in De Livera, Hyndman & Snyder (2011). Parallel processing is used by default to speed up the computations. If \code{tbats}, the \code{forecast} package will be used. Multiple seasonal periods are allowed.
#'		 \item{stR}    	: seasonal-trend decomposition procedure based on regression developed by Dokumentov and Hyndman (2015). If \code{stR}, the \code{stR} package will be used. Multiple seasonal periods are allowed.
#'		 \item{x13}    	: seasonal-trend decomposition procedure based on X13ARIMA/SEATS. If \code{x13}, the \code{seasonal} package will be used.
#'		 \item{x11}    	: seasonal-trend decomposition procedure based on X11. If \code{x11}, the \code{seasonal} package will be used.
#' }
#' @param s.type A one-character string identifying method for the seasonal component framework. If NULL, "M" is default. The letter "A" for additive model, the letter "M" for multiplicative model.
#' @param s.frequency Value(s) of seasonal periodicity. If \code{s.frequency} is not integer, \code{X} must be \code{msts} time series object. c(s1,s2,s3,...) for multiple period. If \code{X} has multiple periodicity, "tbats" or "stR" seasonal model have to be selected.
#' @param seas_attr_set Assign from \code{ATA.SeasAttr} function. Attributes set for unit root and seasonality tests.
#' For example: period of the input data which have one seasonal pattern --> 12 for monthly / 4 for quarterly / 7 for daily / 5 for business days. periods of the input data which have complex/multiple seasonal patterns --> c(7,354.37,365.25).
#'
#' @return Seasonal components of the univariate time series.
#' \code{ATA.Decomposition} is a list containing at least the following elements:
#' \item{AdjustedX}{Deseasonalized data}
#' \item{SeasIndex}{Particular weights of seasonality given cycle/frequency}
#' \item{SeasActual}{Seasonality given original data}
#' \item{SeasType}{Seasonal decomposition technique}
#'
#' @author Ali Sabri Taylan and Hanife Taylan Selamlar
#' @seealso \code{\link[stats]{stl}}, \code{\link[stats]{decompose}}, \code{\link[seasonal]{seas}},
#' \code{\link[forecast]{tbats}}, \code{\link{stlplus}}, \code{\link[stR]{AutoSTR}}.
#'
#' @keywords Ata seasonal decomposition forecast accuracy ts msts mstl
#'
#' @references
#'
#' #'\insertRef{shishkin1967}{ATAforecasting}
#'
#' #'\insertRef{dagum1988}{ATAforecasting}
#'
#' #'\insertRef{cleveland1990stl}{ATAforecasting}
#'
#' #'\insertRef{hafen2010local}{ATAforecasting}
#'
#' #'\insertRef{delivera2011}{ATAforecasting}
#'
#' #'\insertRef{dokumentov2015}{ATAforecasting}
#'
#' #'\insertRef{dokumentov2020str}{ATAforecasting}
#'
#' #'\insertRef{monsell2003toward}{ATAforecasting}
#'
#' #'\insertRef{monsell2007x}{ATAforecasting}
#'
#' #'\insertRef{artseasonal2018}{ATAforecasting}
#'
#'
#' @importFrom forecast mstl msts tbats tbats.components
#' @importFrom stats cycle decompose frequency ts tsp tsp<- stl
#' @importFrom stlplus stlplus
#' @importFrom stR AutoSTR
#' @importFrom seasonal seas series udg
#' @importFrom Rdpack reprompt
#'
#' @export
#'
ATA.Decomposition <- function(input, s.model, s.type, s.frequency, seas_attr_set)
{
  tsp_input <- tsp(input)
  last_seas_type <- s.type
  if (s.model == "none" | min(s.frequency)==1){
    if (s.type=="A"){
      adjX <- input
      SeasActual <- rep(0,times=length(input))
      SeasActual <- ts(SeasActual, frequency = tsp_input[3], start = tsp_input[1])
      s.frequency <- frequency(input)
      SeasIndex <- rep(0,times=s.frequency)
    }else {
      adjX <- input
      SeasActual <- rep(1,times=length(input))
      SeasActual <- ts(SeasActual, frequency = tsp_input[3], start = tsp_input[1])
      s.frequency <- frequency(input)
      SeasIndex <- rep(1,times=s.frequency)
    }
  }else {
    if (class(input)[1]!="ts" & class(input)[1]!="msts"){
      return("The data set must be time series object (ts or msts) ATA Method was terminated!")
    }
    input <- forecast::msts(input, start=tsp_input[1], seasonal.periods = s.frequency)
    tsp_input <- tsp(input)
    if (s.model=="decomp"){									  	# Do classical decomposition
      if (s.type=="A"){
        desX <- stats::decompose(input, type = c("additive"))
        adjX <- forecast::seasadj(desX)
        SeasActual <- desX$seasonal
        SeasIndex <- rep(NA,times=s.frequency)
        for (s in 1:s.frequency){
          SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
        }
      }else {
        desX <- stats::decompose(input, type = c("multiplicative"))
        adjX <- forecast::seasadj(desX)
        SeasActual <- desX$seasonal
        SeasIndex <- rep(NA,times=s.frequency)
        for (s in 1:s.frequency){
          SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
        }
      }
    }else if (s.model=="stl"){									# Do STL decomposition
      if (length(s.frequency)==1){
        stldesX <- stats::stl(input, s.window = "per", robust=TRUE)
        adjX <- forecast::seasadj(stldesX)
        SeasActual <- forecast::seasonal(stldesX)
        SeasIndex <- rep(NA,times=s.frequency)
        for (s in 1:s.frequency){
          SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
        }
      }else {
        stldesX <- forecast::mstl(input, lambda = NULL, s.window = "per")
        nameCol <- colnames(stldesX)
        nameCol <- grep('Season', nameCol, value=TRUE)
        if (length(nameCol)==0){
          if (s.type=="A"){
            adjX <- input
            SeasActual <- forecast::msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
            SeasIndex <- rep(0,times=max(s.frequency))
          }else {
            adjX <- input
            SeasActual <- forecast::msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
            SeasIndex <- rep(1,times=max(s.frequency))
          }
        }else {
          adjX <- forecast::seasadj(stldesX)
          if (length(s.frequency)==1){
            SeasActual <- stldesX[,nameCol]
            SeasIndex <- rep(NA,times=s.frequency)
            for (s in 1:s.frequency){
              SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
            }
          }else {
            SeasActual <- rowSums(stldesX[,nameCol],na.rm=TRUE)
            SeasActual <- forecast::msts(SeasActual, start=tsp_input[1], seasonal.periods = tsp_input[3])
            SeasIndex <- rep(NA,times=max(s.frequency))
            for (s in 1:max(s.frequency)){
              SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
            }
          }
        }
      }
    }else if (s.model=="stlplus"){								# Do STLPlus decomposition
	  stlplusdesX <- stlplus::stlplus(input, s.window = "per", robust=TRUE)
      adjX <- input - stlplusdesX$data$seasonal
      SeasActual <- stlplusdesX$data$seasonal
      SeasActual <- forecast::msts(SeasActual, start=tsp_input[1], seasonal.periods = s.frequency)
      SeasIndex <- rep(NA,times=s.frequency)
      for (s in 1:s.frequency){
        SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
      }
    }else if (s.model=="stR"){									# Do stR decomposition
	  if (length(input)>1600){
        stRdesX <- stR::AutoSTR(input)
      }else {
        stRdesX <- stR::AutoSTR(input, robust=TRUE)
      }
      stRcomp <- stR_components(stRdesX)
      nameCol <- colnames(stRcomp)
      nameCol <- grep('Seasonal', nameCol, value=TRUE)
      if (length(nameCol)==0){
        if (s.type=="A"){
          adjX <- input
          SeasActual <- forecast::msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(0,times=max(s.frequency))
        }else {
          adjX <- input
          SeasActual <- forecast::msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(1,times=max(s.frequency))
        }
      }else {
        adjX <- stR_seasadj(stRdesX)
        if (length(s.frequency)==1){
          SeasActual <- stRcomp[,nameCol]
          SeasIndex <- rep(NA,times=s.frequency)
          for (s in 1:s.frequency){
            SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
          }
        }else {
          SeasActual <- rowSums(stRcomp[,nameCol],na.rm=TRUE)
          SeasActual <- forecast::msts(SeasActual, start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(NA,times=max(s.frequency))
          for (s in 1:max(s.frequency)){
            SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
          }
        }
      }
    }else if (s.model=="tbats"){								# Do tbats decomposition
      tbatsdesX <- forecast::tbats(input, use.box.cox = FALSE)
      tbatscomp <- forecast::tbats.components(tbatsdesX)
      nameCol <- colnames(tbatscomp)
      nameCol <- grep('season', nameCol, value=TRUE)
      if (length(nameCol)==0){
        if (s.type=="A"){
          adjX <- input
          SeasActual <- forecast::msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(0,times=max(s.frequency))
        }else {
          adjX <- input
          SeasActual <- forecast::msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(1,times=max(s.frequency))
        }
      }else {
        adjX <- forecast::seasadj(tbatsdesX)
        if (length(s.frequency)==1){
          SeasActual <- tbatscomp[,nameCol]
          SeasIndex <- rep(NA,times=s.frequency)
          for (s in 1:s.frequency){
            SeasIndex[s] <- as.numeric(SeasActual[cycle(SeasActual)==s][1])
          }
        }else {
          SeasActual <- rowSums(tbatscomp[,nameCol],na.rm=TRUE)
          SeasActual <- forecast::msts(SeasActual, start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(NA,times=max(s.frequency))
          for (s in 1:max(s.frequency)){
            SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
          }
        }
      }
    }else if (s.model=="x13"){									# Do X13ARIMA/SEATS decomposition
	  x13desX <- seasonal::seas(input, transform.function="none", estimate.maxiter=seas_attr_set$x13.estimate.maxiter, estimate.tol=seas_attr_set$x13.estimate.tol)
      SeasActual <- seasonal::series(x13desX,"seats.adjustfac")
      ifelse(seasonal::udg(x13desX, stats = "finmode")=="additive", s.type <- "A", s.type <- "M")
      if (is.null(SeasActual)) {
        if (s.type=="A"){
          adjX <- input
          SeasActual <- forecast::msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(0,times=max(s.frequency))
        }else {
          adjX <- input
          SeasActual <- forecast::msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(1,times=max(s.frequency))
        }
      }else {
        adjX <- seasonal::series(x13desX,"seats.seasonaladj")
        SeasIndex <- rep(NA,times=s.frequency)
        for (s in 1:s.frequency){
          SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
        }
      }
    }else if (s.model=="x11"){									# Do X13ARIMA/SEATS X11 decomposition
	  x11desX <- seasonal::seas(input, x11 = "", transform.function="none", estimate.maxiter=seas_attr_set$x11.estimate.maxiter, estimate.tol=seas_attr_set$x11.estimate.tol)
      SeasActual <- seasonal::series(x11desX,"x11.adjustfac")
      ifelse(seasonal::udg(x11desX, stats = "finmode")=="additive", s.type <- "A", s.type <- "M")
      if (is.null(SeasActual)) {
        if (s.type=="A"){
          adjX <- input
          SeasActual <- forecast::msts(rep(0,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(0,times=max(s.frequency))
        }else {
          adjX <- input
          SeasActual <- forecast::msts(rep(1,times=length(input)), start=tsp_input[1], seasonal.periods = tsp_input[3])
          SeasIndex <- rep(1,times=max(s.frequency))
        }
      }else {
        adjX <- seasonal::series(x11desX,"x11.seasadj")
        SeasIndex <- rep(NA,times=s.frequency)
        for (s in 1:s.frequency){
          SeasIndex[s] <- as.numeric(mean(SeasActual[cycle(SeasActual)==s]))
        }
      }
    }else {
    }
  }
  my_list <- list("AdjustedX" = adjX, "SeasIndex" = SeasIndex, "SeasActual" = SeasActual, "SeasType" = s.type)
  return(my_list)
  gc()
}


### Extract STR components
stR_components <- function(object)
{
  len_y <- length(object$input$data)
  len_x <- length(object$output$predictors) + 2
  str_cmp <- matrix(0, len_y, len_x)

  str_cmp[, 1] <- as.vector(object$input$data)
  str_cmp[, ncol(str_cmp)] <- as.vector(object$output$random$data)
  names <- rep("", ncol(str_cmp))
  names[c(1, ncol(str_cmp))] = c("Data", "Random")

  for(i in seq_along(object$output$predictors)) {
    str_cmp[, i+1] <- object$output$predictors[[i]]$data
    names[i+1] <- object$input$predictors[[i]]$name
  }
  colnames(str_cmp) <- names

  str_cmp <- ts(str_cmp)
  if("ts" %in% class(object$input$data))
    tsp(str_cmp) <- tsp(object$input$data)
  return(str_cmp)
}


### Seasonal adjustment based on STR
stR_seasadj <- function(object, include = c("Trend", "Random"))
{
  str_cmp <- stR_components(object)
  nameTrend <- colnames(str_cmp)[2]
  if(is.null(nameTrend) || is.na(nameTrend) || nchar(nameTrend) == 0) {
    warning("Trend component is not specified by name, using the first component as the Trend component.")
    colnames(str_cmp)[2] <- "Trend"
  }
  for(cmpname in include[!(include %in% colnames(str_cmp))]) {
    warning(paste(cmpname, "is not one of the components of the decomposion, skipping..."))
  }
  result <- NULL
  for(i in include[include %in% colnames(str_cmp)]) {
    if(is.null(result)) {
      result <- str_cmp[,i]
    } else {
      result <- result + str_cmp[,i]
    }
  }
  return(result)
}
