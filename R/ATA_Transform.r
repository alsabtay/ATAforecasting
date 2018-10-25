#' @export ATA.Transform

ATA.Transform <- function(X, tMethod, tLambda, tMethod = c("BoxCox", "sqrt", "inverse", "log"), bcMethod = c("guerrero", "loglik"), bcLower = 0, bcUpper = 1){
	tMethod <- match.arg(tMethod)
	if (is.null(tMethod)){
		trfmX <- X
	}else if (tMethod=="BoxCox"){
		if (is.null(tLambda)){
			bcMethod <- match.arg(bcMethod)
			if (bcMethod == "guerrero") {
				tLambda <- BoxCox.lambda(X, method = "guerrero", lower=bcLower, upper=bcUpper)
			} else {
				tLambda <- BoxCox.lambda(X, method = "loglik", lower=bcLower, upper=bcUpper)
			}
			trfmX <- BoxCox(X,tLambda)
		}else { 
			trfmX <- BoxCox(X,tLambda)
		}
		#if (tLambda == 0){ 
		#	trfmX <- log(X) 
		#}else { 
		#	trfmX <- ((X^tLambda)-1)/tLambda
		#}		
	}else if (tMethod=="sqrt"){
		trfmX <- sqrt(X)	
	}else if (tMethod=="inverse"){
		trfmX <- 1/X	
	}else if (tMethod=="log"){
		trfmX <- log(X)	
	}else{
		trfmX <- X
	}
	my_list <- list("trfmX" = trfmX, "tLambda" = tLambda)
	return(my_list)
}