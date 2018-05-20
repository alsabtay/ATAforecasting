#' @export ATA.Transform

ATA.Transform <- function(X, tMethod, tLambda){
	if (is.null(tMethod)){
		trfmX <- X
	}else if (tMethod=="BoxCox"){
		if (is.null(tLambda)){
			tLambda <- BoxCox.lambda(X)
			trfmX <- BoxCox(X,tLambda)
		}else { 
			if (tLambda=="opt"){
				tLambda <- BoxCox.lambda(X)
			}
			trfmX <- BoxCox(X,tLambda)
		}
		#if (tLambda == 0){ 
		#	trfmX <- log(X) 
		#}else { 
		#	trfmX <- X^(1/tLambda)
		#}		
	}else if (tMethod=="sqrt"){
		trfmX <- sqrt(X)	
	}else if (tMethod=="inverse"){
		trfmX <- 1/X	
	}else if (tMethod=="log10"){
		trfmX <- log10(X)	
	}else if (tMethod=="log"){
		trfmX <- log(X)	
	}else{
		trfmX <- X
	}
	my_list <- list("trfmX" = trfmX, "tLambda" = tLambda)
	return(my_list)
}