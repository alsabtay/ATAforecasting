#' @export ATA.Inv.Transform

ATA.Inv.Transform <- function(X, tMethod, tLambda, tbiasadj=FALSE, tfvar=NULL){
	if (is.null(tMethod)){
		trfmX <- X
	}else if (tMethod == "BoxCox") {
		trfmX <- InvBoxCox(X, lambda=tLambda, biasadj=tbiasadj, fvar=tfvar)
		#if (tLambda == 0){ 
		#	trfmX <- exp(X) 
		#}else { 
		#	trfmX <- (tLambda*X + 1)^(1/tLambda) 
		#}		
	}else if (tMethod == "sqrt"){
		trfmX <- X^2	
	}else if (tMethod == "inverse"){
		trfmX <- 1/X	
	}else if (tMethod == "log"){
		trfmX <- exp(X)	
	}else {
		trfmX <- X
	}
	return(trfmX)
}