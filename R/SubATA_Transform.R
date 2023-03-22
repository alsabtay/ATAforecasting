SubATA.Transform <- function(tX
                          , tMethod = c("Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog")
                          , tType = c("Vanilla", "Back")
						              , tLambda = NULL
                          , tShift = 0)
{
out_transform <- list("tX" = tX, "tType" = tType, "tLambda" = tLambda, "tShift" = tShift)

switch(tType,
	Vanilla = {
		switch(tMethod,
			Box_Cox = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- log(tX + ntShift)
						}else {
							out_transform$tX <- (((tX + ntShift)^tLambda) - 1) / tLambda
						}
			},
			Modulus = {
						mdls <- abs(tX) + 1
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- sign(tX) * log(mdls)
						}else {
							out_transform$tX <- sign(tX) * ((mdls^tLambda) - 1) / tLambda
						}
			},
			BickelDoksum = {
						if (tLambda > 0.00000000001) {
							out_transform$tX <- ((abs(tX)^tLambda) * sign(tX) - 1)/ tLambda
						}else {
							stop("The lambda parameter must be positive for the Bickel-Doksum transformation.")
						}
			},
			Manly = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- tX
						}else {
							out_transform$tX <- (exp(tX * tLambda) - 1) / tLambda
						}
			},
			Dual = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- log(tX + ntShift)
						}else {
							out_transform$tX <- (((tX + ntShift)^tLambda) - ((tX + ntShift)^(-tLambda))) / (2 * tLambda)
						}
			},
			YeoJohnson = {
						lentX <- length(tX)
						tZ <- rep(NA, lentX)
						negtX <- which(tX < 0)
						postX <- which(tX >= 0)
						if (abs(tLambda) <= 0.00000000001) {
							tZ[postX] <- log(tX[postX] + 1)
						}else {
							tZ[postX] <- (((tX[postX] + 1)^tLambda) - 1) / tLambda
						}
						if (abs(tLambda - 2) <= 0.00000000001) {
							tZ[negtX] <- -log(1 - tX[negtX])
						}else {
							tZ[negtX] <- (((1 - tX[negtX])^(2 - tLambda)) - 1)/(tLambda - 2)
						}
						out_transform$tX <- tZ
			},
			NegLog = {
						mdls <- abs(tX) + 1
						out_transform$tX <- sign(tX) * log(mdls)
						out_transform$tLambda <- NULL

			},
			GLog = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						out_transform$tX <- log((tX + ntShift) + sqrt(((tX + ntShift)^2) + 1))
						out_transform$tLambda <- NULL
			},
			GPower = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- log((tX + ntShift) + sqrt(((tX + ntShift)^2) + 1))
						}else {
							out_transform$tX <- ((((tX + ntShift) + (sqrt(tX + ntShift)^2 + 1))^tLambda) - 1) / tLambda
						}
			},
			Sqrt = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						out_transform$tX <- sqrt(tX + ntShift)
						out_transform$tLambda <- NULL
			},
			Log = {
						out_transform$tShift <- ntShift <- calc_shift(min(tX),tShift)
						out_transform$tLambda <- tLambda <- 0
						out_transform$tX <- log(tX + ntShift)
			},
			Reciprocal = {
						out_transform$tLambda <- NULL
						out_transform$tX <- 1/tX
			}
		)
	},
	Back = {
		switch(tMethod,
			Box_Cox = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- exp(tX) - tShift
						}else {
							out_transform$tX <- (tLambda * tX + 1)^(1 / tLambda) - tShift
						}
			},
			Modulus = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- sign(tX) * (exp(abs(tX)) - 1)
						}else {
							out_transform$tX <- sign(tX) * ((abs(tX)*tLambda + 1)^(1/tLambda) - 1)
						}
			},
			BickelDoksum = {
						negtX <- which(tX < 0)
						postX <- which(tX >= 0)
						lentX <- length(tX)
						tZ <- rep(NA, lentX)
						tZ[postX] <- (tLambda * tX[postX] + 1)^(1 / tLambda)
						tZ[negtX] <- (-1) * ((-1) * (tLambda * tX[negtX] + 1))^(1 / tLambda)
						out_transform$tX <- tZ
			},
			Manly = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- tX
						}else {
							out_transform$tX <- log(tLambda * tX + 1) / tLambda
						}
			},
			Dual = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- exp(tX) + tShift
						}else {
							out_transform$tX <- ((tLambda * tX + sqrt(1 + tLambda^2 * tX^2))^(1/tLambda)) - tShift
						}
			},
			YeoJohnson = {
						lentX <- length(tX)
						tZ <- rep(NA, lentX)
						negtX <- which(tX < 0)
						postX <- which(tX >= 0)
						if (abs(tLambda) <= 0.00000000001) {
							tZ[postX] <- exp(tX[postX]) + 1
						}else {
							tZ[postX] <- ((tX[postX] + 1 * tLambda + 1)^(1 / tLambda)) - 1
						}
						if (abs(tLambda - 2) <= 0.00000000001) {
							tZ[negtX] <- (-1) * (exp(-tX[negtX]) - 1)
						}else {
							tZ[negtX] <- (-1) * ((tX[negtX] * (tLambda - 2) + 1)^(1/(2 - tLambda)) - 1)
						}
						out_transform$tX <- tZ
			},
			NegLog = out_transform$tX <- sign(tX) * (exp(abs(tX)) - 1),
			GLog = out_transform$tX <- ((-(1 - exp(tX * 2))) / (2 * exp(tX))) - tShift,
			GPower = {
						if (abs(tLambda) <= 0.00000000001) {
							out_transform$tX <- (-(1 - exp(tX * 2))) / (2 * exp(tX))
						}else {
						    gpX <- (tX * tLambda + 1)^(1 / tLambda)
							out_transform$tX <- (-(1 - gpX^2)) / (2 * gpX)
						}
			},
			Sqrt = out_transform$tX <- tX^2 + tShift,
			Log = out_transform$tX <- exp(tX) + tShift,
			Reciprocal = out_transform$tX <- 1 / tX
		)
	}
)
return(out_transform)
}


calc_shift <- function(mintX, tshft) {
	if (mintX <= 0) {
		if (tshft >=0 & tshft < abs(mintX)) {
			newtShift <- abs(mintX) + 1
			warning("The data has negative values. Besides, the shift parameter is smaller than minimum value of the data. ATAforecasting changed the shift parameter to absolute minimum value of the data.")
		}else if (tshft >=0 & tshft >= abs(mintX)) {
			newtShift <- tshft + 1
		}else {
			newtShift <- abs(mintX) + 1
			warning("The data has negative values. Besides, the shift parameter must be positive value. ATAforecasting changed the shift parameter to 0.")
		}
	}else {
		if (tshft < 0 & abs(tshft) >= mintX) {
			newtShift <- 0
			warning("The shift parameter is negative and bigger than minimum value of the data. ATAforecasting changed the shift parameter to 0.")
		}else {
			newtShift <- tshft
		}
	}
	return(newtShift)
}
