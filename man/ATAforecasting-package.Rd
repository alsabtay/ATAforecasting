\name{ATAforecasting-package}
\alias{ATAforecasting-package}
\alias{ATAforecasting}
\docType{package}
\title{
  Forecasting Time Series by ATA Method
}
\description{
  Based on the modified simple exponential smoothing as described in Yapar, G. (2016).
}
\details{
  Based on the modified simple exponential smoothing as described in Yapar, G. (2016). 
  ATA method is a new univariate time series forecasting method which provides innovative solutions to 
  issues faced during the initialization and optimization stages of existing methods. 
  ATA forecasting performance is superior to existing methods both in terms of easy implementation and accurate forecasting. 
  It can be applied to non-seasonal or deseasonalized time series, where the deseasonalization can be performed via any preferred decomposition method.
  This methodology performed extremely well on the M3 and M4-competition data.
}
\author{
Ali Sabri Taylan <alisabritaylan@gmail.com> and Hanife Taylan Selamlar <hanife.taylan@gmail.com>

Maintainer: Ali Sabri Taylan <alisabritaylan@gmail.com> and Hanife Taylan Selamlar <hanife.taylan@gmail.com>
}
\references{
	Yapar, G., (2016)
	"Modified simple exponential smoothing"
	\emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi:10.15672/HJMS.201614320580

	Yapar, G., Capar, S., Selamlar, H. T., Yavuz, I., (2016) 
	"Modified holt's linear trend method"
	\emph{Hacettepe University Journal of Mathematics and Statistics} Early Access. Doi: 10.15672/HJMS.2017.493
}
\keyword{accuracy}
\keyword{ata}
\keyword{forecast}
\keyword{msts}
\keyword{ts}
\seealso{
  \code{\link{forecast}}, \code{\link{stlplus}}, \code{stR}, \code{\link[stats]{stl}}, \code{\link[stats]{decompose}},
	\code{\link{tbats}}, \code{\link{seasadj}
}
\examples{
  \dontrun{
	fit <- ATA(M3[[1899]]$x,M3[[1899]]$xx)
	plot(ATA.Forecast(fit,h=36)) 
  }
}
