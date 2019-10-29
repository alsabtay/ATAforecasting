#' Specialized Plot Function of The ATA Method
#'
#' @param x an object of \code{ata}
#' @param fcol line color
#' @param flty line type
#' @param flwd line width
#' @param ...
#'
#' @return a graphic output for the components of the ATA Methods
#' @export
#'
plot.ata <- function(x, fcol=4, flty = 2, flwd = 2, ...)
{
  par.default <- par(no.readonly = TRUE)# save default, for resetting...
  caption <- paste(ifelse(x$model.type=="A"," Additive "," Multiplicative "), x$method, sep="")
  xx <- x$actual
  hpred <- length(x$forecast)
  freq <- frequency(xx)
  strt <- start(xx)
  xxx <- ts(c(x$actual, rep(NA,hpred)), end=tsp(xx)[2] + hpred/freq, frequency=freq)
  xxy <- ts(c(x$fitted, rep(NA,hpred)), end=tsp(xx)[2] + hpred/freq, frequency=freq)
  min_y <- min(x$actual, x$fitted, x$out.sample, x$forecast, x$forecast.lower, na.rm=TRUE)
  max_y <- max(x$actual, x$fitted, x$out.sample, x$forecast, x$forecast.upper, na.rm=TRUE)
  range_y <- abs(max_y - min_y)
  min_last <- floor(min_y - range_y * 0.20)
  max_last <- ceiling(max_y + range_y * 0.20)
  range_last <- abs(max_last - min_last)
  dataset <- cbind(xxx,xxy)
  colnames(dataset, do.NULL = FALSE)
  colnames(dataset) <- c("actual","fitted")
  legend_names <- c("actual","fitted","out-sample","forecast")
  tmp <- seq(from = tsp(x$forecast)[1], by = 1/freq, length = hpred)
  if (x$is.season==FALSE){
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(dataset,plot.type="s", ylim=c(min_last, max_last), col=1:ncol(dataset), xlab=NULL, ylab="fitted", yaxt="n")
    axis(side=2,at=seq(min_last, max_last,trunc(range_last/10)), labels=seq(min_last, max_last,trunc(range_last/10)), las=1, lwd=1)
    polygon(x=c(tmp, rev(tmp)), y=c(x$forecast.lower, rev(x$forecast.upper)), col="lightgray", border=NA)
    lines(x$forecast, lty = flty, lwd = flwd, col = fcol)
    lines(x$out.sample, lty = 1, lwd = flwd, col = fcol+2)
    legend("topleft", legend_names, col=c(1,2,fcol+2,fcol), lty=1, cex=.80, box.lty=0, text.font=2, ncol=2,  bg="transparent")
    mtext(caption, side = 3, line = -1.5, outer = TRUE)
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$trend, ylab="trend")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$level, ylab="level")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$residuals, ylab="residuals")
  }else {
    layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow=TRUE))
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(dataset,plot.type="s", ylim=c(min_last, max_last), col=1:ncol(dataset), xlab=NULL, ylab="fitted", yaxt="n")
    axis(side=2,at=seq(min_last, max_last,trunc(range_last/10)), labels=seq(min_last, max_last,trunc(range_last/10)), las=1, lwd=1)
    polygon(x=c(tmp, rev(tmp)), y=c(x$forecast.lower, rev(x$forecast.upper)), col="lightgray", border=NA)
    lines(x$forecast, lty = flty, lwd = flwd, col = fcol)
    lines(x$out.sample, lty = 1, lwd = flwd, col = fcol+2)
    legend("topleft", legend_names, col=c(1,2,fcol+2,fcol), lty=1, cex=.80, box.lty=0, text.font=2, ncol=2, bg="transparent")
    mtext(paste(caption,"with ",ifelse(x$seasonal.type=="A","Additive","Multiplicative"), " Decomposition by '",ifelse(x$seasonal.model=="decomp","classical",x$seasonal.model),"' Method"), side = 3, line = -1.5, outer = TRUE)
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$seasonal.adjusted,ylab="deseasonalized")
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$level, ylab="level")
    par(mar = c(bottom=1, 4.1, top=2, 1.1))
    plot(x$trend,ylab="trend")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$seasonal,ylab="seasonality")
    par(mar = c(bottom=2, 4.1, top=2, 1.1))
    plot(x$residuals, ylab="residuals")
  }
  par(par.default)
}
