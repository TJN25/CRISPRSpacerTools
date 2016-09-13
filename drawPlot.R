#' drawPlot Function
#'
#' @param dfgenes The file produced by setupPlot
#' @export
#' @examples
#' drawPlot()
drawPlot <- function(dfgenes){
  y <- c(barplot(dfgenes$b, dfgenes$a,col=dfgenes$colD,border=dfgenes$colD,ylim=c(0,40), las=2,names.arg=dfgenes$namesDat,cex.names=(0.5)),
         abline(h=25, lty=2),
         abline(h=22, lty=4, col='grey'),
         abline(h=20, lty=3, col='grey'))
  return(y)
}
