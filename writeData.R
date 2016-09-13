#' writeData Function
#'
#' @param Dat The data that will be plotted and writen to a file. 
#' @param folder The directory to write the results to.
#' @export
#' @examples
#' writeData()
writeData <- function(Dat,folder){
  if(length(Dat[,1]) == 0){
    print('No significant results to write to file')
  }else{
    for(g in levels(Dat$Spacer_ID)){
      arrayGenomeDat <- getSubsetBasedOnRows(Dat,1,g)
      for(j in levels(arrayGenomeDat$Protospacer_seq_id)){
        targetDat <- getSubsetBasedOnRows(arrayGenomeDat,2,j)
        for(k in levels(targetDat$Array_Number)){
          arrayDat <- getSubsetBasedOnRows(targetDat,7,k)
          dfgenes <- setupPlot(arrayDat)
          exten <- '.pdf'
          pdf(paste(folder,'protospacer_positions_in_the_target_',j,'_for_array_',k,'_from_',g,exten, sep=''), width=300, height=10)
          y <- c(barplot(dfgenes$b, dfgenes$a,col=dfgenes$colD,border=dfgenes$colD,ylim=c(0,40), las=2,names.arg=dfgenes$namesDat,cex.names=(0.5)),
                 abline(h=25, lty=2),
                 abline(h=22, lty=4, col='grey'),
                 abline(h=20, lty=3, col='grey')
          )
          print(y)
          dev.off()
          write.table(arrayDat, paste(folder,'protospacer_positions_in_the_target_',j,'_for_array_',k,'_from_',g,'.txt', sep=''), 
                      quote=F,col.names=T,row.names=F, sep='\t')
          
        }
      }
    }
  }
}
