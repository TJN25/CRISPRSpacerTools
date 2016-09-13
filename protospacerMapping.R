#' protospacerMapping Function
#'
#' @param nrDat Data frame containing info about all the protospacers that are clusterd
#' @export
#' @examples
#' protospacerMapping()

protospacerMapping <- function(nrDat){
  adjusted.protospacer.positions <- vector(mode="numeric", length=length(nrDat$Score))
  nrDat <- cbind(nrDat, adjusted.protospacer.positions)
  nrDat <- nrDat[order(-nrDat$Spacer_Number),]
  nrDat <- nrDat[order(nrDat$Target_name),]
  for(i in levels(nrDat$Target_name)){
    arrayDat <- getSubsetBasedOnRows(nrDat,10,i)
    pos <- arrayDat[1,3]
    for(j in 1:length(arrayDat$Protospacer_start)){
      row_num <- as.character(rownames(arrayDat[j,]))
      nrDat[row_num,12] <- arrayDat$Protospacer_start[j] - pos
    }
  }
  plottingDat <- nrDat[nrDat[,12]<30000,]
  plottingDat <- plottingDat[plottingDat[,12]>-30000,]
  plottingDat <- plottingDat[plottingDat[,12]!=0,]
  if(length(plottingDat$adjusted.protospacer.positions)!=0){
  plottingDat[,12] <- signif(plottingDat$adjusted.protospacer.positions, 1)
  barplot(prop.table(table(plottingDat$Protospacer_Strand, plottingDat$adjusted.protospacer.positions))*length(plottingDat$adjusted.protospacer.positions), beside=T, legend=T, col=c('blue','red'))
  return(plottingDat)
  }else{
    print('No protospacers within the range')
  }
}



