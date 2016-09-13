#' calcGaps Function
#'
#' @param nrDat Data frame containing info about all the protospacers that are clusterd
#' @export
#' @examples
#' calcGaps()

calcGaps <- function(nrDat){
gap <- vector(mode="numeric", length=length(nrDat$Score))
nrDat <- cbind(nrDat, gap)
nrDat <- nrDat[order(nrDat$Protospacer_start),]
nrDat <- nrDat[order(nrDat$Target_name),]

for(i in levels(nrDat$Target_name)){
  arrayDat <- getSubsetBasedOnRows(nrDat,10,i)
  llm <- length(arrayDat$Score)-1
  next_row <- as.character(rownames(arrayDat[1,]))
  for(l in 1:llm){
    row_num <- as.character(rownames(arrayDat[l,]))
    next_row <- as.character(rownames(arrayDat[l+1,]))
    nrDat[row_num,12] <- abs(arrayDat$Protospacer_start[l+1] - arrayDat$Protospacer_stop[l])
  }
  l <- length(arrayDat$Score)
  row_num <- as.character(rownames(arrayDat[l,]))
  nrDat[row_num,12] <- abs(arrayDat$Protospacer_start[l] - arrayDat$Protospacer_stop[l-1])
}
nrDat <- transform(nrDat, gap = as.character(gap))
nrDat <- transform(nrDat, gap = as.numeric(gap))
return(nrDat)
}
