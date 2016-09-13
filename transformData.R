#' transformData Function
#'
#' @export
#' @examples
#' transformData()
transformData <- function(nrDat){
  nrDat <- transform(nrDat, Score = as.character(Score))
  nrDat <- transform(nrDat, Score = as.numeric(Score))
  nrDat <- transform(nrDat, Spacer_Number = as.character(Spacer_Number))
  nrDat <- transform(nrDat, Spacer_Number = as.numeric(Spacer_Number))  
  nrDat <- transform(nrDat, Protospacer_start = as.character(Protospacer_start))
  nrDat <- transform(nrDat, Protospacer_start = as.numeric(Protospacer_start))
  nrDat <- transform(nrDat, Protospacer_stop = as.character(Protospacer_stop))
  nrDat <- transform(nrDat, Protospacer_stop = as.numeric(Protospacer_stop))
  return(nrDat)
}