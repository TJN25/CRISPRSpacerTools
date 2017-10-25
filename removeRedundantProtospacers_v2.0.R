#' removeRedundantProtospacers_v2.0 Function
#'
#' @export
#' @examples
#' removeRedundantProtospacers_v2.0()
removeRedundantProtospacers_v2.0 <- function(nrDat,ll,single_array){
  single_spacer <- single_array[single_array[,8]==ll,]
  single_spacer <- single_spacer[order(single_spacer$Score),]
  single_spacer <- single_spacer[1,]
  row_num <- as.character(rownames(single_spacer))
  nrDat[row_num,9] <- 'y'
  return(nrDat)
}
