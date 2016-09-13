#' removeRedundantProtospacers Function
#'
#' @export
#' @examples
#' removeRedundantProtospacers()
removeRedundantProtospacers <- function(nrDat,ll,single_array,target_name){
  single_spacer <- single_array[single_array[,8]==ll,]
  single_spacer <- single_spacer[order(single_spacer$Score),]
  single_spacer <- single_spacer[1,]
  row_num <- as.character(rownames(single_spacer))
  nrDat[row_num,9] <- 'y'  
  nrDat[row_num,10] <- target_name
  return(nrDat)
}
