#' getSubsetBasedOnRows Function
#'
#' @export
#' @examples
#' getSubsetBasedOnRows()
getSubsetBasedOnRows <- function(Dat,col_num,g){
  xx <- Dat[Dat[,col_num]==g,]
  xx <- as.matrix(xx)
  xx <- as.data.frame(xx)
  xx <- transform(xx, Score = as.character(Score))
  xx <- transform(xx, Score = as.numeric(Score))
  xx <- transform(xx, Spacer_Number = as.character(Spacer_Number))
  xx <- transform(xx, Spacer_Number = as.numeric(Spacer_Number))
  xx <- transform(xx, Protospacer_start = as.character(Protospacer_start))
  xx <- transform(xx, Protospacer_start = as.numeric(Protospacer_start))
  xx <- transform(xx, Protospacer_stop = as.character(Protospacer_stop))
  xx <- transform(xx, Protospacer_stop = as.numeric(Protospacer_stop))
  return(xx)
  rm(xx)
}
