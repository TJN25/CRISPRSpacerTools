#' identifyClusters Function
#'
#' @export
#' @examples
#' identifyClusters()
identifyClusters <- function(nrDat,threshold){
  protospacer_cluster <- vector(mode="numeric", length=length(nrDat$Score))
  Distance_to_next_spacer <- vector(mode="numeric", length=length(nrDat$Score))
  xx <- cbind(nrDat,Distance_to_next_spacer,protospacer_cluster)
  for(g in levels(nrDat$Spacer_ID)){
    arrayGenomeDat <- getSubsetBasedOnRows(nrDat,1,g)
    for(j in levels(arrayGenomeDat$Protospacer_seq_id)){
      targetDat <- getSubsetBasedOnRows(arrayGenomeDat,2,j)
      for(k in levels(targetDat$Array_Number)){
        cl_num <- 1
        arrayDat <- getSubsetBasedOnRows(targetDat,7,k)
        llm <- length(arrayDat$Score)-1
        next_row <- as.character(rownames(arrayDat[1,]))
        xx[next_row,11] <- cl_num
        for(l in 1:llm){
          row_num <- as.character(rownames(arrayDat[l,]))
          next_row <- as.character(rownames(arrayDat[l+1,]))
          lp <- l + 1
          diff <- abs(arrayDat$Protospacer_start[lp] - arrayDat$Protospacer_stop[l])
          if(diff < threshold){
            xx[next_row,11] <- cl_num
          }else{
            cl_num <- cl_num + 1
            xx[next_row,11] <- cl_num
          }
          xx[row_num,10] <- arrayDat$Protospacer_start[lp] - arrayDat$Protospacer_stop[l]
        }
        l <- length(arrayDat$Score)
        row_num <- as.character(rownames(arrayDat[l,]))
        xx[row_num,10] <- NA

        #print(xx[row_num,])

      }
    }
  }
  return(xx)
}
