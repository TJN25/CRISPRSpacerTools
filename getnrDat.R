#' getnrDat Function
#'
#' @param Dat The CRISPRTarget file that is loaded
#' @param threshold The number of protospacers in a target genome needed to keep the output
#' @export
#' @examples
#' getnrDat()
getnrDat <- function(Dat,threshold){
  bb <- 0
  x <- vector(mode="character", length=length(Dat$Score))
  nrDat <- cbind(Dat,x,x)
  rm(x)
  colnames(nrDat) <- c(colnames(Dat),"Keep",'Target_name')
  nrDat <- transform(nrDat, Target_name = as.character(Target_name))
  levels(nrDat$Keep) <- c('n','y')
  for(g in levels(Dat$Spacer_ID)){
    arrayGenomeDat <- getSubsetBasedOnRows(Dat,1,g)
    if(length(arrayGenomeDat$Spacer_Number) > threshold){
      for(j in levels(arrayGenomeDat$Protospacer_seq_id)){
        targetDat <- getSubsetBasedOnRows(arrayGenomeDat,2,j)
        if(length(targetDat$Spacer_Number) > threshold){
          for(k in levels(targetDat$Array_Number)){
            arrayDat <- getSubsetBasedOnRows(targetDat,7,k)
            if(length(unique(arrayDat$Spacer_Number)) >= threshold){
              bb <- 0
              target_name <- paste(g,j,k, sep='_')
              ww <- nrDat[nrDat[,1]==g,]
              ww <- ww[ww[,2]==j,]
              ww <- getSubsetBasedOnRows(ww,7,k)
              ww <- transform(ww, Spacer_Number = as.factor(Spacer_Number))
              for(ll in levels(ww$Spacer_Number)){
                nrDat <- removeRedundantProtospacers(nrDat,ll,ww,target_name)
              }
              rm(ww)
            }
          }
        }
      }
    }
  }
  
  if(bb==0){
    nrDat <- nrDat[nrDat[,9]=='y',]
    return(nrDat)
  }else{
    print('Nothing found')
  }
}
