#' Function that loads the CRISPRTarget output file
#'
#' @param filename The name of the CRISPRTarget file
#' @param  path_to_wd The path to the working directory
#' @param keep_selfmatch Whether to keep the results that were potentially matched to self. Default is FALSE
#' @keywords import data
#' @export
#' @examples
#' loadCRISPRTargetData()
loadCRISPRTargetData <- function(filename,keep_selfmatch=F,threshold,keep_redundant=F,targets_file_name){

  Dat <- read.table(filename, header=T, sep='\t')
  if(keep_selfmatch!='T'){
    Dat <- Dat[Dat[,14]==0,]
    Dat <- cbind(Dat[,1:5],Dat[,8:9],Dat[,20:21],Dat[,10:13])
    Dat <- Dat[,-2,drop=F]
  }else{
    Dat <- cbind(Dat[,1:5],Dat[,8:9],Dat[,20:21])
    Dat <- Dat[,-2,drop=F]
  }
  bb <- 0
  x <- vector(mode="character", length=length(Dat$Score))
  unique.id<- paste(Dat[,1],Dat[,2],Dat[,7], sep = '_')
  Dat <- cbind(Dat, unique.id)
  nrDat <- cbind(Dat[,1:8],x,unique.id,unique.id,Dat[,9:12])
  colnames(nrDat) <- c(colnames(Dat[,1:8]),"Keep",'Target_name',colnames(Dat[,9:12]))
  #nrDat <- transform(nrDat, Target_name = as.character(Target_name))
  levels(nrDat$Keep) <- c('n','y')
  unique.id <- unique(unique.id)
  unique.id <- as.factor(unique.id)
  for(g in levels(unique.id)){
    arrayDat <- getSubsetBasedOnRows(Dat,13,g)
    if(length(unique(arrayDat$Spacer_Number)) >= threshold){
      bb <- 0
      arrayDat <- transform(arrayDat, Spacer_Number = as.factor(Spacer_Number))
      for(ll in levels(arrayDat$Spacer_Number)){
        nrDat <- removeRedundantProtospacers_v2.0(nrDat,ll,arrayDat)
      }
    }
  }
  if(keep_redundant=='F'){
    nrDat <- nrDat[nrDat[,9]=='y',]
  }
return(nrDat)
  if(missing(targets_file_name)) {
    print('targets_file_name not found: not printing target genome IDs')
  }else{
    ids <- unique(nrDat$Protospacer_seq_id)
  write.table(ids, file=targets_file_name, quote=F, sep = '\t')
  }
}
