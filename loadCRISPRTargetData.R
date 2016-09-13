#' Function that loads the CRISPRTarget output file
#'
#' @param filename The name of the CRISPRTarget file
#' @param  path_to_wd The path to the working directory
#' @param keep_selfmatch Whether to keep the results that were potentially matched to self. Default is FALSE
#' @keywords import data
#' @export
#' @examples
#' loadCRISPRTargetData()
loadCRISPRTargetData <- function(filename,keep_selfmatch,path_to_wd){
  setwd(path_to_wd)
  if(missing(keep_selfmatch)) {
    keep_selfmatch <- 'F'
  }
  Dat <- read.table(filename, header=T, sep='\t')
  if(keep_selfmatch!='T'){
    Dat <- Dat[Dat[,14]==0,]
    Dat <- cbind(Dat[,1:5],Dat[,8:9],Dat[,20:21])
    Dat <- Dat[,-2,drop=F]
  }else{
    Dat <- cbind(Dat[,1:5],Dat[,8:9],Dat[,20:21])
    Dat <- Dat[,-2,drop=F]
  }
  return(Dat)
}
