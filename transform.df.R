#' transform.df Function
#'
#' @export
#' @examples
#' transform.df()
transform.df <- function(dat){
  ##loop that will go across the columns and look at each header
  for(i in 1:length(dat[1,])){

    ##checks if 'num' is in the header and assigns numeric if true
    if(substr(colnames(dat)[i],(nchar(colnames(dat)[i])-3),nchar(colnames(dat)[i]) )=='.num'){
      dat[,i] <- as.character(dat[,i])
      dat[,i] <- as.numeric(dat[,i])

      ##checks if 'y.n' is in the header and assigns logical if true
    }else if(substr(colnames(dat)[i], nchar(colnames(dat)[i])-3, nchar(colnames(dat)[i]))=='.y.n'){
      dat[,i] <- as.logical(dat[,i])

      ##assigns character to all the remaining columns
    }else{
      dat[,i] <- as.character(dat[,i])
    }
  }
  return(dat)
}
