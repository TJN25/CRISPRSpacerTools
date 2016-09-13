#' setupPlot Function
#'
#' @param nrDat The dataframe that contains the target genome and protospacers to draw.
#' @export
#' @examples
#' setupPlot()
setupPlot <- function(nrDat){
  nrDat <- transformData(nrDat)
  a <- c()
  colD <- c()
  namesDat <- c()
  z <- 0
  b <- c()
  lenArray <- length(nrDat$Protospacer_start)
  arrayNames <- nrDat$Protospacer..5.to.3.
  arrayNames <- as.character(arrayNames)
  ##build a list of the lengths of genes and spaces and the corresponding colouring for graphing
  for(j in 1:lenArray){
    aa <- nrDat$Protospacer_start[j]
    a <- c(a, abs(aa-z))
    a <- c(a, abs(nrDat$Protospacer_stop[j] - nrDat$Protospacer_start[j]))
    colD <- c(colD, 'black')
    namesDat <- c(namesDat, NA)
    spacerName <- paste(nrDat$Array_Number[j],nrDat$Spacer_Number[j],sep = '-')
    namesDat <- c(namesDat, spacerName)
    b <- c(b, 0, nrDat$Score[j])
    z <- nrDat$Protospacer_stop[j]
    colD <- c(colD, ifelse(nrDat$Protospacer_Strand[j]=='+','red','blue'))

  }
  nrDat <- transform(nrDat, Genome_length = as.character(Genome_length))
  nrDat <- transform(nrDat, Genome_length = as.numeric(Genome_length))
  a <- c(a, nrDat$Genome_length[1] - nrDat$Protospacer_stop[length(nrDat[,1])])
  colD <- c(colD, 'black')
  namesDat <- c(namesDat, NA)
  b <- c(b,0)
  ##Set up the barplot heights for the spaces and proteins for mapping
  ##Make the data frames for the gene plot and the score plot
  dfgenes <- data.frame(a,b,colD,namesDat)
  dfgenes <- transform(dfgenes, colD = as.character(colD))
  return(dfgenes)
}
