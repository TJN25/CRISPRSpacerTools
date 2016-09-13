#' genomeLength Function
#'
#' @param genome.length.filename The name of the file containing the information about genome length
#' @param nrDat The data.frame with the results to analyse.
#' @export
#' @examples
#' genomeLength()

genomeLength <- function(genome.length.filename, nrDat){
nrDat <- transform(nrDat, Target_name = as.factor(Target_name))
genome.length <- read.table(genome.length.filename, header=F, sep='\t', stringsAsFactors = F)
Genome_length <- vector(mode="numeric", length=length(nrDat$Score))
nrDat <- cbind(nrDat, Genome_length)
nrDat <- transform(nrDat, Genome_length = as.numeric(Genome_length))

for(j in 1:length(nrDat$Spacer_ID)){
  for(i in 1:length(genome.length[,2])){
    x <- genome.length[i,2]
    if(x==nrDat[j,2]){
      xx <- genome.length[i,1]
      nrDat[j,11] <- genome.length[i,1]
    }
  }
}
return(nrDat)
}