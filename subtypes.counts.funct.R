#' subtypes.counts.funct Function
#'
#' @export
#' @examples
#' subtypes.counts.funct()
subtypes.counts.funct <- function(subtype.list, genomes){
  subtype.count <- vector(mode = 'numeric', length = length(subtype.list))

  ##loop through the subtypes
  for(i in 1:length(subtype.list)){
    ##find all genomes containing a subtype
    x <- grep(subtype.list[i], genomes$subtypes)
    subtype.count[i] <- length(x)
  }
  subtypes.dat <- data.frame(subtype = subtype.list, counts = subtype.count)
  return(subtypes.dat)
}
