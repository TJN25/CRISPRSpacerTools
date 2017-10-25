#' calculateXlimits Function
#'
#' @param filename The name of the CRISPRTarget file
#' @param  path_to_wd The path to the working directory
#' @param keep_selfmatch Whether to keep the results that were potentially matched to self. Default is FALSE
#' @keywords import data
#' @export
#' @examples
#' calculateXlimits()
calculateXlimits <- function(dat){
  sdNumbers <- 2
  ##get some of the oldest match information
  pps.dat <- dat%>%filter(spacer_order.num == 1)%>%mutate(ppsToGenomeEnds = genome.length.num - target.start.num)
  xlim.num <- max(c(mean(pps.dat$target.start.num) + sdNumbers*sd(pps.dat$target.start.num),
                    mean(pps.dat$ppsToGenomeEnds) + sdNumbers*sd(pps.dat$ppsToGenomeEnds)))
  return(xlim.num)
}
