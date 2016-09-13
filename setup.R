#' Setup Function
#'
#' This function allows you to set the working directory and install the needed packages.
#' @keywords setup
#' @export
#' @examples
#' setup()
setup <- function(){
  library(ggplot2)
  library(GenomeGraphs)
  library(genoPlotR)
  library(RColorBrewer)
  library(lattice)
  library(genomation)
  library(parallel)
}
