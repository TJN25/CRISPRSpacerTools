#' setup Function
#'
#' This function allows you to set the working directory and install the needed packages.
#' @keywords setup
#' @export
#' @examples
#' setup()
setup <- function(path_to_wd){
  library(ggplot2)
  library(GenomeGraphs)
  library(genoPlotR)
  library(RColorBrewer)
  library(lattice)
  library(genomation)
  library(parallel)
  setwd(path_to_wd)
}
