#' updateCRISPRPackage Function
#'
#' @export
#' @examples
#' updateCRISPRPackage()
updateCRISPRPackage <- function(){
library("devtools")
library(roxygen2)
setwd('~/Desktop/CRISPRSpacerTools/')
document()
}
