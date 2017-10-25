#' plotGenomeLengthDistribution Function
#'
#' @param clust.res data file
#' @export
#' @examples
#' plotGenomeLengthDistribution()
plotGenomeLengthDistribution <- function(Subtype.label = "I-F", dat = targets.dat){
  binwidth.val <- 200
  distance.window <- 10000
  ##store input variable for use later
  input.dat <- dat%>%filter(Subtype == Subtype.label)
  input.subtype <- Subtype.label

  dat <- dat%>%filter(Subtype == Subtype.label)

  dat <- dat%>%filter(spacer_order.num > 1)#%>%filter(protospacer.distance.num != -1)%>%filter(protospacer.distance.num != 1)%>%filter(protospacer.distance.num != -2)%>%filter(protospacer.distance.num != 2)%>%filter(protospacer.distance.num != -3)%>%filter(protospacer.distance.num != 3)

  dat <- dat%>%filter(protospacer.distance.num > -xlim.num)%>%filter(protospacer.distance.num < xlim.num)
  ##get maximum count for the graphs
  targets.dat.n0 <- dat%>%filter(Subtype == Subtype.label)
  aa <- targets.dat.n0%>%filter(strand.plus.direction == "n_3")
  bb <- targets.dat.n0%>%filter(strand.plus.direction == "n_5")
  cc <- targets.dat.n0%>%filter(strand.plus.direction == "t_3")
  dd <- targets.dat.n0%>%filter(strand.plus.direction == "t_5")

  bins.max <- max(c(max(hist(aa$protospacer.distance.num, breaks = (distance.window)/binwidth.val, plot = F)$counts),
                    max(hist(bb$protospacer.distance.num, breaks = (distance.window)/binwidth.val, plot = F)$counts),
                    max(hist(cc$protospacer.distance.num, breaks = (distance.window)/binwidth.val, plot = F)$counts),
                    max(hist(dd$protospacer.distance.num, breaks = (distance.window)/binwidth.val, plot = F)$counts)))

  ##set up data for plotting
  dat <- dat%>%mutate(protospacer.distance.num = ifelse(five.three.prime.dir == "3", protospacer.distance.num + binwidth.val/2,  protospacer.distance.num - binwidth.val/2))
  dat <- dat%>%mutate(protospacer.distance.num = ifelse(target.strand == "t", protospacer.distance.num, protospacer.distance.num*(-1)))

  targets.dat.n0 <- dat%>%filter(Subtype == Subtype.label)


  plot.title1 <- ifelse(distance.window != binwidth.val, paste("RANDOM DATA: Distribution of Subtype ", Subtype.label, " hits (", nrow(targets.dat.n0), " hits.)" , sep = ""),  paste("Quadrant distribution of Subtype ", Subtype.label, " hits (", nrow(targets.dat.n0), " hits)", sep = ""))
  plot.subtitle <- ifelse(distance.window != binwidth.val, paste("Window size = ", distance.window, " nucleotides. \nBinwdith = ", binwidth.val, " nucleotides.", sep = ""),  paste("Window size = ", distance.window, " nucleotides.", sep = ""))

  GenomeEndDensities <- genome.ends.distr(Subtype.label = input.subtype, dat = input.dat)
  GenomeEndDensities <- GenomeEndDensities%>%arrange(distance.breaks.short)

  xlim.num <- calculateXlimits(input.dat)

  p <- ggplot() +
    geom_path(data = GenomeEndDensities, aes(x = distance.breaks.short, y = density.values)) +
    coord_cartesian(xlim = c(-xlim.num, xlim.num))



  return(p)
}
