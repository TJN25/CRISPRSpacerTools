#' plotProtospacerDistribution Function
#'
#' @param Dat The CRISPRTarget file that is loaded
#' @param pval The cutoff pvalue to use for clustering. Default is 0.01
#' @export
#' @examples
#' plotProtospacerDistribution()

plotProtospacerDistribution <- function(Subtype.label = "I-F", dat = targets.dat){
  binwidth.val <- 200
  distance.window <- 10000
  ##store input variable for use later
  input.dat <- dat
  input.subtype <- Subtype.label

  dat <- dat%>%filter(Subtype == Subtype.label)

  dat <- dat%>%filter(spacer_order.num > 1)#%>%filter(protospacer.distance.num != -1)%>%filter(protospacer.distance.num != 1)%>%filter(protospacer.distance.num != -2)%>%filter(protospacer.distance.num != 2)%>%filter(protospacer.distance.num != -3)%>%filter(protospacer.distance.num != 3)

  dat <- dat%>%filter(protospacer.distance.num > -distance.window)%>%filter(protospacer.distance.num < distance.window)
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


  plot.title1 <- ifelse(distance.window != binwidth.val, paste("Distribution of Subtype ", Subtype.label, " hits (", nrow(targets.dat.n0), " hits.)" , sep = ""),  paste("Quadrant distribution of Subtype ", Subtype.label, " hits (", nrow(targets.dat.n0), " hits)", sep = ""))
  plot.subtitle <- ifelse(distance.window != binwidth.val, paste("Window size = ", distance.window, " nucleotides. \nBinwdith = ", binwidth.val, " nucleotides.", sep = ""),  paste("Window size = ", distance.window, " nucleotides.", sep = ""))

  p <- ggplot() +
    geom_histogram(data=subset(targets.dat.n0, strand.plus.direction=="t_5"), binwidth = binwidth.val, aes(protospacer.distance.num, fill = "Target 5' direction", y= ..count..)) +
    geom_histogram(data=subset(targets.dat.n0, strand.plus.direction=="t_3"), binwidth = binwidth.val, aes(protospacer.distance.num, fill = "Target 3' direction", y= ..count..)) +
    geom_histogram(data=subset(targets.dat.n0, strand.plus.direction=="n_5"), binwidth = binwidth.val, aes(protospacer.distance.num, fill = "Non-target 5' direction", y= -..count..)) +
    geom_histogram(data=subset(targets.dat.n0, strand.plus.direction=="n_3"), binwidth = binwidth.val, aes(protospacer.distance.num, fill = "Non-target 3' direction", y= -..count..)) +
    #facet_wrap(~protospacer.distance.num)
    scale_fill_hue("Group") +
    ggtitle(label = plot.title1, subtitle = plot.subtitle) +
    labs(x="Distance from oldest protospacer (nucleotides)",y="Number of hits") +
    coord_cartesian(ylim = c(-bins.max - 3, bins.max + 3)) +
    theme_bw() +
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          plot.title=element_text(size=12, face="bold", color="black"))
  return(p)
}

