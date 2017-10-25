#' genome.ends.distr Function
#'
#' @export
#' @examples
#' genome.ends.distr()
genome.ends.distr <- function(Subtype.label = "I-F", dat = targets.dat){

  sdNumbers <- 2

  dat <- dat%>%mutate(ppsToGenomeEnds = genome.length.num - target.start.num)
  pps.dat <- dat%>%filter(spacer_order.num == 1)
  pps.dat <- pps.dat%>%filter(Subtype == Subtype.label)
  xlim.num <- max(c(mean(pps.dat$target.start.num) + sdNumbers*sd(pps.dat$target.start.num),
                    mean(pps.dat$ppsToGenomeEnds) + sdNumbers*sd(pps.dat$ppsToGenomeEnds)))
  binwidth.val <- mean(pps.dat$genome.length.num)/10

  D <- dat%>%filter(spacer_order.num > 1)%>%select(target.start.num, ppsToGenomeEnds)%>%mutate(data.type = "genome lengths")


  GenomeStart <- D %>% group_by(data.type) %>%
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(target.start.num, from = min(.$target.start.num), to = max(.$target.start.num), n = xlim.num/binwidth.val*2))) %>%
    # make a new data.frame from two density objects
    do(data.frame(distance.breaks.short = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  density.values = .$d[[1]]$y))# %>%    # and subtract the y values
  GenomeStart <- GenomeStart%>%mutate(distance.breaks.short = -distance.breaks.short)

  GenomeEnd <- D %>% group_by(data.type) %>%
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(ppsToGenomeEnds, from = min(.$ppsToGenomeEnds), to = max(.$ppsToGenomeEnds), n = xlim.num/binwidth.val*2))) %>%
    # make a new data.frame from two density objects
    do(data.frame(distance.breaks.short = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  density.values = .$d[[1]]$y))# %>%    # and subtract the y values

  den <- rbind(GenomeStart, GenomeEnd)

  return(den)

}
