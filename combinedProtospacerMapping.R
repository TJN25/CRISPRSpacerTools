#' combinedProtospacerMapping Function
#'
#' @param nrDat Data frame containing info about all the protospacers that are clusterd
#' @export
#' @examples
#' combinedProtospacerMapping()
###Need to edit this script. The naive dist needs to be set up to include a simulation of strand.
combinedProtospacerMapping <- function(nrDat, pos_strand, neg_strand){
  adjusted.protospacer.positions <- vector(mode="numeric", length=length(nrDat$Score))
  nrDat <- cbind(nrDat, adjusted.protospacer.positions)
  nrDat <- nrDat[order(-nrDat$Spacer_Number),]
  nrDat <- nrDat[order(nrDat$Target_name),]
  xx <- nrDat$Target_name
  xx <- as.character(xx)
  xx <- unique(xx)
  xx <- as.factor(xx)
  for(i in levels(xx)){
    arrayDat <- getSubsetBasedOnRows(nrDat,10,i)
    pos <- arrayDat[1,3]
    for(j in 1:length(arrayDat$Protospacer_start)){
      row_num <- as.character(rownames(arrayDat[j,]))
      nrDat[row_num,23] <- arrayDat$Protospacer_start[j] - pos
    }
  }
  plottingDat <- nrDat[nrDat[,23]<30000,]
  plottingDat <- plottingDat[plottingDat[,23]>-30000,]
  plottingDat <- plottingDat[plottingDat[,23]!=0,]
  plottingDat[,23] <- plottingDat[,23]/1000
  if(length(plottingDat$adjusted.protospacer.positions)!=0){
  plottingDat[,23] <- signif(plottingDat$adjusted.protospacer.positions, 2)
  pos.dat <- plottingDat[plottingDat[,6]=='+',]
  neg.dat <- plottingDat[plottingDat[,6]=='-',]
  y1 <- runif(1000000, min= 0, max = 20)
  n1 <- runif(1000000, min = -0.5, max=0.5)
  y2 <- runif(1000000, min= 0, max = 20)
  n2 <- runif(1000000, min = -0.5, max=0.5)
  y <- y1+n1 - y2+n2
  rand_start <- vector(mode = 'numeric', length=100000)
  yy <- vector(mode = 'numeric', length=1000000)
  for(i in 1:100000){
    rand_start[i] <- y[i*10]
  }
  for(i in 1:100000){
    a <- rand_start[i]
    for(j in 1:10){
      yy[i+j-1] <- y[i+j-1]-a
    }
  }
  yy <- yy[yy!=0]

  if(pos_strand == 'y'){
  pos.den <- density(pos.dat$adjusted.protospacer.positions)
  pos.max <- max(pos.den$y)
  }else{
    pos.den <- 0
    pos.max <- 0
  }
  comb.den <- density(plottingDat$adjusted.protospacer.positions)
  comb.max <- max(comb.den$y)
  if(neg_strand == 'y'){
  neg.den <- density(neg.dat$adjusted.protospacer.positions)
  neg.max <- max(neg.den$y)
  }else{
    neg.den <- 0
    neg.max <- 0
  }
  null.den <- density(yy)
  null.max <- max(null.den$y)
  max.list <- c(pos.max,neg.max,comb.max,null.max)
  max.val <- max(max.list)
  plot(density(pos.dat$adjusted.protospacer.positions), col = 'red', main=NA,
       type='n',
       xlab='Distance from the intial protospacer (kb)',
       ylab='Protospacer Density',
       xlim=c(-30,30),
       ylim=c(0,max.val*1.1)
  )
  if(pos_strand == 'y'){
  lines(density(pos.dat$adjusted.protospacer.positions), col = 'red')
  }
  if(neg_strand == 'y'){
  lines(density(neg.dat$adjusted.protospacer.positions), col = 'blue')

  lines(density(plottingDat$adjusted.protospacer.positions), col = 'dark grey', lty=1)
  lines(density(yy), col='black')
  legend("bottomleft",
         text.width=c(15,13.5,14,12.2),
         inset = c(0, -0.25),
         x.intersp=1,
         legend=c('Positive Strand', 'Negative Strand', 'Combined', 'Null Distribution'),
         lwd=c(3,3,3,3),
         col=c('red','blue', 'grey','black'),bty='n', cex=0.8,
         xpd=T,
         horiz=TRUE)
  }else{
    ks.res <- suppressWarnings(ks.test(plottingDat$adjusted.protospacer.positions,yy, exact=F))
    ks.p.val <- ks.res$p.value
    hist(plottingDat$adjusted.protospacer.positions,
         breaks=30, xlab='Distance from the intial protospacer (kb)',
         ylab='Protospacer Density',
         xlim=c(-30,30),
         col='Black',
         freq =F,
         main='')
    lines(density(yy), col='Black', lty=2)
    legend("topleft",
          y.intersp = c(0,1.5,1.5,1.5),
          x.intersp = c(0,1.7,1.7,1.7),
          legend=c('Bacterial I-E','Naive Distribution',paste('p value  = ',round(ks.p.val,2),sep=''), 'n = 128'),
           lwd=c(0,3,0,0),
           col=c('black'),
          lty=3,
          bty='n', cex=c(1.2,0.8,0.8,0.8),
           xpd=T
           )
}
  return(plottingDat)
  }else{
    print('No protospacers within the range')
  }
}




