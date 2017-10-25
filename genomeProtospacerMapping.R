#' genomeProtospacerMapping Function
#'
#' @param clust.res data file
#' @export
#' @examples
#' genomeProtospacerMapping()

genomeProtospacerMapping <- function(target.dat,genbank_filename,plotname,genome.length){
  if(missing(genbank_filename)){
    stop('Genbank file needed')
  }else if(missing(genbank_filename)){
    stop('Genome Length needed')
  }else{
  w <- vector(mode='character', length=length(target.dat[,1])*2+1)
  x <- vector(mode='numeric', length=length(target.dat[,1])*2+1)
  y <- vector(mode='character', length=length(target.dat[,1])*2+1)
  z <- vector(mode='numeric', length=length(target.dat[,1])*2+1)
  target.dat <- transform(target.dat, gap=as.character(gap))
  target.dat <- transform(target.dat, gap=as.numeric(gap))
  for(i in 1:length(target.dat[,1])){
  x[2*i-1]<- target.dat[i,13]
  x [2*i] <- 100
  y[2*i-1] <- 'white'
  y[2*i] <- 'Red'
  z[2*i-1] <- 4
  z[2*i] <- 4
  w[2*i] <- target.dat[i,8]
  w[2*i-1] <- NA
  }
  x[length(target.dat[,1])*2+1] <- genome.length
  y[length(target.dat[,1])*2+1] <- 'white'
  z[length(target.dat[,1])*2+1] <- 4
  w[length(target.dat[,1])*2+1] <- NA

  gb.dat <- read_dna_seg_from_genbank(genbank_filename, tagsToParse=c("CDS"))
  names <- c(1:length(gb.dat[,1]))
  gb.dat$name <- paste(gb.dat$product,names, sep=' ')
  gene.pos <-gb.dat
  a <- vector(mode='numeric', length=length(gene.pos[,1])*2+1)
  b <- vector(mode='numeric', length=length(gene.pos[,1])*2+1)
  c <- vector(mode='numeric', length=length(gene.pos[,1])*2+1)
  d <- vector(mode='character', length=length(gene.pos[,1])*2+1)
  e <- vector(mode='character', length=length(gene.pos[,1])*2+1)

  aa <- 0
  for(i in 1:length(gene.pos[,1])){
  a[2*i-1] <- abs(gene.pos[i,2]-aa)
  a[2*i] <- gene.pos[i,5]
  aa <- gene.pos[i,3]
  b[2*i-1] <- 'white'
  b[2*i] <- 'black'
  c[2*i] <- 1
  c[2*i-1] <- 1
  d[2*i] <- gene.pos[i,1]
  d[2*i-1] <- NA
  e[2*i] <- ifelse(gene.pos[i,4]==1,'+','-')
  e[2*i-1] <- NA
  }
  a[2*length(gene.pos[,1])+1] <- abs(genome.length - aa)
  b[2*length(gene.pos[,1])+1] <- 'white'
  c[2*length(gene.pos[,1])+1] <- 1
  d[2*length(gene.pos[,1])+1] <- NA
  e[2*length(gene.pos[,1])+1] <- NA
  plotting.dat <- data.frame(x=x, y=y,z=z,spacer.names=w)
  plotting.dat <- transform(plotting.dat, spacer.names=as.character(spacer.names))
  plotting.dat <- transform(plotting.dat, y=as.character(y))
  gene.dat <- data.frame(x=a, y=b,z=c,gene.names=d,strand=e)
  gene.dat <- transform(gene.dat, y=as.character(y))
  gene.dat <- transform(gene.dat, gene.names=as.character(gene.names))
  gene.dat <- transform(gene.dat, strand=as.character(strand))
  x.len <- c()
  y.len <- c()
  x.len.2 <- c()
  y.len.2 <- c()
  col.val <- c()
  names.dat <- c()
  col.dat <- vector(mode='character',length=25)
  blank.names <- vector(mode='character',length=25)
  for(k in 1:25){
    blank.names[k] <- NA
  }
  for(k in 1:25){
    col.dat[k] <- 'black'
  }
  cc <- vector(mode='numeric',length=25)
  cc <- cc+4
  dd <- c(50:26)
  dd <- dd/50
  dd.2 <- c(26:50)
  dd.2 <- dd.2/50
  ee <- c(1:25)
  ee <- ee/50
  ee.2 <- c(25:1)
  ee.2 <- ee.2/50
  for(k in 1:length(gene.dat[,1])){
    if(is.na(gene.dat[k,5])){
      if(k < length(gene.dat[,])){
        bb <- ifelse(gene.dat[k+1,1]>50,gene.dat[k,1],gene.dat[k,1]-49)
      }else{
        bb <- gene.dat[k,1]
      }
      x.len <- c(x.len, bb)
      y.len <- c(y.len, 1)
      x.len.2 <- c(x.len.2,bb)
      y.len.2 <- c(y.len.2, 1)
      col.val <- c(col.val, 'white')
      names.dat <- c(names.dat,NA)
    }else if(gene.dat[k,5]=='+'){
      bb <- gene.dat[k,1]
      bb <- ifelse(bb -50>0,bb-50,1)
      x.len <- c(x.len,bb,cc)
      y.len <- c(y.len, 1,dd)
      x.len.2 <- c(x.len.2,bb,cc)
      y.len.2 <- c(y.len.2, 0,ee)
      col.x <- ifelse(gene.dat[k,1]>100,'black',"grey")
      for(j in 1:25){
        col.dat[j] <- col.x
      }
      col.val <- c(col.val,col.x,col.dat)
      names.dat <- c(names.dat,gene.dat[k,4],blank.names)
    }else{
      bb <- gene.dat[k,1]
      bb <- ifelse(bb -50>0,bb-50,1)
      x.len <- c(x.len,cc,bb)
      y.len <- c(y.len,dd.2,1)
      x.len.2 <- c(x.len.2,cc,bb)
      y.len.2 <- c(y.len.2,ee.2, 0)
      col.x <- ifelse(gene.dat[k,1]>100,'black',"grey")
      for(j in 1:25){
        col.dat[j] <- col.x
      }
      col.val <- c(col.val,col.x,col.dat)
      names.dat <- c(names.dat,blank.names,gene.dat[i,4])
    }
  }
  y.len.3 <- vector(mode='numeric', length = length(y.len))
  y.len.3 <- y.len.3 + 1.5
  positions.protospacers <- vector(mode='numeric', length=length(plotting.dat[,1]))
  pos.a <- 0
  for(i in 1:length(plotting.dat[,1])){
    positions.protospacers[i] <- pos.a + 0.5*plotting.dat[i,1]
    pos.a <- pos.a + plotting.dat[i,1]
  }
  plotting.dat <- cbind(plotting.dat, positions.protospacers)
  if(missing(plotname)){
    par(las=2)
    par(mar = c(7, 4, 2, 2))
    barplot(plotting.dat$z,plotting.dat$x, col=plotting.dat$y, space = 0,border=NA)
    barplot(y.len.3, x.len, col='white',space=0, border='white',add=T)
    barplot(y.len, x.len, col=col.val,space=0, border=col.val, names.arg=names.dat,add=T)
    barplot(y.len.2, x.len.2, col='white',space=0, border='white',add=T)
    text(plotting.dat$positions.protospacers, 2.75,col='black', plotting.dat$spacer.names, cex = 0.7)
    }else{
      pdf(paste(plotname,'genes_and_protospacers.pdf',sep='_'), width=length(gene.pos$name)*12, height=7)
      par(las=2)
      par(mar = c(7, 4, 2, 2))
      barplot(plotting.dat$z,plotting.dat$x, col=plotting.dat$y, space = 0,border=NA)
      barplot(y.len.3, x.len, col='white',space=0, border='white',add=T)
      barplot(y.len, x.len, col=col.val,space=0, border=col.val, names.arg=names.dat,add=T)
      barplot(y.len.2, x.len.2, col='white',space=0, border='white',add=T)
      text(plotting.dat$positions.protospacers, 2.75, plotting.dat$spacer.names, cex = 5)
      dev.off()

  }
}
}
