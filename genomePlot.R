#' genomePlot Function
#'
#' @param clust.res data file
#' @export
#' @examples
#' genomePlot()

genomePlot <- function(target.dat,genbank_filename,plotname,genome.length,plot_scores,score.val.input){
  if(missing(genbank_filename)){
    stop('Genbank file needed')
  }else if(missing(genbank_filename)){
    stop('Genome Length needed')
  }else if(missing(plot_scores)){
    plot_scores='F'
  }else{
    w <- vector(mode='character', length=length(target.dat[,1])*2+1)
    x <- vector(mode='numeric', length=length(target.dat[,1])*2+1)
    y <- vector(mode='character', length=length(target.dat[,1])*2+1)
    z <- vector(mode='numeric', length=length(target.dat[,1])*2+1)
    target.dat <- transform(target.dat, gap=as.character(gap))
    target.dat <- transform(target.dat, gap=as.numeric(gap))
    for(i in 1:length(target.dat[,1])){
      x[2*i-1]<- target.dat[i,13] - (100-abs(target.dat[i,4]-target.dat[i,3]))
      x [2*i] <- 100
      y[2*i-1] <- 'white'
      y[2*i] <- 'Red'
      z[2*i-1] <- 4
      z[2*i] <- 4
      w[2*i] <- target.dat[i,8]
      w[2*i-1] <- NA
    }
    x[length(target.dat[,1])*2+1] <- genome.length - sum(x[1:(length(target.dat[,1])*2)])
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
    score.val.dat <- vector(mode='numeric', length=length(gene.pos[,1])*2+1)
    aa <- 0
    for(i in 1:length(gene.pos[,1])){
      if(gene.pos[i,2]< gene.pos[i,3]){
      a[2*i-1] <- abs(gene.pos[i,2] - aa)
      a[2*i] <- gene.pos[i,3] - gene.pos[i,2]
      aa <- gene.pos[i,3]
      }else{
        a[2*i-1] <- abs(gene.pos[i,3] - aa)
        a[2*i] <- gene.pos[i,2] - gene.pos[i,3]
        aa <- gene.pos[i,2]

      }
      b[2*i-1] <- 'white'
      b[2*i] <- 'black'
      c[2*i] <- 1
      c[2*i-1] <- 1
      d[2*i] <- gene.pos[i,1]
      d[2*i-1] <- NA
      e[2*i] <- ifelse(gene.pos[i,4]==1,'+','-')
      e[2*i-1] <- NA
      score.val.dat[2*i-1] <- 1
      score.val.dat[2*i] <- score.val.input[i]
    }
    a[2*length(gene.pos[,1])+1] <- genome.length - sum(a[1:(2*length(gene.pos[,1]))])
    b[2*length(gene.pos[,1])+1] <- 'white'
    c[2*length(gene.pos[,1])+1] <- 1
    d[2*length(gene.pos[,1])+1] <- NA
    e[2*length(gene.pos[,1])+1] <- NA
    score.val.dat[2*length(gene.pos[,1])+1] <- 1
    plotting.dat <- data.frame(x=x, y=y,z=z,spacer.names=w)
    plotting.dat <- transform(plotting.dat, spacer.names=as.character(spacer.names))
    plotting.dat <- transform(plotting.dat, y=as.character(y))
    gene.dat <- data.frame(x=a, y=b,z=c,gene.names=d,strand=e, score.val=score.val.dat)
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
    cc.2 <- cc + 1
    cc <- cc+2
    dd <- c(50:26)
    dd <- dd/50
    dd.2 <- c(26:50)
    dd.2 <- dd.2/50
    ee <- c(1:25)
    ee <- ee/50
    ee.2 <- c(25:1)
    ee.2 <- ee.2/50
    score.val <- c()
    for(k in 1:length(gene.dat[,1])){
      if(is.na(gene.dat[k,5])){
        bb <- gene.dat[k,1]
        x.len <- c(x.len, bb)
        y.len <- c(y.len, 1)
        x.len.2 <- c(x.len.2,bb)
        y.len.2 <- c(y.len.2, 1)
        col.val <- c(col.val, 'white')
        names.dat <- c(names.dat,NA)
        score.val <- c(score.val, 0)
      }else if(gene.dat[k,5]=='+'){
        bb <- gene.dat[k,1]
        bb.2 <- ifelse(bb -50>0,bb-50,1)
        x.len <- c(x.len,bb.2,cc)
        y.len <- c(y.len, 1,dd)
        x.len.2 <- c(x.len.2,bb.2,cc)
        y.len.2 <- c(y.len.2, 0,ee)
        col.x <- ifelse(gene.dat[k,1]>100,'black',"grey")
        score.val <- c(score.val, score.val.dat[k],cc.2)
        for(j in 1:25){
          col.dat[j] <- col.x
        }
        col.val <- c(col.val,col.x,col.dat)
        names.dat <- c(names.dat,gene.dat[k,4],blank.names)
      }else{
        bb <- gene.dat[k,1]
        bb.2 <- ifelse(bb -50>0,bb-50,1)
        x.len <- c(x.len,cc,bb.2)
        y.len <- c(y.len,dd.2,1)
        x.len.2 <- c(x.len.2,cc,bb.2)
        y.len.2 <- c(y.len.2,ee.2, 0)
        col.x <- ifelse(gene.dat[k,1]>100,'black',"grey")
        score.val <- c(score.val, cc.2,score.val.dat[k])

        for(j in 1:25){
          col.dat[j] <- col.x
        }
        col.val <- c(col.val,col.x,col.dat)
        names.dat <- c(names.dat,blank.names,gene.dat[k,4])
      }
    }
    y.len.3 <- vector(mode='numeric', length = length(y.len))
    y.len.3 <- y.len.3 + 1.5
    genes.dat <- data.frame(x.dists=x.len,y.heights=y.len,x.dists.arrows=x.len.2,y.heights.arrows=y.len.2,
                            gene.colours=col.val, gene.names=names.dat, y.white.space=y.len.3,score=score.val)
    genes.dat <- transform(genes.dat, gene.colours=as.character(gene.colours))
    positions.protospacers <- vector(mode='numeric', length=length(plotting.dat[,1]))
    positions.genes <- vector(mode='numeric', length=length(genes.dat[,1]))

    pos.a <- 0
    for(i in 1:length(plotting.dat[,1])){
      positions.protospacers[i] <- pos.a + 0.5*plotting.dat[i,1]
      pos.a <- pos.a + plotting.dat[i,1]
    }
    pos.a <- 0
    for(i in 1:length(genes.dat[,1])){
      positions.genes[i] <- pos.a + 0.5*genes.dat[i,1]
      pos.a <- pos.a + genes.dat[i,1]
    }
    plotting.dat <- cbind(plotting.dat, positions.protospacers)
    genes.dat <- cbind(genes.dat, positions.genes)
    if(plot_scores=='F'){
      if(missing(plotname)){
        barplot(plotting.dat$z,plotting.dat$x, col=plotting.dat$y, space = 0,border=NA)
        barplot(genes.dat$y.white.space, genes.dat$x.dists, col='white',space=0, border='white',add=T)
        barplot(genes.dat$y.heights, genes.dat$x.dists, col=genes.dat$gene.colours,space=0,
                border=genes.dat$gene.colours, add=T)
        barplot(genes.dat$y.heights.arrows, genes.dat$x.dists.arrows, col='white',space=0, border='white',add=T)
        text(plotting.dat$positions.protospacers, 2.75,col='black', plotting.dat$spacer.names, cex = 0.7)
      }else{
        pdf(paste(plotname,'genes_and_protospacers.pdf',sep='_'), width=length(gene.pos$name)*12, height=7)
        barplot(plotting.dat$z,plotting.dat$x, col=plotting.dat$y, space = 0,border=NA)
        barplot(genes.dat$y.white.space, genes.dat$x.dists, col='white',space=0, border='white',add=T)
        barplot(genes.dat$y.heights, genes.dat$x.dists, col=genes.dat$gene.colours,space=0,
                border=genes.dat$gene.colours, add=T)
        barplot(genes.dat$y.heights.arrows, genes.dat$x.dists.arrows, col='white',space=0, border='white',add=T)
        text(plotting.dat$positions.protospacers, 2.75, plotting.dat$spacer.names, cex = 5)
        text(genes.dat$positions.genes, 0.5, genes.dat$gene.names, cex = 1, col = 'white')
        dev.off()

      }
    }else
      if(missing(plotname)){
        barplot(genes.dat$score,genes.dat$x.dists, col=genes.dat$gene.colours, space = 0,border=NA)
        barplot(genes.dat$y.white.space, genes.dat$x.dists, col='white',space=0, border='white',add=T)
        barplot(genes.dat$y.heights, genes.dat$x.dists, col=genes.dat$gene.colours,space=0,
                border=genes.dat$gene.colours, names.arg=genes.dat$gene.names,add=T)
        barplot(genes.dat$y.heights.arrows, genes.dat$x.dists.arrows, col='white',space=0, border='white',add=T)

        #text(plotting.dat$positions.protospacers, 2.75,col='black', plotting.dat$spacer.names, cex = 0.7)
      }else{
        pdf(paste(plotname,'genes_and_protospacers.pdf',sep='_'), width=length(gene.pos$name)*12, height=7)
        barplot(genes.dat$score,genes.dat$x.dists, col=genes.dat$gene.colours, space = 0,border=NA)
        barplot(y.len.3, x.len, col='white',space=0, border='white',add=T)
        barplot(y.len, x.len, col=col.val,space=0, border=col.val, add=T)
        barplot(y.len.2, x.len.2, col='white',space=0, border='white',add=T)
        text(genes.dat$positions.genes, 0.5, genes.dat$gene.names, cex = 1, col = 'white')
        dev.off()

      }
  }
  return(list(genes.dat,plotting.dat))
}
