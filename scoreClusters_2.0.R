#' scoreClusters_2.0 Function
#'
#' @param Dat The CRISPRTarget file that is loaded
#' @export
#' @examples
#' scoreClusters_2.0()

scoreClusters_2.0 <- function(nrDat){
  x <- vector(mode="numeric", length=length(nrDat$Score))
  x1<- x+1
  clust.res <- cbind(nrDat, x1,x,x,x,x,x,x1)
  colnames(clust.res) <- c(colnames(nrDat),'p-value','gap','cl1','cl2','cl3','cl4','adj.p.value')
  clust.res <- clust.res[order(clust.res$Protospacer_start),]
  clust.res <- clust.res[order(clust.res$Target_name),]

  for(i in levels(clust.res$Target_name)){
    arrayDat <- getSubsetBasedOnRows(clust.res,10,i)
    llm <- length(arrayDat$Score)-1
    next_row <- as.character(rownames(arrayDat[1,]))
    for(l in 1:llm){
      row_num <- as.character(rownames(arrayDat[l,]))
      next_row <- as.character(rownames(arrayDat[l+1,]))
      clust.res[row_num,13] <- abs(arrayDat$Protospacer_start[l+1] - arrayDat$Protospacer_stop[l])
    }
    l <- length(arrayDat$Score)
    row_num <- as.character(rownames(arrayDat[l,]))
    clust.res[row_num,13] <- abs(arrayDat$Protospacer_start[l] - arrayDat$Protospacer_stop[l-1])
  }
  clust.res <- transform(clust.res, gap = as.character(gap))
  clust.res <- transform(clust.res, gap = as.numeric(gap))
  for(i in levels(clust.res$Target_name)){
    phage.dat <- getSubsetBasedOnRows(clust.res, 10, i)
    phage.dat <- transform(phage.dat, Genome_length = as.character(Genome_length))
    phage.dat <- transform(phage.dat, Genome_length = as.numeric(Genome_length))
    phage.dat <- phage.dat[order(phage.dat$Protospacer_start),]
    phage.dat <- transform(phage.dat, gap = as.character(gap))
    phage.dat <- transform(phage.dat, gap = as.numeric(gap))
    yy <- phage.dat$Genome_length[1]
    if(yy > 0){
      row_names <- vector(mode="character", length=length(phage.dat$Spacer_ID))
      x <- vector(mode="numeric", length=length(phage.dat$Spacer_ID))
      for(j in 1:length(phage.dat$Spacer_ID)){
          row_names[j] <- as.character(rownames(phage.dat[j,]))
          x[j] <- phage.dat$gap[j]
      }
          y <- vector(mode="numeric", length=10000)
          for(jj in 1:10000){
            y1 <- runif(2, min=0, max=yy)
            y[jj] <- abs(y1[1]-y1[2])
          }
          ks.res <- suppressWarnings( ks.test(x,y,alternative = 't', exact = F))
          for(k in 1:length(phage.dat$gap)){

            clust.res[row_names[k],12] <- ks.res$p.value
          if(ks.res$p.value < 0.01){
            clust.res[row_names[k],14] <- 1
          }
          if(ks.res$p.value < 0.05){
            clust.res[row_names[k],15] <- 1
          }
          if(ks.res$p.value < 0.1){
            clust.res[row_names[k],16] <- 1
          }
          if(ks.res$p.value < 0.15){
            clust.res[row_names[k],17] <- 1
          }
          }
          }
          }





clust.res[,18] <- p.adjust(clust.res$p.value, method = 'BH', n=length(clust.res$p.value))

return(clust.res)
}

