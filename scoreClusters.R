#' scoreClusters Function
#'
#' @param Dat The CRISPRTarget file that is loaded
#' @export
#' @examples
#' scoreClusters()

scoreClusters <- function(nrDat){
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
    clust.num <- c(1,1,1,1)
    phage.dat <- getSubsetBasedOnRows(clust.res, 10, i)
    phage.dat <- transform(phage.dat, Genome_length = as.character(Genome_length))
    phage.dat <- transform(phage.dat, Genome_length = as.numeric(Genome_length))
    phage.dat <- phage.dat[order(phage.dat$Protospacer_start),]
    phage.dat <- transform(phage.dat, gap = as.character(gap))
    phage.dat <- transform(phage.dat, gap = as.numeric(gap))
    yy <- phage.dat$Genome_length[1]
    if(yy > 0){
      xx <- length(phage.dat$Spacer_ID)-3
      for(j in 1:xx){
        if(j > 0){
          rn0 <- as.character(rownames(phage.dat[j,]))
          rn1 <- as.character(rownames(phage.dat[j+1,]))
          rn2 <- as.character(rownames(phage.dat[j+2,]))
          rn3 <- as.character(rownames(phage.dat[j+3,]))
          #print(phage.dat[j+3,1])
          x <- c(phage.dat$gap[j],phage.dat$gap[j + 1],phage.dat$gap[j + 2])
          y1 <- runif(10000, min=0, max=yy)
          y <- vector(mode="numeric", length=10000)
          for(i in 1:(length(y1)-3)){
            a1 <- abs(y1[i]-y1[i+1])
            a2 <- abs(y1[i+1]-y1[i+2])
            a3 <- abs(y1[i+2]-y1[i+3])
            y[i] <- mean(a1 + a2 +a3)
          }
          ks.res <- suppressWarnings( ks.test(x,y1,alternative = 't', exact = F))
          ks.pval <- ks.res$p.value

          x0 <- clust.res[rn0,12]
          x1 <- clust.res[rn1,12]
          x2 <- clust.res[rn2,12]
          x3 <- clust.res[rn3,12]
          if(x0 > ks.pval){
            clust.res[rn0,12] <- ks.pval
          }
          if(x1 > ks.pval){
            clust.res[rn1,12] <- ks.pval
          }
          if(x2 > ks.pval){
            clust.res[rn2,12] <- ks.pval
          }
          #print(x3)
          #print(ks.pval)
          #if(x3 > ks.pval){
            #clust.res[rn3,12] <- ks.pval
          #}
          for(k in 1:4){
            if(ks.pval < 0.01+0.04*k){
              if(clust.res[rn0, 13+k]==0){
                clust.res[rn0, 13+k] <- clust.num[k]
                clust.res[rn1, 13+k] <- clust.num[k]
                clust.res[rn2, 13+k] <- clust.num[k]
                clust.res[rn3, 13+k] <- clust.num[k]
                clust.num[k] <- clust.num[k] + 1
              }else{
                clust.res[rn1, 13+k]  <- clust.res[rn0, 13+k]
                clust.res[rn2, 13+k]  <- clust.res[rn0, 13+k]
                clust.res[rn3, 13+k]  <- clust.res[rn0, 13+k]
              }
            }
          }
        }
      }
    }
  }

  clust.res[,18] <- p.adjust(clust.res$p.value, method = 'BH', n=length(clust.res$p.value))

  return(clust.res)
}

