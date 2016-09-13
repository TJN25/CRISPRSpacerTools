#' clust.data.frame Function
#'
#' @param clust.res Data frame containing info about all the protospacers that are clusterd
#' @export
#' @examples
#' clust.data.frame()

clust.data.frame <- function(clust.res){
tt <- clust.res
##Assign unique identifier to each cluster
unique.cluster.names <- vector(mode="numeric", length=length(tt$Score))
for(i in 1:length(unique.cluster.names)){
  if(tt$cl1[i]>0){
  unique.cluster.names[i] <- paste(tt$Target_name[i],tt$cl1[i], sep='_')
  }else{
    unique.cluster.names[i] <- NA
  }
}
tt <- cbind(tt, unique.cluster.names)

##Create new data frame which will include a single row for each cluster
x <- vector(mode="numeric", length=length(levels(tt$unique.cluster.names)))
clust.dat <- data.frame(target.name=x,cluster=x,spacer.number=x,mean.gap=x,p.value=x,adj.p.value=x,seq.score=x,alternating.score=x,seq.strand.score=x,alt.strand.score=x,gap.score=x,sequential.score=x,alt.score=x)
a <- 0
for(i in levels(tt$Target_name)){
  xx <- getSubsetBasedOnRows(tt,10,i)
  xx <- transform(xx, cl1 = as.factor(cl1))
  for(j in levels(xx$cl1)){
    if(j != '0'){
      a <- a +1
      clusterDat <- getSubsetBasedOnRows(xx,14,j)
      clusterDat <- transform(clusterDat, gap = as.character(gap))
      clusterDat <- transform(clusterDat, gap = as.numeric(gap))
      clusterDat <- transform(clusterDat, p.value = as.character(p.value))
      clusterDat <- transform(clusterDat, p.value = as.numeric(p.value))
      clusterDat <- transform(clusterDat, adj.p.value = as.character(adj.p.value))
      clusterDat <- transform(clusterDat, adj.p.value = as.numeric(adj.p.value))
      spacers <- c()
      for(k in 1:length(clusterDat[,1])){
        spacers <- c(spacers,clusterDat[k,8])
      }
      ##score the cluster based on the order of the spacers, the length of the cluster and the strand
      scoringDat <- clusterDat
      scoringDat <- scoringDat[order(scoringDat$Spacer_Number),]


      seq.score.val <- 0
      alt.score.val <- 0
      strand.score.val <- 0
      alt.strand.val <- 0
      gaps <- scoringDat[order(scoringDat$gap),13]
      gap.score.val <- 1/(gaps[length(gaps)]-gaps[1])*8000
      #seq.strand.val calculation
      seq.score.val <- ifelse(scoringDat$Protospacer_start[1] < scoringDat$Protospacer_start[2],seq.score.val + 1, seq.score.val - 1)
      #alt.strand.val calculation
      if(scoringDat$Protospacer_start[1] < scoringDat$Protospacer_start[2]){
        alt.score.val <- alt.score.val + 1
        dir.val <- 1
        print.val <- paste('a < b and first dir is +ve. Score is ',alt.score.val,' after adding 1')
        #print(print.val)
      }else{
        alt.score.val <- alt.score.val + 1
        dir.val <- -1
        print.val <- paste('a > b and first dir is -ve. Score is ',alt.score.val,' after adding 1')
        #print(print.val)
      }
      ##seq.strand.score.val calcuation
      strand.score.val <- ifelse(scoringDat$Protospacer_Strand[1]=='+', strand.score.val +1, strand.score.val -1)
      ##alt.strand.val calculation
      if(scoringDat$Protospacer_Strand[1] == '+'){

          alt.strand.val <- alt.strand.val + 1
          strand.val <- 1
      }else{
          alt.strand.val <- alt.strand.val + 1
          strand.val <- -1
      }

      for(k in 2:(length(scoringDat$Spacer_Number)-1)){
        #seq.strand.val calculation
        seq.score.val <- ifelse(scoringDat$Protospacer_start[k] < scoringDat$Protospacer_start[k+1],seq.score.val + 1, seq.score.val - 1)
        #alt.strand.val calculation
        if(scoringDat$Protospacer_start[k] < scoringDat$Protospacer_start[k+1]){
          if(dir.val==-1){
            alt.score.val <- alt.score.val + 1
            dir.val <- 1
            print.val <- paste('a < b and previous dir was -ve. Score is ',alt.score.val,' after adding 1')
            #print(print.val)
          }else{
            alt.score.val <- alt.score.val - 1
            dir.val <- 1
            print.val <- paste('a < b and previous dir was +ve. Score is ',alt.score.val,' after subtracting 1')
            #print(print.val)
          }
        }else{
          if(dir.val==1){
            alt.score.val <- alt.score.val + 1
            dir.val <- -1
            print.val <- paste('a > b and previous dir was +ve. Score is ',alt.score.val,' after adding 1')
            #print(print.val)
          }else{
            alt.score.val <- alt.score.val - 1
            dir.val <- -1
            print.val <- paste('a > b and previous dir was -ve. Score is ',alt.score.val,' after subtracting 1')
            #print(print.val)
          }

        }
        ##seq.strand.score.val calcuation
        strand.score.val <- ifelse(scoringDat$Protospacer_Strand[k]=='+', strand.score.val +1, strand.score.val -1)
        ##alt.strand.val calculation

        if(scoringDat$Protospacer_Strand[k] == '+'){
          if(strand.val==-1){
            alt.strand.val <- alt.strand.val + 1
            strand.val <- 1
            #print.val <- paste('a < b and previous dir was -ve. Score is ',alt.strand.val,' after adding 1')
            #print(print.val)
          }else{
            alt.strand.val <- alt.strand.val - 1
            strand.val <- 1
            #print.val <- paste('a < b and previous dir was +ve. Score is ',alt.score.val,' after subtracting 1')
            #print(print.val)
          }
        }else{
          if(strand.val==1){
            alt.strand.val <- alt.strand.val + 1
            strand.val <- -1
            print.val <- paste('a > b and previous dir was +ve. Score is ',alt.score.val,' after adding 1')
            #print(print.val)
          }else{
            alt.strand.val <- alt.strand.val - 1
            strand.val <- -1
            print.val <- paste('a > b and previous dir was -ve. Score is ',alt.score.val,' after subtracting 1')
            #print(print.val)
          }
        }
      }
      ##seq.strand.score.val calcuation
      strand.score.val <- ifelse(scoringDat$Protospacer_Strand[length(scoringDat$Protospacer_Strand)]=='+', strand.score.val +1, strand.score.val -1)
      ##alt.strand.val calculation
      if(scoringDat$Protospacer_Strand[length(scoringDat$Protospacer_Strand)] == '+'){
        if(strand.val==-1){
          alt.strand.val <- alt.strand.val + 1
          strand.val <- 1
          #print.val <- paste('a < b and previous dir was -ve. Score is ',alt.strand.val,' after adding 1')
          #print(print.val)
        }else{
          alt.strand.val <- alt.strand.val - 1
          strand.val <- 1
          #print.val <- paste('a < b and previous dir was +ve. Score is ',alt.score.val,' after subtracting 1')
          #print(print.val)
        }
      }else{
        if(strand.val==1){
          alt.strand.val <- alt.strand.val + 1
          strand.val <- -1
          print.val <- paste('a > b and previous dir was +ve. Score is ',alt.score.val,' after adding 1')
          #print(print.val)
        }else{
          alt.strand.val <- alt.strand.val - 1
          strand.val <- -1
          print.val <- paste('a > b and previous dir was -ve. Score is ',alt.score.val,' after subtracting 1')
          #print(print.val)
        }
      }

      seq.score.val <- abs(seq.score.val)
      strand.score.val <- abs(strand.score.val)

      spacers <- toString(spacers)
      clust.dat[a,1] <- i
      clust.dat[a,2] <- j
      clust.dat[a,3] <- spacers
      clust.dat[a,4] <- mean(clusterDat$gap)
      clust.dat[a,5] <- mean(clusterDat$p.value)
      clust.dat[a,6] <- mean(clusterDat$adj.p.value)
      clust.dat[a,7] <- seq.score.val
      clust.dat[a,8] <- alt.score.val
      clust.dat[a,9] <- strand.score.val
      clust.dat[a,10] <- alt.strand.val
      clust.dat[a,11] <- gap.score.val

      seq.score.total <- seq.score.val + strand.score.val + gap.score.val
      alt.score.total <- alt.strand.val + alt.score.val + gap.score.val
      clust.dat[a,12] <- seq.score.total
      clust.dat[a,13] <- alt.score.total
    }
  }
}
clust.dat <- transform(clust.dat, target.name = as.factor(target.name))
return(clust.dat)
}
