#' plotProtospacerPoints Function
#'
#' @param singleDat Dataframe for plotting
#' @export
#' @examples
#' plotProtospacerPoints()

plotProtospacerPoints <- function(singleDat){
  xx <- as.numeric(as.character(singleDat$Genome_length))/3.85
plot(singleDat$Protospacer_start,singleDat$Spacer_Number,
      main=NA,
     xlab='Target genome',
     ylab='Spacer number',
     type='p',
     pch=ifelse(is.na(singleDat$Matching_pam_in_3p_of_forward_strand),17,16),
     col=ifelse(singleDat$Protospacer_Strand=='+','red','blue'),
     xlim=c(0,as.numeric(as.character(singleDat$Genome_length[1]))))
legend("bottomleft",
       text.width=c(xx,xx,xx,xx),
       inset = c(0, -0.25),
       x.intersp=1,
       legend=c('Includes a PAM', 'Missing a PAM', 'Positive Strand', 'Negative Strand'),
       pch=c(16,17,15,15),
       col=c('black','black', 'red','blue'),bty='n', cex=0.8,
       xpd=T,
       horiz=TRUE)
text(singleDat$Protospacer_start,singleDat$Spacer_Number, labels=singleDat$Score, cex= 0.7, pos=4)
text(singleDat$Protospacer_start,singleDat$Spacer_Number, labels=singleDat$Spacer_Number, cex= 0.7, pos=2)

}
