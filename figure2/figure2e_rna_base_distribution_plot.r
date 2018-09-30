#module load R/3.4.0
#R
#### Base composition plots ####

base <- read.table("qc_sample_wt_rna_batch1_rep1_reads_overlap_CDS.base",sep="\t",header=T)

base <- base[1:75,2:5]

base_norm <- 100* base/rowSums(base)

base_norm_gc <- rowSums(base_norm[,3:4])

gc_perc <- mean(base_norm_gc [10:65])

gc_perc <- format(round(gc_perc, 2), nsmall = 2)

pdf("qc_sample_wt_rna_batch1_rep1_reads_overlap_CDS.base.pdf",width=6,height=6)
par(pin=c(4,3.3))
  
tick1=seq(0, 75, 15)
x1=seq(0, 75, 15)

ymax = 2*mean(as.matrix(base_norm))

matplot(base_norm, type = c("l","l","l","l"),lty=c("solid","solid","solid","solid"),lwd=2.5,col = c( "#FC9272", "#99000D","#9ECAE1","#084594"),ylim=c(0,ymax),xlim=c(0,75),ylab="Percentage",xlab="Nucleotide position",xaxt = "n",cex.lab=1.5,cex.axis=1.5);axis(1, at = tick1, labels = as.character(x1), cex.axis = 1.5);legend("topright", legend = c("A","T","C","G"), col=c( "#FC9272", "#99000D","#9ECAE1","#084594"), lty = c("solid","solid","solid","solid"),bty = "n",cex=1.2,lwd=2.5);  
title(main=paste("GC%=",gc_perc,sep=""), line = 0.5,cex.main=1.5,font.main = 1)

dev.off()
