#module load R/3.4.0
#R
#### Scatter plots between monosomal RNA amounts and GC contents or peak RPF lengths ####

base <- read.table("ipsc_lower_rnase_gc.txt",sep="\t",header=T)
log2_rna = log(base[,"rna"])/log(2)

pdf("ipsc_lower_rnase_gc_scatter.pdf",width=7,height=6)
par(pin=c(4,3.3))
  
r1 <- format(round(cor(log2_rna, base[,"gc"]),3), nsmall = 2)
r2 <- format(round(cor(log2_rna, base[,"length"]),3), nsmall = 2)

plot(log2_rna, base[,"gc"], axes=FALSE, pch=17,cex.lab=1.5,cex.axis=1.5,col="black", xlab="", ylab="",xlim=c(1.8,4),ylim=c(49,56))
axis(2, ylim=c(49,56),col="black",las=1,cex.axis=1.2)
mtext("GC%",side=2,line=2.5,cex=1.5)

abline(lm(base[,"gc"]~log2_rna), col="black")
box()

par(new=TRUE)

plot(log2_rna, base[,"length"], axes=FALSE, pch=17,cex.lab=1.5,cex.axis=1.5,col="red", xlab="", ylab="",xlim=c(1.8,4),ylim=c(27,34))
axis(4, ylim=c(27,34),at=27:34,col="red",col.axis="red",las=1,cex.axis=1.2)
mtext("RPF length (nt)",side=4,col="red",line=3,cex=1.5)
abline(lm(base[,"length"]~log2_rna), col="red")

x <- (-1):4
base_num <- 2
axis(1,at=x,labels = base_num^(x), col="black",las=1,cex.axis=1.2)
mtext("RNA amount (ug)",side=1,col="black",line=2.5,cex=1.5)  

legend("topright", legend = c(paste("r=",r1,sep=""),paste("r=",r2,sep="")),text.col=c("black","red"),bty = "n",cex=1.5, text.font=3)
dev.off()
