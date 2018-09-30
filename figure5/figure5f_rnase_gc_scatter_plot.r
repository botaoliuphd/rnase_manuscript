#module load R/3.4.0
#R
#### Scatter plots between monosomal RNA amounts and GC contents or peak RPF lengths ####

base <- read.table("hip_rnase_titration_gc.txt",sep="\t",header=T)

pdf("hip_rnasese_titration_gc_scatter.pdf",width=7,height=6)
par(pin=c(4,3.3))
  
r1 <- format(round(cor(base[,"rnase"], base[,"gc"]),3), nsmall = 2)
r2 <- format(round(cor(base[,"rnase"], base[,"length"]),3), nsmall = 2)

plot(base[,"rnase"], base[,"gc"], ylim=c(48,60), axes=FALSE, pch=17,cex.lab=1.5,cex.axis=1.5,col="black", xlab="", ylab="")
axis(2, ylim=c(48,60),col="black",las=1,cex.axis=1.2)
mtext("GC%",side=2,line=2.5,cex=1.5)
abline(lm(base[,"gc"]~base[,"rnase"]), col="black")
box()

par(new=TRUE)

plot(base[,"rnase"], base[,"length"], ylim=c(27,34),axes=FALSE, pch=17,cex.lab=1.5,cex.axis=1.5,col="red", xlab="", ylab="")
axis(4, ylim=c(27,34),at=27:34,col="red",col.axis="red",las=1,cex.axis=1.2)
mtext("RPF length (nt)",side=4,col="red",line=3,cex=1.5)
abline(lm(base[,"length"]~base[,"rnase"]), col="red")

axis(1,xlim=c(-2,6),col="black",las=1,cex.axis=1.2)
mtext("log5 RNase concentration",side=1,col="black",line=2.5,cex=1.5)  

legend("topright", legend = c(paste("r=",r1,sep=""),paste("r=",r2,sep="")),text.col=c("black","red"),bty = "n",cex=1.5, text.font=3)
dev.off()
