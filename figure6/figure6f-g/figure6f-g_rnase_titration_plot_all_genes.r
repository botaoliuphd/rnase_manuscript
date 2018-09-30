# The most abundant isoform in rnase4 sample estimated by RSEM is selected as the representative transcript.
# perl select_abundant_isoform.pl
# Calculate the RPF offsets and frame preferences with Plastid.
# Extract the counts at each nucleotide position using P-sites of RPFs with variable offsets across the transcriptome.
# Follow the workflow described in psite_frame_get_count_plastid.bash

### Plot RPF distributon on the representivie mRNA isoforms ###
# module load R/3.2.2
# bsub -o log_plot -q long -n 1 -W 1000 -R rusage[mem=20000] "module load R/3.2.2; R CMD BATCH rnase_titration_plot_all_genes.r"

library(zoo)

#If file doesn't exist or no line, store NULL and print out warning but do not terminate
read.table.safe<-function(file,as.is=T,warning=F,...){
 if(file.exists(file) && length(system(paste("cat ",file,sep=""),intern=T))>0) {
   return(read.table(file,as.is=as.is,...))
 }
 else {
   if(warning==T) warning(paste("file ",file," is read as NULL variable.",sep=""))
   return(NULL)
 }
}

r=50

# Chosse the most abundant isoform in rnase4 sample as the representive transcript.
  
  proteinpos_file="rnase4_isoform_filtered_summary.txt"
  proteinpos<-read.table(proteinpos_file,sep="\t",header=T)
  row.names(proteinpos) = as.character(proteinpos[,1])

  filedir1="rnase1"
  fileprefix1="rnase1_"
  filesuffix1="txt"
  
  filedir2="rnase2"
  fileprefix2="rnase2_"
  filesuffix2="txt"
  
  filedir3="rnase3"
  fileprefix3="rnase3_"
  filesuffix3="txt"
  
  filedir4="rnase4"
  fileprefix4="rnase4_"
  filesuffix4="txt"

  filedir5="rnase5"
  fileprefix5="rnase5_"
  filesuffix5="txt"

pdf("rnase_titration.r50.all.plots.pdf",width=8.5,height=12,pointsize=15)

	for(mRNA in rownames(proteinpos)[order(proteinpos[,2])]){
    startcodon = as.numeric(proteinpos[mRNA,3])+1
	stopcodon = as.numeric(proteinpos[mRNA,4])
	
    f1=paste(filedir1,"/",fileprefix1,mRNA,".",filesuffix1,sep="");
    bucket1 = read.table.safe(f1,header=F)
    if(is.null(bucket1)) next
    if(!is.na(r)){ if (sum(bucket1[startcodon:stopcodon,1])<r) next } 
	
	f2=paste(filedir2,"/",fileprefix2,mRNA,".",filesuffix2,sep="");
    bucket2 = read.table.safe(f2,header=F)
    if(is.null(bucket2)) next
    if(!is.na(r)){ if (sum(bucket2[startcodon:stopcodon,1])<r) next } 
	
	f3=paste(filedir3,"/",fileprefix3,mRNA,".",filesuffix3,sep="");
    bucket3 = read.table.safe(f3,header=F)
    if(is.null(bucket3)) next
    if(!is.na(r)){ if (sum(bucket3[startcodon:stopcodon,1])<r) next } 

	f4=paste(filedir4,"/",fileprefix4,mRNA,".",filesuffix4,sep="");
    bucket4 = read.table.safe(f4,header=F)
    if(is.null(bucket4)) next
    if(!is.na(r)){ if (sum(bucket4[startcodon:stopcodon,1])<r) next } 

	f5=paste(filedir5,"/",fileprefix5,mRNA,".",filesuffix5,sep="");
    bucket5 = read.table.safe(f5,header=F)
    if(is.null(bucket5)) next
    if(!is.na(r)){ if (sum(bucket5[startcodon:stopcodon,1])<r) next } 

	# Normalize the counts to the mean count on the CDS.
	gene=proteinpos[as.character(mRNA),2]
	norm.bucket1 = as.matrix(bucket1[,1]/mean(bucket1[startcodon:stopcodon,1]))
	norm.bucket2 = as.matrix(bucket2[,1]/mean(bucket2[startcodon:stopcodon,1]))
	norm.bucket3 = as.matrix(bucket3[,1]/mean(bucket3[startcodon:stopcodon,1]))
	norm.bucket4 = as.matrix(bucket4[,1]/mean(bucket4[startcodon:stopcodon,1]))
	norm.bucket5 = as.matrix(bucket5[,1]/mean(bucket5[startcodon:stopcodon,1]))


  plot.col=4
  axis.col=rgb(0.7,0.7,0.7)
  
  layout(matrix(c(0,1:7,0),nrow=9,ncol=1),width=c(1),height=c(0.1,0.5,1.2,1.2,1.2,1.2,1.2,0.5,0.1))
  par(mar=c(0,8,0,5),lend=2)

  plot.new(); plot.window(c(0,1),c(0,1)); text(0.5,0.5,paste(gene,mRNA,sep=", "),cex=2)
  plot(1:nrow(norm.bucket1),as.numeric(norm.bucket1[,1]),xlim=c(1,nrow(norm.bucket1)),type="l",col=plot.col,axes=F,ylab="RNase1",cex.lab=2,ylim=c(0,max(as.numeric(norm.bucket1[,1]))*1.2));box(lwd=1/3,col=axis.col);axis(2,lwd=1/3,col=axis.col)
  plot(1:nrow(norm.bucket2),as.numeric(norm.bucket2[,1]),xlim=c(1,nrow(norm.bucket2)),type="l",col=plot.col,axes=F,ylab="RNase2",cex.lab=2,ylim=c(0,max(as.numeric(norm.bucket2[,1]))*1.2));box(lwd=1/3,col=axis.col);axis(2,lwd=1/3,col=axis.col)
  plot(1:nrow(norm.bucket3),as.numeric(norm.bucket3[,1]),xlim=c(1,nrow(norm.bucket3)),type="l",col=plot.col,axes=F,ylab="RNase3",cex.lab=2,ylim=c(0,max(as.numeric(norm.bucket3[,1]))*1.2));box(lwd=1/3,col=axis.col);axis(2,lwd=1/3,col=axis.col)
  plot(1:nrow(norm.bucket4),as.numeric(norm.bucket4[,1]),xlim=c(1,nrow(norm.bucket4)),type="l",col=plot.col,axes=F,ylab="RNase4",cex.lab=2,ylim=c(0,max(as.numeric(norm.bucket4[,1]))*1.2));box(lwd=1/3,col=axis.col);axis(2,lwd=1/3,col=axis.col)
  plot(1:nrow(norm.bucket5),as.numeric(norm.bucket5[,1]),xlim=c(1,nrow(norm.bucket5)),type="l",col=plot.col,axes=F,ylab="RNase5",cex.lab=2,ylim=c(0,max(as.numeric(norm.bucket5[,1]))*1.2));box(lwd=1/3,col=axis.col);axis(2,lwd=1/3,col=axis.col)
  axis(1,lwd=1/3,col=axis.col)
  plot(startcodon,0.5,cex=2,type="p",pch=24,ylim=c(0,1),axes=F,col=rgb(0,0.5,0),bg=rgb(0,0.5,0),xlim=c(1,nrow(norm.bucket4)),ylab="")
  points(stopcodon,0.5,cex=2,col=rgb(0.5,0,0),pch=24,bg=rgb(0.5,0,0))
  }
  dev.off()
  