# module load R/3.4.0

library("DESeq2")

#########################
####### Batch1-2 ########
#########################
all_data <-read.table("fmr1.slice.rpf.all.rsem.gene.summary.counts.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)

dataRPF <- data.frame(all_data[,c(
"wt_rpf_batch1_rep1", "wt_rpf_batch2_rep1", 
"wt_rpf_batch1_rep2", "wt_rpf_batch2_rep2", 
"ko_rpf_batch1_rep1",  
"ko_rpf_batch1_rep2", "ko_rpf_batch2_rep1"
)])
colsRPF <- c(1:7)
dataRPF[,colsRPF] <- apply(dataRPF[,colsRPF], 2, function(x) as.numeric(as.integer(x)))
genotype <- factor( c (
"WT","WT",
"WT","WT",
"KO",
"KO","KO"))
batch <- factor( c (
"Batch1","Batch2",
"Batch1","Batch2",
"Batch1",
"Batch1","Batch2"))
colDataRPF <- as.data.frame((colnames(dataRPF)))
colDataRPF <- cbind(colDataRPF, genotype,batch)
colnames(colDataRPF) <- c("libs","genotype","batch")
ddsRPF <- DESeqDataSetFromMatrix(countData=as.matrix(dataRPF),colData=colDataRPF, design = ~ genotype)
ddsRPF$genotype = relevel(ddsRPF$genotype,"WT")
ddsRPF$batch = relevel(ddsRPF$batch,"Batch1")

# pre-filtering
ddsRPF <- estimateSizeFactors(ddsRPF)
mnc <- rowMeans(counts(ddsRPF, normalized=TRUE))
ddsRPF <- ddsRPF[mnc > 5,]

rld <- rlogTransformation(ddsRPF, blind=FALSE)

#pca plot
library(ggplot2)
library(stats)
library(ggbiplot)

ntop=1000

rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(ntop)]
pcaData = prcomp(t(assay(rld)[select,]),center = TRUE,scale. = TRUE)

#summary(pcaData)

output_name="fmr1ko.slice.batch1-2.rpf.pca.pdf"
pdf(output_name,width=5,height=5)
par(mar=c(5,5,5,5),cex=20)

g <- ggbiplot(pcaData, choices = 1:2, obs.scale = 1, var.scale = 1, 
              groups = colDataRPF$genotype, varname.size=0, var.axes = F, ellipse = TRUE, circle = FALSE)+
			  geom_point(aes(colour=genotype, shape=batch), size = 3)+
			  scale_color_manual(name="Genotype", values=c("#FC9272","#9ECAE1")) +  
			  scale_shape_manual(name="Batch", values=c(15,17)) +
              theme(legend.direction = 'horizontal',legend.position = "top", panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
			  scale_y_continuous(limits = c(-40,40))+
			  scale_x_continuous(limits = c(-50,50))
print(g)
dev.off()

#### full model batch effect correction ####
library(sva)
rawCounts <- as.matrix(dataRPF)
mode(rawCounts)<-"numeric"
meanc <- rowMeans(rawCounts)
filtc <- subset(rawCounts, meanc > 5)
filtcnames <- rownames(filtc)
gc <- as.matrix(read.table("wt_rpf_batch1_rep1.genes.results", row.names=1, header=T))
mode(gc)<-"numeric"
filtc_len <- gc[filtcnames, "length"]
normTPM <- apply(filtc,2,function(x) (x/filtc_len)/(sum(x/filtc_len))*1000000) 
total_perbase <- apply(filtc,2,function(x) sum(x/filtc_len))
log2filtc <- log2(normTPM+1)
meta <- colDataRPF
model <- model.matrix(~as.factor(genotype),data=meta) # used full model
combat <- ComBat(dat=log2filtc,batch=meta$batch,mod=model,par.prior=TRUE,prior.plots=FALSE)
tpm_batchN <- 2^(combat)
mult <- sweep(tpm_batchN,MARGIN=2,total_perbase,'*')
counts_batchN <- apply(mult,2,function(x) (x/1000000)*filtc_len)
write.table(counts_batchN, "2batchN.full_model.fmr1ko.slice.batch1-2.rpf.counts.txt", sep="\t")
write.table(tpm_batchN, "2batchN.full_model.fmr1ko.slice.batch1-2.rpf.tpm.txt", sep="\t")


#########after batch correction#########
all_data <-read.table("2batchN.full_model.fmr1ko.slice.batch1-2.rpf.counts.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)

dataRPF <- data.frame(all_data[,c(
"wt_rpf_batch1_rep1", "wt_rpf_batch2_rep1", 
"wt_rpf_batch1_rep2", "wt_rpf_batch2_rep2", 
"ko_rpf_batch1_rep1",  
"ko_rpf_batch1_rep2", "ko_rpf_batch2_rep1"
)])
colsRPF <- c(1:7)
dataRPF[,colsRPF] <- apply(dataRPF[,colsRPF], 2, function(x) as.numeric(as.integer(x)))
genotype <- factor( c (
"WT","WT",
"WT","WT",
"KO",
"KO","KO"))
batch <- factor( c (
"Batch1","Batch2",
"Batch1","Batch2",
"Batch1",
"Batch1","Batch2"))
colDataRPF <- as.data.frame((colnames(dataRPF)))
colDataRPF <- cbind(colDataRPF, genotype,batch)
colnames(colDataRPF) <- c("libs","genotype","batch")
ddsRPF <- DESeqDataSetFromMatrix(countData=as.matrix(dataRPF),colData=colDataRPF, design = ~ genotype)
ddsRPF$genotype = relevel(ddsRPF$genotype,"WT")
ddsRPF$batch = relevel(ddsRPF$batch,"Batch1")

# pre-filtering
ddsRPF <- estimateSizeFactors(ddsRPF)
mnc <- rowMeans(counts(ddsRPF, normalized=TRUE))
ddsRPF <- ddsRPF[mnc > 5,]

rld <- rlogTransformation(ddsRPF, blind=FALSE)

#pca plot
library("ggplot2")
library("stats")
library(ggbiplot)

ntop=1000

rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(ntop)]
pcaData = prcomp(t(assay(rld)[select,]),center = TRUE,scale. = TRUE)

#summary(pcaData)

output_name="2batchN.fmr1ko.slice.batch1-2.rpf.pca.pdf"
pdf(output_name,width=5,height=5)
par(mar=c(5,5,5,5),cex=20)

g <- ggbiplot(pcaData, choices = 1:2, obs.scale = 1, var.scale = 1, 
              groups = colDataRPF$genotype, varname.size=0, var.axes = F, ellipse = TRUE, circle = FALSE)+
			  geom_point(aes(colour=genotype, shape=batch), size = 3)+
			  scale_color_manual(name="Genotype", values=c("#FC9272","#9ECAE1")) +  
			  scale_shape_manual(name="Batch", values=c(15,17)) +
              theme(legend.direction = 'horizontal',legend.position = "top", panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
			  scale_y_continuous(limits = c(-40,40))+
			  scale_x_continuous(limits = c(-50,50))
print(g)
dev.off()


## differential gene expression analysis ##

dds_all <- DESeq(ddsRPF)


##################################
##			WT vs KO	        ##
##################################  

res <- results(dds_all,contrast=c("genotype","KO", "WT"),alpha=0.1)

library("AnnotationDbi")
library("org.Mm.eg.db")

res$symbol <- rownames(res)

res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
					 
res$genename <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")
					 
res<-res[c(7,8,9,1,2,3,4,5,6)]
res_ordered <- res[order(res$log2FoldChange),]

write.table(res_ordered, "2batchN.fmr1ko.slice.batch1-2.ko.vs.wt.rpf.all.txt", sep = "\t", row.names = FALSE)

res_ordered_tmp <- res_ordered[!is.na(res_ordered$padj)& !is.na(res_ordered$log2FoldChange),]

res_ordered_sig <- res_ordered_tmp[(res_ordered_tmp$padj<0.1),]

write.table(res_ordered_sig, "2batchN.fmr1ko.slice.batch1-2.ko.vs.wt.rpf.adjp.0.1.txt", sep = "\t", row.names = FALSE)

res_ordered_sig <- res_ordered_tmp[(res_ordered_tmp$padj<0.05),]

write.table(res_ordered_sig, "2batchN.fmr1ko.slice.batch1-2.ko.vs.wt.rpf.adjp.0.05.txt", sep = "\t", row.names = FALSE)


############ volcano plot ##############
output_name="2batchN.fmr1ko.slice.batch1-2.rpf.volcano.pdf"
pdf(output_name,width=10,height=10)
par(mar=c(5,5,5,5), cex=1.7, cex.main=1.7, cex.axis=1.7, cex.lab=1.7)

topT <- as.data.frame(res)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, cex=1.0, ylim=c(0,20), xlim=c(-1,1),xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~adjusted~p~value)))
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#Add lines P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
dev.off()
###########################################



#########################
####### Batch3-4 ########
#########################

all_data <-read.table("fmr1.slice.rpf.all.rsem.gene.summary.counts.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)

dataRPF <- data.frame(all_data[,c(
"wt_rpf_batch3_rep1", "wt_rpf_batch4_rep1", 
"wt_rpf_batch3_rep2", "wt_rpf_batch4_rep2", 
"ko_rpf_batch3_rep1",  "ko_rpf_batch4_rep1",
"ko_rpf_batch3_rep2", "ko_rpf_batch4_rep2"
)])
colsRPF <- c(1:8)
dataRPF[,colsRPF] <- apply(dataRPF[,colsRPF], 2, function(x) as.numeric(as.integer(x)))
genotype <- factor( c (
"WT","WT",
"WT","WT",
"KO","KO",
"KO","KO"))
batch <- factor( c (
"Batch3","Batch4",
"Batch3","Batch4",
"Batch3","Batch4",
"Batch3","Batch4"))
colDataRPF <- as.data.frame((colnames(dataRPF)))
colDataRPF <- cbind(colDataRPF, genotype,batch)
colnames(colDataRPF) <- c("libs","genotype","batch")
ddsRPF <- DESeqDataSetFromMatrix(countData=as.matrix(dataRPF),colData=colDataRPF, design = ~ genotype)
ddsRPF$genotype = relevel(ddsRPF$genotype,"WT")
ddsRPF$batch = relevel(ddsRPF$batch,"Batch3")

# pre-filtering
ddsRPF <- estimateSizeFactors(ddsRPF)
mnc <- rowMeans(counts(ddsRPF, normalized=TRUE))
ddsRPF <- ddsRPF[mnc > 5,]

rld <- rlogTransformation(ddsRPF, blind=FALSE)

#pca plot
library(ggplot2)
library(stats)
library(ggbiplot)

ntop=1000

rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(ntop)]
pcaData = prcomp(t(assay(rld)[select,]),center = TRUE,scale. = TRUE)

#summary(pcaData)

output_name="fmr1ko.slice.batch3-4.rpf.pca.pdf"
pdf(output_name,width=5,height=5)
par(mar=c(5,5,5,5),cex=20)

g <- ggbiplot(pcaData, choices = 1:2, obs.scale = 1, var.scale = 1, 
              groups = colDataRPF$genotype, varname.size=0, var.axes = F, ellipse = TRUE, circle = FALSE)+
			  geom_point(aes(colour=genotype, shape=batch), size = 3)+
			  scale_color_manual(name="Genotype", values=c("#FC9272","#9ECAE1")) +  
			  scale_shape_manual(name="Batch", values=c(15,17)) +
              theme(legend.direction = 'horizontal',legend.position = "top", panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
			  scale_y_continuous(limits = c(-40,40))+
			  scale_x_continuous(limits = c(-50,50))
print(g)
dev.off()


#### full model batch effect correction ####
library(sva)
rawCounts <- as.matrix(dataRPF)
mode(rawCounts)<-"numeric"
meanc <- rowMeans(rawCounts)
filtc <- subset(rawCounts, meanc > 5)
filtcnames <- rownames(filtc)
gc <- as.matrix(read.table("wt_rpf_batch3_rep1.genes.results", row.names=1, header=T))
mode(gc)<-"numeric"
filtc_len <- gc[filtcnames, "length"]
normTPM <- apply(filtc,2,function(x) (x/filtc_len)/(sum(x/filtc_len))*1000000) 
total_perbase <- apply(filtc,2,function(x) sum(x/filtc_len))
log2filtc <- log2(normTPM+1)
meta <- colDataRPF
model <- model.matrix(~as.factor(genotype),data=meta) # used full model
combat <- ComBat(dat=log2filtc,batch=meta$batch,mod=model,par.prior=TRUE,prior.plots=FALSE)
tpm_batchN <- 2^(combat)
mult <- sweep(tpm_batchN,MARGIN=2,total_perbase,'*')
counts_batchN <- apply(mult,2,function(x) (x/1000000)*filtc_len)
write.table(counts_batchN, "2batchN.full_model.fmr1ko.slice.batch3-4.rpf.counts.txt", sep="\t")
write.table(tpm_batchN, "2batchN.full_model.fmr1ko.slice.batch3-4.rpf.tpm.txt", sep="\t")


#########after batch correction#########
all_data <-read.table("2batchN.full_model.fmr1ko.slice.batch3-4.rpf.counts.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)

dataRPF <- data.frame(all_data[,c(
"wt_rpf_batch3_rep1", "wt_rpf_batch4_rep1", 
"wt_rpf_batch3_rep2", "wt_rpf_batch4_rep2", 
"ko_rpf_batch3_rep1",  "ko_rpf_batch4_rep1",
"ko_rpf_batch3_rep2", "ko_rpf_batch4_rep2"
)])
colsRPF <- c(1:8)
dataRPF[,colsRPF] <- apply(dataRPF[,colsRPF], 2, function(x) as.numeric(as.integer(x)))
genotype <- factor( c (
"WT","WT",
"WT","WT",
"KO","KO",
"KO","KO"))
batch <- factor( c (
"Batch3","Batch4",
"Batch3","Batch4",
"Batch3","Batch4",
"Batch3","Batch4"))
colDataRPF <- as.data.frame((colnames(dataRPF)))
colDataRPF <- cbind(colDataRPF, genotype,batch)
colnames(colDataRPF) <- c("libs","genotype","batch")
ddsRPF <- DESeqDataSetFromMatrix(countData=as.matrix(dataRPF),colData=colDataRPF, design = ~ genotype)
ddsRPF$genotype = relevel(ddsRPF$genotype,"WT")
ddsRPF$batch = relevel(ddsRPF$batch,"Batch3")

# pre-filtering
ddsRPF <- estimateSizeFactors(ddsRPF)
mnc <- rowMeans(counts(ddsRPF, normalized=TRUE))
ddsRPF <- ddsRPF[mnc > 5,]

rld <- rlogTransformation(ddsRPF, blind=FALSE)

#pca plot
library("ggplot2")
library("stats")
library(ggbiplot)

ntop=1000

rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(ntop)]
pcaData = prcomp(t(assay(rld)[select,]),center = TRUE,scale. = TRUE)

#summary(pcaData)

output_name="2batchN.fmr1ko.slice.batch3-4.rpf.pca.pdf"
pdf(output_name,width=5,height=5)
par(mar=c(5,5,5,5),cex=20)

g <- ggbiplot(pcaData, choices = 1:2, obs.scale = 1, var.scale = 1, 
              groups = colDataRPF$genotype, varname.size=0, var.axes = F, ellipse = TRUE, circle = FALSE)+
			  geom_point(aes(colour=genotype, shape=batch), size = 3)+
			  scale_color_manual(name="Genotype", values=c("#FC9272","#9ECAE1")) +  
			  scale_shape_manual(name="Batch", values=c(15,17)) +
              theme(legend.direction = 'horizontal',legend.position = "top", panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
			  scale_y_continuous(limits = c(-50,50))+
			  scale_x_continuous(limits = c(-50,50))
print(g)
dev.off()

## differential gene expression analysis ##
dds_all <- DESeq(ddsRPF)

##################################
##			WT vs KO	        ##
##################################  

res <- results(dds_all,contrast=c("genotype","KO", "WT"),alpha=0.1)

library("AnnotationDbi")
library("org.Mm.eg.db")

res$symbol <- rownames(res)

res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
					 
res$genename <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")
					 
res<-res[c(7,8,9,1,2,3,4,5,6)]
res_ordered <- res[order(res$log2FoldChange),]

write.table(res_ordered, "2batchN.fmr1ko.slice.batch3-4.ko.vs.wt.rpf.all.txt", sep = "\t", row.names = FALSE)

res_ordered_tmp <- res_ordered[!is.na(res_ordered$padj)& !is.na(res_ordered$log2FoldChange),]

res_ordered_sig <- res_ordered_tmp[(res_ordered_tmp$padj<0.1),]

write.table(res_ordered_sig, "2batchN.fmr1ko.slice.batch3-4.ko.vs.wt.rpf.adjp.0.1.txt", sep = "\t", row.names = FALSE)

res_ordered_sig <- res_ordered_tmp[(res_ordered_tmp$padj<0.05),]

write.table(res_ordered_sig, "2batchN.fmr1ko.slice.batch3-4.ko.vs.wt.rpf.adjp.0.05.txt", sep = "\t", row.names = FALSE)

############ volcano plot ##############
output_name="2batchN.fmr1ko.slice.batch3-4.rpf.volcano.pdf"
pdf(output_name,width=10,height=10)
par(mar=c(5,5,5,5), cex=1.7, cex.main=1.7, cex.axis=1.7, cex.lab=1.7)

topT <- as.data.frame(res)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, cex=1.0, ylim=c(0,20), xlim=c(-1,1),xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~adjusted~p~value)))
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#Add lines P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
dev.off()
###########################################