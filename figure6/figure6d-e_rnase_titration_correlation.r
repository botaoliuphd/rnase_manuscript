#module load R/3.4.0

#R
library("DESeq2")
rsemRPF <-read.table("rnase.titration.ribo.all.rsem.gene.summary.counts.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)
dataRPF <- data.frame(rsemRPF[,c("rnase1", "rnase2", "rnase3", "rnase4", "rnase5")])
colsRPF <- c(1:5)
dataRPF[,colsRPF] <- apply(dataRPF[,colsRPF], 2, function(x) as.numeric(as.integer(x)))
colDataRPF <- as.data.frame(c("RNase1", "RNase2", "RNase3", "RNase4", "RNase5"))
colnames(colDataRPF) <- c("Libs")
sumdRPF <- rowSums(dataRPF)
filtdRPF <- subset(dataRPF, sumdRPF > 10)
ddsRPF <- DESeqDataSetFromMatrix(countData=as.matrix(filtdRPF),colData=colDataRPF, design = ~1)


rld <- rlogTransformation(ddsRPF, blind=FALSE)


#scatter plot
output_name="conc_titration_scatterplot_of_transformed_counts.pdf"
pdf(output_name,width=10,height=10)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(assay(rld)[,1:5], c("Conc.1", "Conc.2", "Conc.3", "Conc.4", "Conc.5"),pch = 16, lower.panel = panel.smooth, upper.panel = panel.cor, cex.labels = 3, font.labels = 1,xaxt = "n", yaxt = "n")

dev.off()


#pca plot
library(ggplot2)
library(stats)
library(ggbiplot)
library(RColorBrewer)


ntop=1000

rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(ntop)]
pcaData = prcomp(t(assay(rld)[select,]),center = TRUE,scale. = TRUE)

#summary(pcaData)

sequential <- brewer.pal(6, "Blues")

output_name="rnase_titration_pca.pdf"
pdf(output_name,width=5,height=5)
par(mar=c(5,5,5,5),cex=20)

g <- ggbiplot(pcaData, choices = 1:2, obs.scale = 1, var.scale = 1,
              varname.size=0, var.axes = F, ellipse = TRUE, circle = TRUE)+
			  scale_color_manual(name="Sample", values=sequential[2:6]) +  
			  geom_point(aes(colour=colDataRPF$Libs),shape=15,size = 3)+
              theme(legend.direction = 'horizontal',legend.position = c(0.5, 1.13), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
			  scale_y_continuous(limits = c(-15,15))+
			  scale_x_continuous(limits = c(-50,30))
print(g)
dev.off()