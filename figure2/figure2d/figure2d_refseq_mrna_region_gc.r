# Transcript annotation and GC content are extracted with plastid package from mm10 fastq sequences and bed annotation file download from UCSC genome browser. 
# python extract_feature_anonotation_seq_genome.py &

# The most abundant transcript isoform is selected based on the RSEM estimation to represent each gene in the genome.
# perl extract_gc_isoform.pl &
# The output file is used for plotting with the following R scripts. 

#module load R/3.4.0
#R
## Box plot of gc content in differnt mRNA regions ##

library(genefilter)
require(reshape2)
require(ggplot2)
library(ggsignif)
library(ggpubr)

all_data <-read.table("wt_rna_batch1_rep1_most_abundant_isoform_annotation.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = TRUE)

all_data[,c("Full")] <- all_data[,c("gc_transcript")]
all_data[,c("CDS")] <- all_data[,c("gc_cds")]
all_data[,c("5'UTR")] <- all_data[,c("gc_5utr")]
all_data[,c("3'UTR")] <- all_data[,c("gc_3utr")]

all_box <- all_data[,c("Full","5'UTR","CDS","3'UTR")]
all_melt <- melt (all_box)
combined_melt <- all_melt


stats <- boxplot.stats(combined_melt[combined_melt[,"variable"]=="CDS","value"])$stats
df1 <- data.frame(x="label1", ymin=stats[1], lower=stats[2], middle=stats[3], 
                 upper=stats[4], ymax=stats[5])
				 
stats <- boxplot.stats(combined_melt[combined_melt[,"variable"]=="3'UTR","value"])$stats
df2 <- data.frame(x="label1", ymin=stats[1], lower=stats[2], middle=stats[3], 
                 upper=stats[4], ymax=stats[5])

ymax <- 1.1*df1$ymax				 
ymin <- 0.9*df2$ymin

			  
p <- ggplot(data = combined_melt, aes(x=variable, y=value)) 
p <- p + theme_classic()
p <- p + theme(plot.margin = unit(c(1,1,1,1), "cm"),strip.background = element_blank(),strip.text.x = element_text(colour="black",size=20),panel.spacing = unit(2, "lines"),axis.text.x = element_text(colour="black",size=20),axis.text.y = element_text(colour="black",size=20),axis.title.y = element_text(colour="black",size=20),axis.title.x = element_blank(),legend.text=element_blank(),legend.title=element_blank())
p <- p + coord_cartesian(ylim=c(ymin,ymax))
p <- p + geom_boxplot(aes(fill=variable),position = position_dodge(width = .9), outlier.shape = NA)
p <- p + scale_fill_manual(values=c("white","white","white","white"))
p <- p + stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Full", label.y = c(0.95*ymax))
p <- p + ylab("GC%")
p <- p + guides(fill=FALSE)

output_name="wt_rna_batch1_rep1.mm10.gc.boxplots.pdf"
pdf(output_name,width=5,height=6)
p
dev.off()
