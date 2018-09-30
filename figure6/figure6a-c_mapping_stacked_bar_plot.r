# module load R/3.4.0
#### Mapping stacked bar plots ####

#### Mapping stats ######
library(RColorBrewer)
sequential <- brewer.pal(6, "Blues")

dataset <- read.table("hip_rnase_titration_mapping.txt",sep="\t",header=T)

pdf("hip_rnase_titration_mapping_bar.pdf",width=7,height=6)
par(mar=c(5.1,6.1,4.1,2.1))

barplot(t(dataset[,2:7]),
names.arg = dataset$Sample, # x-axis labels
cex.names = 1, # makes x-axis labels small enough to show all
col = sequential, # colors
xlab = "",
ylab = "Percentage",
xlim = c(0,8), # these two lines allow space for the legend
width = 1,cex.lab=1.5,cex.axis=1.2) # these two lines allow space for the legend
legend("bottomright", 
legend = c("Duplicates","Multimapped", "Unmapped","tRNA", "rRNA","Unique"), #in order from top to bottom
fill = sequential[6:1], # 6:1 reorders so legend order matches graph
bty = "n",cex=1)
dev.off()


#### mRNA region ####
library(RColorBrewer)
sequential <- brewer.pal(5, "Blues")

dataset <- read.table("hip_rnase_titration_region.txt",sep="\t",header=T)

pdf("hip_rnase_titration_region_bar.pdf",width=7,height=6)
par(mar=c(5.1,6.1,4.1,2.1))

barplot(t(dataset[,2:6]),
names.arg = dataset$Sample, # x-axis labels
cex.names = 1, # makes x-axis labels small enough to show all
col = sequential, # colors
xlab = "",
ylab = "Percentage",
xlim = c(0,8), # these two lines allow space for the legend
width = 1,cex.lab=1.5,cex.axis=1.2) # these two lines allow space for the legend
legend("bottomright", 
legend = c("Intergenic", "Intron","3'UTR", "5'UTR","CDS"), #in order from top to bottom
fill = sequential[5:1], # 5:1 reorders so legend order matches graph
bty = "n",cex=1)
dev.off()

#### frame ####
library(RColorBrewer)
sequential <- brewer.pal(3, "Blues")

dataset <- read.table("hip_rnase_titration_frame.txt",sep="\t",header=T)

pdf("hip_rnase_titration_frame_bar.pdf",width=7,height=6)
par(mar=c(5.1,6.1,4.1,2.1))

barplot(t(dataset[,2:4]),
names.arg = dataset$Sample, # x-axis labels
cex.names = 1, # makes x-axis labels small enough to show all
col = sequential, # colors
xlab = "",
ylab = "Percentage",
xlim = c(0,8), # these two lines allow space for the legend
width = 1,cex.lab=1.5,cex.axis=1.2) # these two lines allow space for the legend
legend("bottomright", 
legend = c("Frame3", "Frame2","Frame1"), #in order from top to bottom
fill = sequential[3:1], # 5:1 reorders so legend order matches graph
bty = "n",cex=1)
dev.off()
