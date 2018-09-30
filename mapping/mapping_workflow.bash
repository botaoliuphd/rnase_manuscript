# This file exhibits the general reads processing and mapping workflow.
# The pre-processing steps and post-mapping steps are performed in the UMMS HPCC enviroment.
# The mapping steps are performed on the Dolphin platform at the UMMS Bioinformatics Core based on the UMMS HPCC.
# The dependent software and packages are loaded on cluster with the "module load" function.
# Similar analyses could be perfomed in other enviroment after installing all the dependencies.
# Substitue SAMPLE with your run name and SAMPLE_1 etc with individual sample name. Used a human sample as an example. Substitute hg19 with the references for your sample.

##################################
remove the spaces in Nextseq500 read names
###################################
perl remove_space_read_name.pl SAMPLE.fastq SAMPLE.fq &

##################################
#          debarcoding           #
##################################
# make a map_barcodes file similar as the following example to match barcodes with corresponding samples.

GTGAT   SAMPLE_1
CATCG   SAMPLE_2
AGCTA   SAMPLE_3
GGAGC   SAMPLE_4
ACTGT   SAMPLE_5
CGGGA   SAMPLE_6
TAGCC   SAMPLE_7
TGACT   SAMPLE_8
ACAAG   SAMPLE_9
ATTCA   SAMPLE_10
TCTAC   SAMPLE_11

####################################################################
# The splitter package is developed by Garber lab and UMMS Biocore #

module load java/1.8.0_31
bsub -o log_split -q short -W 120 -R rusage[mem=20000] "java -Xmx1g -jar splitter_06.17.15_10.59.jar F1=SAMPLE.fq B=NNNNNNNNSSSSSOO M=map_barcodes HD=2"

## The processing steps above have been done for the fastq files deposited to GEO ##

#################################
#        Remove adaptor         #
#################################
module load cutadapt/1.7.1
bsub -o log_cutadapt -q short -W 240 -R rusage[mem=20000] "cutadapt -a TGGAATTCTCGGGTGCCAAGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -O 2 -q 15 GTGAT_SAMPLE_1.fq -o cutadapt_GTGAT_SAMPLE_1.fq"

##################################################
#                 Map Dolphin                    #
##################################################
#Genome Build   hg19
#Is this library mate-paired?   Single-end
#Is this a fresh run or would you like continue from an old run?        Fresh
#Input Directory(Full path in the cluster)      ### path to your cutadapt output files
#Output Directory(Full path in the cluster)     ### path you want for the output
#Barcode Separation     No
#Input Files (See below for example)    SAMPLE_1 cutadapt_GTGAT_SAMPLE_1.fq
#FastQC reports         Yes
#Adapter Removal        No
#Quality Filtering      Yes
#Window Size    5
#Required Quality       20
#leading        5
#trailing       20
#minlen         25
#Trimming       No
#Common RNAs    Yes
#ERCC   No
#rRNA   Yes
#miRNA  No
#tRNA   Yes
#snRNA  No
#rmsk   No
#Genome         No
#Change Parameters      No
#Split FastQ    No
#Clean Intermediate Files       Yes
#Add Pipelines  Tophat
#Tophat parameters      --segment-length 12
#IGV/TDF Conversion     No
#BigWig Conversion      No

###################################################
##          Mark and Remove Duplicates            #
###################################################
#####################################################################################
# The MarkDuplicatesFromNameUMI package is developed by Garber lab and UMMS Biocore #
module load java/1.8.0_31
bsub -o log_dup -q long -n 5 -W800 -R "select[hname!='ghpcc-sgi'] rusage[mem=2000] span[hosts=1]" "java -Xmx10g -jar MarkDuplicatesFromNameUMI.jar -in SAMPLE_1.sorted.bam -out marked_SAMPLE_1.sorted.bam -bcSize 8 -bcHistogram SAMPLE_1.bcHistogram -bcCountsFile SAMPLE_1.bcCountsFile"

module load samtools/0.0.19
samtools view -F 1024 marked_SAMPLE_1.sorted.bam | awk -F "_|\t" '{if ($2!~/N/) print $0}' | grep -w NH:i:1 | samtools view -bt hg19.fa.fai - > SAMPLE_1.sorted.wodup.uniq.bam &

###################################################
##               Region Analysis                  #
###################################################
# Download the annotatin bed files from UCSC Genome Browser #

bsub -o log -q short -W 200 -R "select[hname!='ghpcc-sgi'] rusage[mem=5000]" "module load bedtools/2.22.0;bedtools intersect -u -f 1.0 -abam SAMPLE_1.sorted.wodup.uniq.bam -b hg19_ucsc.cds.bed > SAMPLE_1_reads_overlap_CDS.bam"
bsub -o log -q short -W 200 -R "select[hname!='ghpcc-sgi'] rusage[mem=5000]" "module load bedtools/2.22.0;bedtools intersect -u -f 1.0 -abam SAMPLE_1.sorted.wodup.uniq.bam -b hg19_ucsc.utr.bed > SAMPLE_1_reads_overlap_UTR.bam"
bsub -o log -q short -W 200 -R "select[hname!='ghpcc-sgi'] rusage[mem=5000]" "module load bedtools/2.22.0;bedtools intersect -u -f 1.0 -abam SAMPLE_1.sorted.wodup.uniq.bam -b hg19_ucsc.intron.bed > SAMPLE_1_reads_overlap_INTRON.bam"
bsub -o log -q short -W 200 -R "select[hname!='ghpcc-sgi'] rusage[mem=5000]" "module load bedtools/2.22.0;bedtools intersect -u -f 1.0 -abam SAMPLE_1.sorted.wodup.uniq.bam -b hg19_ucsc.intergenic.bed > SAMPLE_1_reads_overlap_INTER.bam"

####################################
##base composition footprint length#
####################################
# Extract the read sequences from reads mapped to the CDS region.
samtools view SAMPLE_1_reads_overlap_CDS.bam|awk '{OFS="\t"; print $10}' > SAMPLE_1_reads_overlap_CDS.txt &

# Calculate the base composition at each nucleotide position and the read length distribution.
perl length_base_composition.pl SAMPLE_1_reads_overlap_CDS.txt SAMPLE_1_reads_overlap_CDS &

###################################################
##               RSEM Quantification              #
###################################################
# Convert bam files containing uniquely mapped reads back to fastq files # 
module load bedtools/2.22.0
bsub -o log_bamtofq -q short -W200 -R "select[hname!='ghpcc-sgi'] rusage[mem=10000]" "bamToFastq -i SAMPLE_1.sorted.wodup.uniq.bam -fq SAMPLE_1.sorted.wodup.uniq.fq"

# Map the cleaned reads to Refseq CDS reference with Bowtie2 and quantify the expression with RSEM #
bsub -o log_rsem -n 5 -q long -W 2000 -R "select[hname!='ghpcc-sgi'] rusage[mem=6000] span[hosts=1]" "module load bowtie2/2-2.1.0; module load RSEM/1.2.11; rsem-calculate-expression --bowtie2 -p 5 --output-genome-bam --seed-length 12 --bowtie2-path /share/pkg/bowtie2/2-2.1.0 SAMPLE_1.sorted.wodup.uniq.fq refseq_hg19_ref_rsem_CDS SAMPLE_1"

###################################################
##     Organize count tables for DESEq2           #
###################################################
################################################################################
# The rsem.to.table.perl script is writtern by Alper Kucukural at UMMS Biocore #
perl rsem.to.table.perl -out rsem.gene.summary.tpm.txt -indir . -gene_iso genes -quantType tpm
perl rsem.to.table.perl -out rsem.gene.summary.counts.txt -indir . -gene_iso genes -quantType expected_count
