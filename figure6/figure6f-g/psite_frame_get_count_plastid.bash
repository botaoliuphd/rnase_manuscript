module load python/2.7.9
module load python/2.7.9_packages/numpy/1.9.2
module load python/2.7.9_packages/pysam/0.8.4
module load python/2.7.9_packages/cython/0.23.4

# pip install --user plastid
# Add plastid path to the $PATH variable

# Generate region of interest. 200nt downstream of annotated start codon for mm10 ucsc gtf file.
metagene generate mm10_cds_start_200 --landmark cds_start \
                                     --annotation_files ucsc.gtf \
                                     --downstream 200

# Generate region of interest. 200nt upstream of annotated stop codon for mm10 ucsc gtf file.									 
metagene generate mm10_cds_stop_200 --landmark cds_stop \
                                     --annotation_files ucsc.gtf \
                                     --upstream 200


# Estimate the P-site offset based on the annotated start codons.
									 
psite mm10_cds_start_200_rois.txt rnase1 --min_length 25 --max_length 35 --require_upstream --count_files rnase1.sorted.wodup.uniq.bam --aggregate &
						 
# offset for RNase1 is 12 for 33mer peak.

psite mm10_cds_start_200_rois.txt rnase2 --min_length 25 --max_length 35 --require_upstream --count_files rnase2.sorted.wodup.uniq.bam --aggregate &
							 
# offset for RNase2 is 12 for 33mer peak.
							 
psite mm10_cds_start_200_rois.txt rnase3 --min_length 25 --max_length 35 --require_upstream --count_files rnase3.sorted.wodup.uniq.bam --aggregate &

# offset for RNase3 is 12 for 31mer peak.
							 
psite mm10_cds_start_200_rois.txt rnase4 --min_length 25 --max_length 35 --require_upstream --count_files rnase4.sorted.wodup.uniq.bam --aggregate &
							 
# offset for RNase4 is 11 for 29mer peak.
							 
psite mm10_cds_start_200_rois.txt rnase5 --min_length 25 --max_length 35 --require_upstream --count_files rnase5.sorted.wodup.uniq.bam --aggregate &
							 
# offset for RNase5 is 11 for 28mer peak.

# Manually correct the incorrectly assigned offset. The default should be 12 or 11 depending on the RNase concentration.
							 
# Estimate the frame preference.
# variable offset after manual correction is better than fixed offset.
# Refine the offset to maximize the percentage of reads on frame 1.

phase_by_size mm10_cds_start_200_rois.txt rnase1_os_var --count_files rnase1.sorted.wodup.uniq.bam --fiveprime_variable --offset rnase1_p_offsets.txt --codon_buffer 5 --min_length 25 --max_length 35 &

phase_by_size mm10_cds_start_200_rois.txt rnase2_os_var --count_files rnase2.sorted.wodup.uniq.bam --fiveprime_variable --offset rnase2_p_offsets.txt --codon_buffer 5 --min_length 25 --max_length 35 &
				
phase_by_size mm10_cds_start_200_rois.txt rnase3_os_var --count_files rnase3.sorted.wodup.uniq.bam --fiveprime_variable --offset rnase3_p_offsets.txt --codon_buffer 5 --min_length 25 --max_length 35 &

phase_by_size mm10_cds_start_200_rois.txt rnase4_os_var --count_files rnase4.sorted.wodup.uniq.bam --fiveprime_variable --offset rnase4_p_offsets.txt --codon_buffer 5 --min_length 25 --max_length 35 &			

phase_by_size mm10_cds_start_200_rois.txt rnase5_os_var --count_files rnase5.sorted.wodup.uniq.bam --fiveprime_variable --offset rnase5_p_offsets.txt --codon_buffer 5 --min_length 25 --max_length 35 &				
						
# Extract the counts at each nucleotide position using P-sites of RPFs with variable offsets.

bsub -o log_get_count -q short -n 1 -W 240 -R rusage[mem=20000] "get_count_vectors --annotation_files mm10_genome.bed --annotation_format BED --count_files rnase1.sorted.wodup.uniq.bam --out_prefix rnase1_ --min_length 25 --max_length 35 --fiveprime_variable --offset rnase1_p_offsets.txt rnase1"

bsub -o log_get_count -q short -n 1 -W 240 -R rusage[mem=20000] "get_count_vectors --annotation_files mm10_genome.bed --annotation_format BED --count_files rnase2.sorted.wodup.uniq.bam --out_prefix rnase2_ --min_length 25 --max_length 35 --fiveprime_variable --offset rnase2_p_offsets.txt rnase2"

bsub -o log_get_count -q short -n 1 -W 240 -R rusage[mem=20000] "get_count_vectors --annotation_files mm10_genome.bed --annotation_format BED --count_files rnase3.sorted.wodup.uniq.bam --out_prefix rnase3_ --min_length 25 --max_length 35 --fiveprime_variable --offset rnase3_p_offsets.txt rnase3"

bsub -o log_get_count -q short -n 1 -W 240 -R rusage[mem=20000] "get_count_vectors --annotation_files mm10_genome.bed --annotation_format BED --count_files rnase4.sorted.wodup.uniq.bam --out_prefix rnase4_ --min_length 25 --max_length 35 --fiveprime_variable --offset rnase4_p_offsets.txt rnase4"

bsub -o log_get_count -q short -n 1 -W 240 -R rusage[mem=20000] "get_count_vectors --annotation_files mm10_genome.bed --annotation_format BED --count_files rnase5.sorted.wodup.uniq.bam --out_prefix rnase5_ --min_length 25 --max_length 35 --fiveprime_variable --offset rnase5_p_offsets.txt rnase5"						