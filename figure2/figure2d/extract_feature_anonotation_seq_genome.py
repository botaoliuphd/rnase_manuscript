#module load python/2.7.9
#module load python/2.7.9_packages/numpy/1.9.2
#module load python/2.7.9_packages/pysam/0.8.4
#module load python/2.7.9_packages/cython/0.23.4
#module load python/2.7.9_packages/biopython/1.66

import plastid
import re
from Bio import SeqIO
from plastid import BAMGenomeArray, FivePrimeMapFactory, BED_Reader, Transcript
from Bio.SeqUtils import GC

# open the mm10 fasta file
genome = SeqIO.to_dict(SeqIO.parse(open("refseq_mm10_ref_rsem_riboprofiling.fa"),"fasta"))

# open the mm10 bed annotation file
transcripts = list(BED_Reader("mm10_genome.bed",return_type=Transcript))

output = open("mm10_refseq_annotation_seq.txt","w")
output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("transcript_id","cds_start","cds_stop","transcript_length","length_cds","length_5utr","length_3utr","gc_transcript","gc_cds","gc_5utr","gc_3utr","seq_transcript","seq_cds","seq_5utr","seq_3utr"))

# extract annotation information for very transcript
for transcript in transcripts:
	utr5 = transcript.get_utr5()
	cds = transcript.get_cds()
	utr3 = transcript.get_utr3()
	
	match = re.search('NM_', transcript.get_name())
	if (match is not None) and transcript.get_name() and transcript.cds_start and transcript.cds_end and utr5 and cds and utr3:
		
		transcript_seq = transcript.get_sequence(genome)
		utr5_seq = utr5.get_sequence(genome)
		cds_seq = cds.get_sequence(genome)
		utr3_seq = utr3.get_sequence(genome)
		transcript_gc = GC(transcript_seq)
		utr5_gc = GC(utr5_seq)
		cds_gc = GC(cds_seq)
		utr3_gc = GC(utr3_seq)
		if transcript_seq and utr5_seq and cds_seq and utr3_seq:
			output.write("%s\t%6i\t%6i\t%6i\t%6i\t%6i\t%6i\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n"%(transcript.get_name(),transcript.cds_start,transcript.cds_end,transcript.length,cds.length,utr5.length,utr3.length,transcript_gc,utr5_gc,cds_gc,utr3_gc,transcript_seq,cds_seq,utr5_seq,utr3_seq))

output.close()
