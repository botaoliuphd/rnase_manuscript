#!/usr/bin/perl -w

use Data::Dumper;
use strict; use warnings;

open (RSEM,"<wt_rna_batch1_rep1.isoforms.results");
open (ANNOTATION,"<mm10_refseq_annotation_seq.txt");
open (OUT,">wt_rna_batch1_rep1_most_abundant_isoform_annotation.txt");

my %annotation;

while ( <ANNOTATION> ){
	chomp;
	my @info=split(/\t/,$_);
	for(my $i=0; $i < @info; $i++) {
    $info[$i]=~s/\s+//g;
	}
	if($info[0] ne "transcript_id" && scalar @info == 15){
	$annotation{$info[0]}=\@info;
}
}

# Select the most abudant transcript isoform based on the RSEM estimation.
my %gc;
while ( <RSEM> ){
	chomp;
	my @info=split(/\t/,$_);
	for(my $i=0; $i < @info; $i++) {
    $info[$i]=~s/\s+//g;
	}
if($info[0] ne "transcript_id"){	
	my @record=($info[0],$info[2],$info[4],$info[5],$info[7]);
	if($info[0] ne "transcript_id" && $info[4]+0.0 && $info[5]+0.0 && exists $annotation{$info[0]}){
	if(!exists $gc{$info[1]}){
	$gc{$info[1]}=\@record;
	}
	elsif($info[7] > $gc{$info[1]}->[4]){
	$gc{$info[1]}=\@record;
	}
	elsif($info[7] == $gc{$info[1]}->[4] && $annotation{$info[0]}->[4] > $annotation{$gc{$info[1]}->[0]}->[4]){
	$gc{$info[1]}=\@record;
	}
	elsif($info[7] == $gc{$info[1]}->[4] && $annotation{$info[0]}->[4] == $annotation{$gc{$info[1]}->[0]}->[4] && $annotation{$info[0]}->[5] > $annotation{$gc{$info[1]}->[0]}->[5]){
	$gc{$info[1]}=\@record;
	}
	elsif($info[7] == $gc{$info[1]}->[4] && $annotation{$info[0]}->[4] == $annotation{$gc{$info[1]}->[0]}->[4] && $annotation{$info[0]}->[5] == $annotation{$gc{$info[1]}->[0]}->[5] && $annotation{$info[0]}->[6] > $annotation{$gc{$info[1]}->[0]}->[6]){
	$gc{$info[1]}=\@record;
	}
}
}
}

print OUT "symbol\ttranscript_id\tisopct\texpected_gc\ttpm\tcds_start\tcds_stop\tlength_transcript\tlength_cds\tlength_5utr\tlength_3utr\tgc_transcript\tgc_cds\tgc_5utr\tgc_3utr\n";

foreach my $key (keys %gc){

print OUT "$key\t$gc{$key}->[0]\t$gc{$key}->[4]\t$gc{$key}->[2]\t$gc{$key}->[3]\t$annotation{$gc{$key}->[0]}->[1]\t$annotation{$gc{$key}->[0]}->[2]\t$annotation{$gc{$key}->[0]}->[3]\t$annotation{$gc{$key}->[0]}->[4]\t$annotation{$gc{$key}->[0]}->[5]\t$annotation{$gc{$key}->[0]}->[6]\t$annotation{$gc{$key}->[0]}->[7]\t$annotation{$gc{$key}->[0]}->[8]\t$annotation{$gc{$key}->[0]}->[9]\t$annotation{$gc{$key}->[0]}->[10]\n";
}

close RSEM;
close ANNOTATION;
close OUT;