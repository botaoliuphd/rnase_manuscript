#!/usr/bin/perl -w
# select_abundant_isoform.pl

use Data::Dumper;
use strict; use warnings;

open (COUNTS,"<combined_rnase_4.isoforms.results");
open (ANNOTATION,"<mm10_refseq_cds_utr_annotation.txt");
open (OUT,">rnase4_isoform_filtered_summary.txt");

my %counts;
my %annotation;

while ( <ANNOTATION> ){
	chomp;
	my @info=split(/\t/,$_);
	$info[0]=~s/\s+//g;
	$info[1]=~s/\s+//g;
	$info[2]=~s/\s+//g;
	$info[3]=~s/\s+//g;
	$info[4]=~s/\s+//g;
	$info[5]=~s/\s+//g;
	$info[6]=~s/\s+//g;
	my @record=($info[1],$info[2],$info[3],$info[4],$info[5],$info[6]);
	if($info[0] ne "transcript_id"){
	$annotation{$info[0]}=\@record;
}
}

while ( <COUNTS> ){
	chomp;
	my @info=split(/\t/,$_);
	$info[0]=~s/\s+//g;
	$info[1]=~s/\s+//g;
	$info[2]=~s/\s+//g;
	$info[4]=~s/\s+//g;
	$info[7]=~s/\s+//g;
	my @record=($info[0],$info[2],$info[4],$info[5],$info[7]);
	if($info[0] ne "transcript_id" && $info[4]+0.0 && $info[5]+0.0 && exists $annotation{$info[0]}){
	if(!exists $counts{$info[1]}){
	$counts{$info[1]}=\@record;
	}
	elsif($info[7] > $counts{$info[1]}->[4]){
	$counts{$info[1]}=\@record;
	}
	elsif($info[7] == $counts{$info[1]}->[4] && $annotation{$info[0]}->[3] > $annotation{$counts{$info[1]}->[0]}->[3]){
	$counts{$info[1]}=\@record;
	}
	elsif($info[7] == $counts{$info[1]}->[4] && $annotation{$info[0]}->[3] == $annotation{$counts{$info[1]}->[0]}->[3] && $annotation{$info[0]}->[4] > $annotation{$counts{$info[1]}->[0]}->[4]){
	$counts{$info[1]}=\@record;
	}
	elsif($info[7] == $counts{$info[1]}->[4] && $annotation{$info[0]}->[3] == $annotation{$counts{$info[1]}->[0]}->[3] && $annotation{$info[0]}->[4] == $annotation{$counts{$info[1]}->[0]}->[4] && $annotation{$info[0]}->[5] > $annotation{$counts{$info[1]}->[0]}->[5]){
	$counts{$info[1]}=\@record;
	}
}
}

print OUT "transcript_id\tsymbol\tcds_start\tcds_stop\ttranscript_length\tcds_length\tlength_5utr\tlength_3utr\tisopct\texpected_counts\ttpm\n";

foreach my $key (keys %counts){
print OUT "$counts{$key}->[0]\t$key\t$annotation{$counts{$key}->[0]}->[0]\t$annotation{$counts{$key}->[0]}->[1]\t$annotation{$counts{$key}->[0]}->[2]\t$annotation{$counts{$key}->[0]}->[3]\t$annotation{$counts{$key}->[0]}->[4]\t$annotation{$counts{$key}->[0]}->[5]\t$counts{$key}->[4]\t$counts{$key}->[2]\t$counts{$key}->[3]\n";
}

close COUNTS;
close ANNOTATION;
close OUT;