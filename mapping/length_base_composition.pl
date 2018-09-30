#!/usr/bin/perl -w

use Data::Dumper;
use strict;
use Getopt::Long;

my $inputfile = shift @ARGV;
my $outprefix = shift @ARGV;

# A, T ,C ,G

my %hash_length;
my %hash_base;

open IN, $inputfile or die "Can't open input file $inputfile\n";

while (<IN>) {
		chomp;
		my $read_len = length($_);
		if(exists $hash_length{$read_len}) {
			$hash_length{$read_len}++;}
		else{
			$hash_length{$read_len}=1;
		}
		my @data=split(//,$_);
		for my $i (0 .. $#data)
		{
	
		if(exists $hash_base{$i}) {
			if($data[$i]eq'A'){
			$hash_base{$i}->[0]++;
			}
			elsif($data[$i]eq'T'){
			$hash_base{$i}->[1]++;
			}
			elsif($data[$i]eq'C'){
			$hash_base{$i}->[2]++;
			}
			elsif($data[$i]eq'G'){
			$hash_base{$i}->[3]++;
			}
			}
		else{
			$hash_base{$i}=[0,0,0,0];
		}	
		
	}

}

open OUT_LEN, ">$outprefix.length";
open OUT_BASE, ">$outprefix.base";
print OUT_BASE "position\tA\tT\tC\tG\n";

for my $key ( sort {$a<=>$b} keys %hash_length) {
           print OUT_LEN "$key\t$hash_length{$key}\n";
}

for my $key ( sort {$a<=>$b} keys %hash_base) {
			my $position=$key+1;
		   print OUT_BASE "$position\t$hash_base{$key}->[0]\t$hash_base{$key}->[1]\t$hash_base{$key}->[2]\t$hash_base{$key}->[3]\n";
}

close IN;
close OUT_LEN;
close OUT_BASE;