#!/usr/bin/perl


if(@ARGV<2) { &help; exit; }

my $inputfile = shift @ARGV;
my $outfile = shift @ARGV;

open OUT, ">$outfile";
open IN, $inputfile or die "Can't open input file $inputfile\n";
while($line1=<IN>){
 chomp($line1);
 $line1=~s/\s+/:/;
 chomp($line2=<IN>);
 chomp($line3=<IN>);
 chomp($line4=<IN>);

 print OUT "$line1\n$line2\n$line3\n$line4\n";
}
close IN;
close OUT;
