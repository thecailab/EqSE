#!/usr/bin/perl

use strict;
use warnings;

open IN, "samtools view -h $ARGV[0]|" or die $!;

my $outfile4 = "H3K4me1_".$ARGV[0];
my $outfile27 = "H3K27ac_".$ARGV[0];

open OUT4, "|samtools view - -b > $outfile4" or die $!;
open OUT27, "|samtools view - -b > $outfile27" or die $!;

while(<IN>){
	chomp;
	if(m/^\@/){
		print OUT4 $_."\n";
		print OUT27 $_."\n";
	}
	else{
		my @sp = split/\s+/, $_;
		my @sp1 = split/:/, $sp[0];
		if($sp1[3]<=6){
			print OUT4 $_."\n";
		}
		else{
			print OUT27 $_."\n";
		}
	}
}

close IN;
close OUT4;
close OUT27;