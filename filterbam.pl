#!/usr/bin/perl

use strict;
use warnings;

my $prefix = $ARGV[1];
$prefix =~ s/_.+//;

my $outfile = $prefix."_".$ARGV[0];

open OUT, "|samtools view - -b > $outfile" or die $!;

open IN, "<$ARGV[1]" or die $!;
my %keep;
while(<IN>){
        chomp;
        $keep{$_}=1;
}
close IN;

open IN, "samtools view -h $ARGV[0]|" or die $!;
while(<IN>){
        chomp;
        if(m/^\@/){
                print OUT $_."\n";
        }
        else{
                my @sp = split/\s+/, $_;
                my @sp1=split/:/, $sp[0];
                my $key=$sp1[1].':'.$sp1[2].':'.$sp1[3];
                if(exists($keep{$key})){
                        print OUT $_."\n";
                }
        }
}
close IN;

close OUT;