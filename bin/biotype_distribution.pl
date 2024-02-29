#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
##########################################################################
# Input: 
#   arg0: file containing the list of biotypes to report on, one per line.
#   arg1: 2-column TSV file listing unique pairs of fragment ID and biotype.
#         Should be the output of "sort|uniq". 
#          Column #1: fragment ID; 
#          Column #2: gene_type (aka biotype). Example:
#               id19718_TCGCGCAA     misc_RNA
#               id11168_GAACCTCT     lncRNA
#               id11168_GAACCTCT     misc_RNA
# 
# Output:
#   Tab-separated output, one line per biotype.
#     Column #1: biotype (as listed in arg0)
#     Column #2: number of fragments in arg1 that overlap biotype. Count can 
#                be fractional, when a single fragment overlap several distinct 
#                biotypes
#
##########################################################################

open BIOTYPES, "$ARGV[0]" or die $!; #list of biotypes to include in report, one per line. This needs to be provided to guarantee that all biotypes will be present in the output, even those that are not covered.

# Initialize %biotype_counts with zeroes
my %biotype_counts=();
while (<BIOTYPES>){
    chomp;
    my $biotype=$_;
    $biotype_counts{$biotype}=0;
}
close BIOTYPES;
open INTERSECT, "$ARGV[1]" or die $!; #list of intersects between fragments and biotypes. Tab-separated, one unique intersect pair per line (ie should be "sort|uniq"'ed beforehand). 

#initialize previous fragment to random string value. This is just to process the first line of the file.
my $previous_fragment='null'; 
my %fragment_to_biotype=(); #this associates a fragment ID to a list of biotypes. This hash is re-initialized each time a new fragment ID is encountered in the file (that's why arg1 needs to be sorted), in order to save memory

while(<INTERSECT>){
    chomp;
    my @line=split "\t";
    my $current_fragment=$line[0];
    my $biotype=$line[1];
    if($current_fragment ne $previous_fragment){ # this means we switched to next fragment
        updateBiotypeCount($previous_fragment);
        $previous_fragment=$current_fragment;
        %fragment_to_biotype=(); #empty hash to save memory
    }
    push(@{$fragment_to_biotype{$current_fragment}}, $biotype);
    if(eof(INTERSECT)){ #last line reached, we need to write last record
        updateBiotypeCount($current_fragment);
    }
}
#print tab-separated output
foreach my $bt (keys %biotype_counts){
    print "$bt\t$biotype_counts{$bt}\n"
}

sub updateBiotypeCount{
    my $fragment=$_[0];
    my $fragment_to_nr_biotypes=$#{$fragment_to_biotype{$fragment}}+1; #count of biotypes for fragment is last index of array +1
    foreach my $bt (@{$fragment_to_biotype{$fragment}}){
        $biotype_counts{$bt}+=1/$fragment_to_nr_biotypes;
    }
}
