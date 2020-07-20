#!/usr/bin/perl
# Script to convert WhatHap polyphase (Schrinner et al. 2020) output (in VCF format) to haplogenerator format, which is accepted by Hapcompare.py
# As an example, one uses ./WhatsHap2Haplogenerator.pl test.vcf 4, where test.vcf is the input WhatsHap phasing
# and 4 is the ploidy level. The output will be save in test.vcf.haplogen.
# Written by Ehsan Motazedi, 20-07-2020
 
use warnings;
use strict;

open(my $FILE, "<", $ARGV[0]) || die "cannot open ${ARGV[0]}!";
open(my $OUT, ">", "${ARGV[0]}.haplogen") || die "cannot write to ${ARGV[0]}.haplogen!";
my $ploidy = $ARGV[1];
my @haps = map { 'hap_'.$_ } 1..$ploidy;
my @header = ('var_id', 'contig', 'varpos', 'ref_allele', 'alt_allele');
push(@header, @haps);
printf $OUT "%s\n", join "\t", @header;

my $count = 0; 
my $hline_char="";
my @tokens=();
my $genos="";
my @alleles=();
my $snp_id="";

foreach my $line (<$FILE>){   
    chop $line;
    $hline_char = substr($line, 0, 1);
	if($hline_char ne "#"){
		@tokens = split /\t/, $line;
		$genos= (split /:/, $tokens[9])[0];
		if ($genos !~ m/\//) {
			$count++;
			@alleles = split /\|/, $genos;
			$genos = join "\t", @alleles;
			$snp_id = "SNP_${count}";
			printf $OUT "%s\t%s\t%s\t%s\t%s\t%s\n", $snp_id, $tokens[0], $tokens[1], $tokens[3], $tokens[4], $genos
		}
	}
}
	
close($FILE);
close($OUT)
