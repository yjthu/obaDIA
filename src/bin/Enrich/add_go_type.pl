#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my %type;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	$type{$tmp[0]}=$tmp[2];
}
close IN;

print "#Term\tDatabase\tID\tTerm_type\tInput_number\tBackground_number\tP_Value\tCorrected_P_Value\tInput\tHyperlink\n";
open IN, "$ARGV[1]" or die $!;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $type='na';
	if(exists $type{$tmp[2]}){
                $type=$type{$tmp[2]};
        }
	next if($type eq 'na');
	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$type";
	$tmp[-1] =~ s/.*GO://;
	$tmp[-1] = 'http://amigo.geneontology.org/amigo/search/ontology?q=GO:'.$tmp[-1];
	for(my $i=3;$i<=$#tmp; $i++){
		print "\t$tmp[$i]";
	}

	print "\n";
}
close IN;
	
	

