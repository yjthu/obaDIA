#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
open IN2, "$ARGV[1]" or die $!;

my %hash;
while(<IN>){
	chomp;
	$hash{$_}='';
}
close IN;

print "##ko    KEGG Orthology\n\n";
while(<IN2>){
	chomp;
	next if(/^#/);
	next if(/^$/);
	last if(/^---/);
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}){
		print "$_\n";
	}
}
close IN2;
