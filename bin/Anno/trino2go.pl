#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my %hash;
my %go;
<IN>;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $sw='';
	my $go = $tmp[5].'^'.$tmp[6];
	my @go = $go=~/(GO:\d+)/g;
	for my  $each(@go){
		$go{$tmp[0]}{$each}='';
		$go{$tmp[0]}{$each}='';
	}
}
close IN;

foreach my $key1(sort keys %go){
	print "$key1";
	foreach my $key2(sort keys %{$go{$key1}}){
		print "\t$key2";
	}
	print "\n";
}

