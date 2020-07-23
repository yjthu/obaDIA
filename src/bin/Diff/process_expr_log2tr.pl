#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my $head=<IN>;
chomp($head);

print  "$head\n";
while(<IN>){
	chomp;
	my @tmp=split;
	print  "$tmp[0]";
	for(my $i=1; $i<=$#tmp; $i++){
		unless($tmp[$i]){
			$tmp[$i]=0;
		}
		my $num = 2**$tmp[$i];
		print "\t$num";
	}
	print "\n";
}
close IN;

