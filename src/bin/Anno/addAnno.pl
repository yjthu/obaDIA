#!/usr/bin/perl 
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my %nog;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	$nog{$tmp[0]}=$tmp[2];
}
close IN;

open IN, "$ARGV[1]" or die $!;
my %kegg;
while(<IN>){
        chomp;
        my @tmp=split /\t/;
        $kegg{$tmp[0]}=$tmp[2];
}
close IN;

open IN, "$ARGV[2]" or die $!;
my $head=<IN>;
print "$head";
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	if(exists $nog{$tmp[0]}){
		$tmp[3] .= '|'.$nog{$tmp[0]};
	}
	if(exists $kegg{$tmp[0]}){
                $tmp[4] .= '|'.$kegg{$tmp[0]};
        }
	for(my $i=0; $i<$#tmp; $i++){
		print "$tmp[$i]\t";
	}
	print "$tmp[$#tmp]\n";

}
close IN;

