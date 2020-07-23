#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;
my %hash;
for(my $i=1; $i<=$#head; $i++){
	my $id=$head[$i];
	$hash{$id}=0;
}

	
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	for(my $i=1; $i<=$#tmp; $i++){
		my $id=$head[$i];
		if($tmp[$i] > 0){
			$hash{$id}++;
		}	
	}
}
close IN;

print "Sample\tCount\n";
for(my $i=1; $i<=$#head; $i++){
	my $id=$head[$i];
	print "$id\t$hash{$id}\n";
}

