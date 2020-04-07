#!/usr/bin/perl
use strict;
use warnings;

# this script is used to generate a pair of matched input files for obaDIA pipeline
# if you want to use it, you should modify it if needed

die "perl $0 <ab> <seq>" unless @ARGV==2;

# read raw abundance file
open IN, "$ARGV[0]" or die $!;
my %hash;
my $head=<IN>;

while(<IN>){
	chomp;
	my @tmp = split /\t/;
	$hash{$tmp[0]}=$_;;
}
close IN;

open OUT1, ">$ARGV[0].tsv" or die $!;
print OUT1 "$head";

open OUT2, ">$ARGV[1].fa" or die $!;

# read raw sequence file
open IN, "$ARGV[1]" or die $!;
$/ = '>';
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\n/;

	# seperater is depend on your fasta format!
	# my @name = split /\s+/, $tmp[0];
	my @name = split /\|/, $tmp[0];
	$name[0] =~ s/\>//;

	# id is depend on your fasta fromat !
	# my $id = $name[0];
	my $id = $name[1];
	if(exists $hash{$id}){
		shift(@tmp);
		my $seq = join('',@tmp);
		print OUT1 "$hash{$id}\n";
		print OUT2 ">$id\n$seq\n";
	}
}
$/ = "\n";

close IN;
close OUT1;
close OUT2;
