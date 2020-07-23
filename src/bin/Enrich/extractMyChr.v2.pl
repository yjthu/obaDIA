#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <input> <type> <list>" unless @ARGV==3;
open IN1, "$ARGV[2]" or die $!;
my %hash;
while(<IN1>){
	chomp;
	my @tmp = split;
	$hash{$tmp[0]}='';
}

if($ARGV[1] eq 'fa'){
open IN, "$ARGV[0]" or die $!;
$/ = '>';
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\n/;
	my @name = split /\s+/, $tmp[0];
	$name[0] =~ s/\>//;
	
	if(exists $hash{$name[0]}){
		shift(@tmp);
		my $seq = join('',@tmp);
		print ">$name[0]\n$seq\n";
	}
}
$/ = "\n";
close IN;
}

if($ARGV[1] eq 'gtf'){
open IN, "$ARGV[0]" or die $!;
while(<IN>){
	chomp;
	my @tmp=split;
	if(exists $hash{$tmp[0]}){
		print "$_\n";
	}
}
close IN;
}
