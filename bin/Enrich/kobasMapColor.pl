#!/usr/bin/perl
use strict;
use warnings;


die "perl $0 <anno> <enrich>" unless @ARGV==2;
open IN, "$ARGV[0]" or die $!;
open IN2, "$ARGV[1]" or die $!;
my %anno;
while(<IN>){
	chomp;
	next if(/^#/);
	next if(/^$/);
	next if(/\---/);
	last if(/\/\/\//);
	my @tmp = split /\s+/;
	if($tmp[1]){
		my @tmp1 = split /\|/,$tmp[1];
		$anno{$tmp[0]}=$tmp1[0];
	}
}
close IN;


my $sharp=0;
while(<IN2>){
	chomp;
	if(/^#/){
		$sharp++;
		if($sharp < 5){
			print "$_\n";
		}
		next;
	}
	next if(/^$/);
	next if(/^---/);
	my @tmp=split /\t/;
	$tmp[-1] =~ s/red/pink/g;
	
	for(my $i=0; $i<$#tmp; $i++){
		print "$tmp[$i]\t";
	}
	my $h=$tmp[-1];
	if(length($tmp[-1]) > 2000){
			$tmp[-1] = substr($tmp[-1], 0, 2000);
			my @temp=split '/', $tmp[-1];
			pop @temp;
			$h = join('/', @temp);
	}

	print "$h\n";


	
}
close IN2;
				
