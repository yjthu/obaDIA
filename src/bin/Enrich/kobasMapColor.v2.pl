#!/usr/bin/perl
use strict;
use warnings;


die "perl $0 <anno> <diff> <enrich>" unless @ARGV==3;
open IN, "$ARGV[0]" or die $!;
open IN2, "$ARGV[1]" or die $!;
open IN3, "$ARGV[2]" or die $!;
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


my %color;
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	if(exists $anno{$tmp[0]}){
		my $ko=$anno{$tmp[0]};
		if(not exists $color{$ko}){
			if($tmp[1] eq 'up'){
				$color{$ko}='red';
			}
			else{
				$color{$ko}='cyan';
			}
		}
		else{
			if($tmp[1] eq 'up' and $color{$ko} eq 'cyan'){
				 $color{$ko}='yellow';
			}
			if($tmp[1] eq 'down' and $color{$ko} eq 'red'){
                                 $color{$ko}='yellow';
                        }
		}
	}
}
close IN2;
my $sharp=0;
while(<IN3>){
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
	my $hyperlink = $tmp[-1];
	$hyperlink =~ s/.*\?//;
	$hyperlink =~ s/red/green/g;
	my @arr=split /\//,$hyperlink;
	shift @arr;
	foreach my $term(@arr){
		my @ko = split /\%09/, $term;
		if(exists $color{$ko[0]}){
			my $str = $ko[0].'%09'.$color{$ko[0]};
			$tmp[-1] =~ s/$term/$str/;
		}
	}
	for(my $i=0; $i<$#tmp; $i++){
		print "$tmp[$i]\t";
	}
	print "$tmp[-1]\n";
}
close IN3;
				
