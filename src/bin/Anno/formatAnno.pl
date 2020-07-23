#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
<IN>;
print "id\tPfam\tSignalP\teggNOG\tKEGG\tGO\n";
while(<IN>){
	chomp;
	
	my @tmp=split /\t/;
	my $id=$tmp[0];
	my $pfam='.';
	if($tmp[1] ne '.'){
		my @temp=split /\`/, $tmp[1];
		my $anno;
		foreach my $e(@temp){
			my @each=split /\^/, $e;
			$anno=$each[0].'|'.$each[1].'|'.$each[2];
			$pfam .= $anno.'; ';
		}
		$pfam =~ s/; $//g;
        	$pfam =~ s/^\.//g;
	}
	my $sigP=$tmp[2];
	$sigP=~s/\^/\|/g;
	$sigP=~s/\`/; /g;
	my $eggN=$tmp[3];
	$eggN=~s/\^/\|/g;
        $eggN=~s/\`/; /g;
	my $kegg=$tmp[4];
	$kegg=~s/\^/\|/g;
        $kegg=~s/\`/; /g;
	my $go='.';
	if($tmp[5].'`'.$tmp[6] ne '.`.'){
		my @temp=split /\`/, $tmp[5].'`'.$tmp[6];
		my %hash=();
		foreach my $e(@temp){
			$hash{$e}='';
		}

		foreach my $k(sort keys %hash){
			$k=~s/\^/\|/g;
			$go .= $k.'; ';
		}
		$go=~ s/; $//g;
		$go=~ s/^\.//g;
		$go=~ s/\.;//g;
		$go=~ s/biological_process/BP/g;
		$go=~ s/cellular_component/CC/g;
		$go=~ s/molecular_function/MF/g;	
	}
	print "$id\t$pfam\t$sigP\t$eggN\t$kegg\t$go\n";
}
close IN;	
	
	
	
