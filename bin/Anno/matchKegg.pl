#!/usr/bin/perl
use strict;
use warnings;

my $db_dir="/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db";
my $mapfile="$db_dir/ko00001.keg";
my %hash;
open IN, "$mapfile" or die $!;
while(<IN>){
	chomp;
	next unless(/^D/);
	if(/^D\s+(K\d+?)\s+.*;(.*?)\[/){
		my $ko=$1;
		my $anno=$2;
		$anno =~ s/^\s+//;
		$anno =~ s/\s+$//;
		$anno =~ s/\W+//;
		if(length($anno) > 100){
			$anno = substr($anno, 0, 100);
		}
#		print "$ko\t$anno\n";
		unless($anno){
			$anno=$ko;
		}		
		$hash{$ko}=$anno;
	}
}
close IN;

open IN, "$ARGV[0]" or die $!;
<IN>;
while(<IN>){
	chomp;
#	next if(/^#/);
	my @tmp=split /\t/;
	next unless ($tmp[1] =~ /KO:/);
	my $id=(split /KO:/, $tmp[1])[1];
	if(exists $hash{$id}){
		$hash{$id} =~ s/\s+/_/g;
		print "$tmp[0]\t$id\t$hash{$id}\n";
	}
}
close IN;
