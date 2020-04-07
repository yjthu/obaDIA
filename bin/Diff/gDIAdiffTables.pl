#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <mapDIAout> <GROUPname> <COM> <OUTDIR> <FC_cut> <FDR_cut>" unless @ARGV==6;
my $groupname = $ARGV[1];
my $compare = $ARGV[2];
my $outdir = $ARGV[3];
my $fc_cutoff = $ARGV[4];
my $fdr_cutoff = $ARGV[5];

my @groupname = split /,/, $groupname;
my @com;

my @comp = split /,/, $compare;
for(my $i=0; $i<=$#comp; $i++){
	my $comp1 = int((split /:/, $comp[$i])[0]) - 1;
	my $comp2 = int((split /:/, $comp[$i])[1]) - 1;
	push @com, $groupname[$comp1]."vs".$groupname[$comp2];
}

foreach my $e(@com){
	`mkdir -p $outdir/$e`;
}

open IN, "$ARGV[0]" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;
for(my $i=1; $i<=$#head; $i++){
	$head[$i] =~ s/\//vs/;
}

my $head2=<IN>;
chomp($head2);
my @head2=split /\t/, $head2;

my %hash;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $pid=$tmp[0];
	for(my $i=1; $i<=$#tmp; $i++){
		my $grp = $head[$i];
		my $headn = $head2[$i];
		if($tmp[$i] ne ''){
			$hash{$grp}{$pid}{$headn}=$tmp[$i];
		}
	}
}
close IN;

foreach my $k(sort keys %hash){
	open OUT, ">$outdir/$k/$k.diff.xls" or die $!;
	print OUT "Protein\tlog2FC\tscore\tFDR\tsignificant\n";
	foreach my $k2(sort keys %{$hash{$k}}){
		my $sig="FALSE";
		if(abs($hash{$k}{$k2}{'log2FC'}) >= $fc_cutoff and $hash{$k}{$k2}{'FDR'} <= $fdr_cutoff){
			$sig="TRUE";
		}
		print OUT "$k2\t$hash{$k}{$k2}{'log2FC'}\t$hash{$k}{$k2}{'score'}\t$hash{$k}{$k2}{'FDR'}\t$sig\n";
	}
	close OUT;
}

		
	

