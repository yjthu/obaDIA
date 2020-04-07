#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <mapDIAout> <expression> <GROUPname> <COM> <OUTDIR> <FC_cut> <FDR_cut>" unless @ARGV==7;
my $exp = $ARGV[1];
my $groupname = $ARGV[2];
my $compare = $ARGV[3];
my $outdir = $ARGV[4];
my $fc_cutoff = $ARGV[5];
my $fdr_cutoff = $ARGV[6];

open IN, "$exp" or die $!;
my $header=<IN>;
chomp($header);
my @header=split /\t/, $header;
my %exp;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	for(my $i=1; $i<=$#tmp; $i++){
		my $sn=$header[$i];
		$exp{$tmp[0]}{$sn}=$tmp[$i];
	}
}
close IN;

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
		else{
			$hash{$grp}{$pid}{$headn}='NA';
		}
	}
}
close IN;

foreach my $k(sort keys %hash){
	open OUT, ">$outdir/$k/$k.diff.xls" or die $!;
	my @temp = split /vs/, $k;
	my $s1 = $temp[0];
	my $s2 = $temp[1];
	print OUT "Protein\t$s1\t$s2\tlog2FC\tscore\tFDR\tsignificant\n";
	foreach my $k2(sort keys %{$hash{$k}}){
		my $sig="FALSE";
		my $fc=999;
		if($hash{$k}{$k2}{'FDR'} eq 'NA'){
			$hash{$k}{$k2}{'score'} = 0;
			$hash{$k}{$k2}{'FDR'} = 1;

		}
		if($exp{$k2}{$s2} > 0){
			$fc=&log2($exp{$k2}{$s1}/$exp{$k2}{$s2});
		}
		if(abs($fc) >= $fc_cutoff and $hash{$k}{$k2}{'FDR'} <= $fdr_cutoff){
			$sig="TRUE";
		}
		if($exp{$k2}{$s1}==0 and $exp{$k2}{$s2}==0){
			$fc=0;
		}
		print OUT "$k2\t$exp{$k2}{$s1}\t$exp{$k2}{$s2}\t$fc\t$hash{$k}{$k2}{'score'}\t$hash{$k}{$k2}{'FDR'}\t$sig\n";
	}
	close OUT;
}

sub log2 { my $n = shift; return log($n)/log(2); }