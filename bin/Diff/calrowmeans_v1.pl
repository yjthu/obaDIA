#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($rpkm,$group,$groupname,$out);
GetOptions(
	"rpkm=s"=>\$rpkm,
	"group=s"=>\$group,
	"groupname=s"=>\$groupname,
	"out-rowmeans=s"=>\$out,
);

if(!defined($rpkm) ||  !defined($group) ||  !defined($groupname) ||  !defined($out)){
	die "Usage:perl $0 <-rpkm> <-group> <-groupname> <-out-rowmeans>";
	exit 0;
}

open RPKM, "$rpkm" or die $!;

my $line=<RPKM>;
chomp($line);
my @name = split /\t/, $line;
my %rp;
while(<RPKM>){
        chomp;
        my @tmp1 = split /\t/, $_;
        for(my $i=1;$i<@tmp1;$i++) {
        	$rp{$tmp1[0]}{$name[$i]}=$tmp1[$i];
	}
}
close RPKM;

my @gn=split /,/,$groupname;
my @sn=split /,/,$group;

unless(@gn == @sn){
	die "Error input!";
}

my %score;
foreach my $key1 (sort keys %rp) {
	for(my $j=0;$j<@gn;$j++){
		my @tmp2=split /\//,$sn[$j];
		my $number=0;
		for(my $k=0;$k<@tmp2;$k++){
			$score{$key1}{$gn[$j]} += $rp{$key1}{$tmp2[$k]};
			$number++;
		}
		$score{$key1}{$gn[$j]}=$score{$key1}{$gn[$j]}/$number;
	}
}

open OUT,">$out";
#print head
foreach my $keyt1(sort keys %score) {
        print OUT "ID";
        foreach my $keyt2(sort keys %{$score{$keyt1}}) {
                print OUT "\t$keyt2";
        }
        print OUT "\n";
        last;
}
#print OUT
foreach my $Key1(sort keys %score) {
        print OUT "$Key1";
        foreach my $Key2(sort keys %{$score{$Key1}}) {
                print OUT "\t$score{$Key1}{$Key2}";
        }
        print OUT "\n";
}

close OUT;
