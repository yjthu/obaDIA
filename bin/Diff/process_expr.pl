#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <table1> <output_prefix> <order>" unless @ARGV>1;

#A2/A4/A9,P10/P2/P6,A2H/A4H/A9H,P1H/P3H/P10H
my %order;
my $order=$ARGV[2];
$order=~s/\//,/g;
my @order=split /,/, $order;
my %hash;
my %sum;
my $total=0;
my $num=@order;

foreach my $each(@order){
	@{$hash{$each}}=();
	$sum{$each}=0;
}

open IN, "$ARGV[0]" or die $!;
open OUT, ">$ARGV[1].normalized.xls" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;
my @prot;

while(<IN>){
	chomp;
	my @tmp=split /\t/;
	push @prot, $tmp[0];
	for(my $i=1; $i<=$#tmp; $i++){
		my $k=$head[$i];
		push @{$hash{$k}}, $tmp[$i];
		$sum{$k}+=$tmp[$i];
		$total+=$tmp[$i];
	 }	
}
close IN;

print OUT "ID";
foreach my $k(@order){
        print OUT "\t$k";
}
print OUT "\n";

for(my $i=0; $i<=$#prot; $i++){
	my $id = $prot[$i];
	print OUT "$id";
	foreach my $k(@order){
		my @arr = @{$hash{$k}};
		my $val = $arr[$i]*$total/$sum{$k}/$num;
		print OUT "\t$val";
	}
	print OUT "\n";
}




