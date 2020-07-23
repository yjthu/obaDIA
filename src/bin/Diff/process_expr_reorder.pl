#!/usr/bin/perl
use strict;
use warnings;


#A2/A4/A9,P10/P2/P6,A2H/A4H/A9H,P1H/P3H/P10H

die "perl $0 <table1> <output_prefix> <order>" unless @ARGV>1;

my %order;
my $order=$ARGV[2];
$order=~s/\//,/g;
my @order=split /,/, $order;
my %hash;


open IN, "$ARGV[0]" or die $!;
open OUT, ">$ARGV[1].order.xls" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;

print OUT "$head\n";

while(<IN>){
        chomp;
        my @tmp=split /\t/;
        for(my $i=1; $i<=$#tmp; $i++){
			unless($tmp[$i]){
				$tmp[$i]=0;
			}

			my $id=$head[$i];
			$hash{$id}=$tmp[$i];
        }

        print OUT "$tmp[0]";
        foreach my $e(@order){
			print OUT "\t$hash{$e}";
        }

        print OUT "\n";

}
close IN;




