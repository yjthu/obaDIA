#!/usr/bin/perl
use strict;
use warnings;


die "perl $0 <anno> <diff> <bg enrich> <enrich> <annmap>" unless @ARGV==5;
open IN, "$ARGV[0]" or die $!;
open IN2, "$ARGV[1]" or die $!;
open IN3, "$ARGV[2]" or die $!;
open IN4, "$ARGV[3]" or die $!;
open OUT, ">$ARGV[4]" or die $!;

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

foreach my $k(keys %anno){
	print OUT "$k\t$anno{$k}\n";
}
close OUT;

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


my %bg;
while(<IN3>){
	chomp;
	if(/^#/){
		next;
	}
	next if(/^$/);
	next if(/^---/);
	my @tmp=split /\t/;
	$bg{$tmp[0]}=$tmp[8];
}
close IN3;
	

my $sharp=0;
while(<IN4>){
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
	my @arr=split /\//,$hyperlink;
	shift @arr;
	foreach my $e(@arr){
		my @ko = split /\%09/, $e;
		if(exists $color{$ko[0]}){
			my $str = $ko[0].'%09'.$color{$ko[0]};
			my $term = $ko[0].'%09red';
			$tmp[-1] =~ s/$term/$str/;
		}
	}
	
	# add bg ko
	if(exists $bg{$tmp[0]}){
		my @ko1=split /\|/, $tmp[8];
		my %ko1;
		foreach my $e(@ko1){
			$ko1{$e}='';
		}
		my $kos=$bg{$tmp[0]};
		my @kos=split /\|/, $kos;
		my %kos;
		foreach my $e(@kos){
			unless(exists $ko1{$e}){
				
				$kos{$e}='';
			}
		}
		foreach my $e(%kos){
			my $k='/'.$e.'%09pink';
			$tmp[-1] .= $k;
			$tmp[-1] =~ s/\/%09pink//;
		}

	}
	
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
close IN4;
				
