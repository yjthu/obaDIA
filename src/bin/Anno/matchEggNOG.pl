#!/usr/bin/perl
use strict;
use warnings;

my $db_dir="/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db";
my $funcfile="$db_dir/NOG.annotations.tsv"; 
my %code=("J"=>"J:Translation, ribosomal structure and biogenesis",
"A"=>"A:RNA processing and modification",
"K"=>"K:Transcription",
"L"=>"L:Replication, recombination and repair",
"B"=>"B:Chromatin structure and dynamics",
"D"=>"D:Cell cycle control, cell division, chromosome partitioning",
"Y"=>"Y:Nuclear structure",
"V"=>"V:Defense mechanisms",
"T"=>"T:Signal transduction mechanisms",
"M"=>"M:Cell wall/membrane/envelope biogenesis",
"N"=>"N:Cell motility",
"Z"=>"Z:Cytoskeleton",
"W"=>"W:Extracellular structures",
"U"=>"U:Intracellular trafficking, secretion, and vesicular transport",
"O"=>"O:Posttranslational modification, protein turnover, chaperones",
"C"=>"C:Energy production and conversion",
"G"=>"G:Carbohydrate transport and metabolism",
"E"=>"E:Amino acid transport and metabolism",
"F"=>"F:Nucleotide transport and metabolism",
"H"=>"H:Coenzyme transport and metabolism",
"I"=>"I:Lipid transport and metabolism",
"P"=>"P:Inorganic ion transport and metabolism",
"Q"=>"Q:Secondary metabolites biosynthesis, transport and catabolism",
"R"=>"R:General function prediction only",
"S"=>"S:Function unknown",
);

my %hash;
open IN, "$funcfile" or die $!;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	if(exists $code{$tmp[4]}){
		$hash{$tmp[1]}=$code{$tmp[4]};
	}
}
close IN;


open IN, "$ARGV[0]" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	if($tmp[1] =~ /^(.*?)\^/){
		my $id = $1;
		if(exists $hash{$id}){
			my $nog = $hash{$id};
			$nog =~ s/\s+/_/g;
			my $nog_s = (split /:/, $nog)[0];
			print "$tmp[0]\t$nog_s\t$nog\n";
		}
	}
}
close IN;

		
	
