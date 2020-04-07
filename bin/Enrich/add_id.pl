#!/usr/bin/perl
die "perl $0 <ko annotation from annotate.py> <enrichment result from identify.py> > <output>\n" unless @ARGV==2;
my $annot_result=shift;
my $path_result=shift;
my %ko;
my %Entrez;
my %Ensembl;
my %symbol;

open AN, $annot_result or die "cannot open $annot_result !";
#print "open the annot file the first time\n";
while (<AN>) {
	chomp;
	if (/\/\/\/\//) {
#		print "finished the first time reading annot file\n";
		last;
	}
	if (/^\w/) {
		my @tmp=split /\t/;
		my @IDs=split /\|/,$tmp[1];
		if ($IDs[0] eq "None") {
			$ko{$tmp[0]}=" ";
			$Entrez{$tmp[0]}=" ";
			$Ensembl{$tmp[0]}=" ";
			$symbol{$tmp[0]}=" ";
		}
		else
		{
			$ko{$tmp[0]}=$IDs[0];
			$Entrez{$tmp[0]}=" ";
			if ($IDs[0]=~/\:(.+)/) {
				$Ensembl{$tmp[0]}=$1;
			}
			if ($IDs[1] eq "") {
				$symbol{$tmp[0]}=" ";
			}
			else
			{$symbol{$tmp[0]}=$IDs[1];}
		}
	}
}
close AN;

open AN, $annot_result or die "cannot open $annot_result !";
#print "open the annot file the second time\n";
$/="\/\/\/\/\n";
#print "change the line break symbol\n";
<AN>;
#print "read once\n";
while (<AN>) {
	chomp;
	my @line=split /\n/;
	if (@line >=3) {
		my @tmp1=split /\t/,$line[0];
#		my @tmp3=split /\t/,$line[2];
		if($line=~/^Entrez Gene ID:\s+(\d+)/){
			$Entrez{$tmp1[1]}=$1;
		}
	}
	else
	{next;}
}
close AN;

$/="\n";
#print "change the line break symbol back!\n";
open PA, $path_result or die "cannot open $path_result !";
while (<PA>) {
	chomp;
	if (/^\#Term/) {
#		print "Oh, find the title\n";
		my @tmp=split /\t/;
		print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\tKO\t$tmp[8]\n";
	}
	elsif (/^[\w@]/) {
		my @tmp=split /\t/;
		my @id=split/\|/,$tmp[7];
		my $ko_id="";
		my $Entrez_id="";
		my $Ensembl_id="";
		my $symbol_id="";
		foreach my $i (@id) {
			$ko_id .=$ko{$i}."\|";
			$Entrez_id .=$Entrez{$i}."\|";
			$Ensembl_id .=$Ensembl{$i}."\|";
			$symbol_id .=$symbol{$i}."\|";
		}
	$ko_id =~ s/\|$//;
   if ( $tmp[2] !~ /ko01100|ko01110|ko01120/ ) {
       print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$ko_id\t$tmp[8]\n";
     }
	}
	else
	{
		print "$_\n";
	}
}
close PA;
