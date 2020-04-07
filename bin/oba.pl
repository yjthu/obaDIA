#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description:	one-step biological analysis pipeline for quantitative proteomics data 
Author:		yjthu
Copyright:	www.e-omics.com, MDZK, Inc.
Version:	v1.0

Options

	Required:
		-fas		<s> : protein squence file 
		-exp		<s> : protein abundance file
		-out		<s> : output directory. [An absolute path is required!]
		-group 		<s> : group that the samples belong to: s1/s2/s3,s4/s5/s6
		-name	 	<s> : group names: g1,g2	
		-spe 		<s> : kobas3.0 KEGG enrichment species abbr. [refer to "species_abbr.txt" file]
	Optional:	
		-fc 		<f> : default: 1
		-fdr 		<f> : default: 0.1
		-comp 		<s> : how to compare the groups, just use numbers, 1:2		
		-alt	 	<s> : kobas3.0 GO/Rectome enrichment species abbr. [default: the same as "-spe"]
		
		-h|?|help   :  Show this help
=========================================================================
USAGE
}

my ($help, $fa, $expr, $outdir, $group, $groupname, $compare, $fc_cutoff, $fdr_cutoff, $species, $spe_go);
GetOptions(
		"h|?|help"=>\$help,
		"fas=s"=>\$fa,
		"exp=s"=>\$expr,
		"out=s"=>\$outdir,
		"group=s"=>\$group,
		"name=s"=>\$groupname,
		"comp=s"=>\$compare,
		"fc=f"=>\$fc_cutoff,
		"fdr=f"=>\$fdr_cutoff,
		"spe=s"=>\$species,
		"alt=s" =>\$spe_go,

);

if(!defined($fa) || !defined($expr) || !defined($outdir) || !defined($group) || !defined($groupname) || !defined($species) || defined($help)){
        &usage;
        exit 0;
}

#=================================
#configure
my $base_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA';
my $bin="$base_dir/bin";

#=================================
$compare ||= '1:2';
$fc_cutoff ||= 1;
$fdr_cutoff ||= 0.1;
$spe_go ||= $species;
my $ko = 'no';
if($species eq 'ko'){
	$ko = 'yes';
}

my @g=split /,/, $groupname;
my @c=split /,/, $compare;
my @compares;
foreach my $each(@c){
	my $g1=(split /:/, $each)[0]-1;
	my $g2=(split /:/, $each)[1]-1;
	my $g1vsg2=$g[$g1].'vs'.$g[$g2];
	push (@compares,$g1vsg2) ;
}
#check input files
my %fa;
open IN, "$fa" or die $!;
while(<IN>){
	chomp;
	if(/^>/){
		my $id=$_;
		$id =~ s/>//;
		$fa{$id}='';
	}
}
close IN;
my $fa_id='';
foreach my $k(sort keys %fa){
	$fa_id .= ','.$k;
}

my %ab;
open IN, "$expr" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $s=0;
	for(my $i=1; $i<=$#tmp; $i++){
		$s += $tmp[$i];
	}
	if($s == 0){
		die "abundance of $tmp[0] is zero!";
	}
	$ab{$tmp[0]}='';
}
close IN;
my $ab_id='';
foreach my $k(sort keys %ab){
	$ab_id .= ','.$k;
}
if($fa_id ne $ab_id){
	die "input files are not consistent!";
}


`mkdir -p $outdir`;
open OUT, ">$outdir/.step0.sh" or die $!;
print OUT "mkdir -p $outdir/Diff $outdir/Anno $outdir/Enrich\n";

#Diff
print OUT "cd $outdir/Diff\n";
print OUT "perl $bin/Diff/g_mapDIAconf.pl -group $group -groupname $groupname -compare $compare >$outdir/Diff/mapDIA.input\n";
print OUT "perl $bin/Diff/diffP.pl -expr $expr -conf $outdir/Diff/mapDIA.input -o $outdir/Diff -group $group -groupname $groupname -compare $compare -fc_cutoff $fc_cutoff -fdr_cutoff $fdr_cutoff\n";

#Anno
print OUT "cd $outdir/Anno\n";
print OUT "cp $bin/Anno/annoP.sh $outdir/Anno/.annoP.sh\n";

#Enrich
print OUT "cd $outdir/Enrich\n";
print OUT "mkdir -p $outdir/Enrich/All\n";
print OUT "sh $bin/Enrich/fa_len.sh $fa $outdir/Enrich/All/seq.len\n";
print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich/All -samplename All -species $species -diffGeneSeq $fa -goanno $outdir/Anno/GO.anno -pfanno $outdir/Anno/PFAM.anno -summary $outdir/Enrich/All/seq.len -ko no -spe_go $spe_go\n";

foreach my $each(@compares){
	print OUT "mkdir -p $outdir/Enrich/$each\n";
		print OUT "echo \"perl $bin/Enrich/extractMyChr.v2.pl $fa fa $outdir/Diff/diff_results/$each/$each.diff.id >$outdir/Diff/diff_results/$each/$each.diff.seq\" >>$outdir/Diff/.diffP.sh\n";
		print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich/$each -samplename $each -species $species -bg $outdir/Enrich/All/All -diffGeneSeq $outdir/Diff/diff_results/$each/$each.diff.seq -goanno $outdir/Anno/GO.anno -pfanno $outdir/Anno/PFAM.anno -summary $outdir/Enrich/All/seq.len -diffid $outdir/Diff/diff_results/$each/$each.diff.id -updown $outdir/Diff/diff_results/$each/$each.diff.updown -difftb $outdir/Diff/diff_results/$each/$each.diff.xls -ko $ko -spe_go $spe_go\n";
}
close OUT;
`sh $outdir/.step0.sh`;


# all in one step pipeline
open OUT, ">$outdir/OneStep.sh" or die $!;
print OUT "\necho Begin Diff Analysis :\ndate\n";
print OUT "cd $outdir/Diff\n";
print OUT "sh .diffP.sh >.diffP.sh.o 2>.diffP.sh.e\n";
print OUT "\necho Begin Annotaion:\ndate\n";
print OUT "cd $outdir/Anno\n";
print OUT "sh .annoP.sh $fa >.annoP.sh.o 2>.annoP.sh.e\n";
print OUT "\necho Begin Enrichment:\ndate\n";
print OUT "cd $outdir/Enrich/All\n";
print OUT "sh .All.enrichP.sh >.All.enrichP.sh.o 2>.All.enrichP.sh.e\n";
foreach my $each(@compares){
	print OUT "cd $outdir/Enrich/$each\n";
	print OUT "sh .$each.enrichP.sh >.$each.enrichP.sh.o 2>.$each.enrichP.sh.e\n";
}
print OUT "\necho All done!\ndate\n";
close OUT;

`chmod +x $outdir/OneStep.sh`;

