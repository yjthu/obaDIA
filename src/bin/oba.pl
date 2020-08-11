#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description:	obaDIA: one-step biological analysis pipeline for data-independent acquisition 
		and other quantitative proteomics data
Author:		github.com/yjthu
Copyright:	www.e-omics.com, MDZK, Inc.
Version:	v1.1

Options

	Required:
		-fas		<s> : protein squence file, fasta format.
		-exp		<s> : protein abundance file, tsv format.
					protein-level, peptide-level or fragment-level file is acceptable
		-out		<s> : output directory. [An absolute path is required!]
		-group 		<s> : group that the samples belong to. [eg: s1/s2/s3,s4/s5/s6] 
					groups should be seperated by ',', 
					samples within a group should be seperated by '/'.
		-name	 	<s> : the name of each group. [eg: group1,group2]	
					should be seperated by ',' and ordered the same as '-group'.
		-spe 		<s> : KOBAS3.0 KEGG enrichment species abbr. [refer to "db/species_abbr.txt" file]
					dowload KOBAS database from 'ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD'
	Optional:	
		-fc 		<f> : log2(foldChange) cutoff for DE proteins. [default: 1]
		-fdr 		<f> : FDR cutoff for DE proteins by mapDIA. [default: 0.1]
		-comp 		<s> : how to compare the groups, use indexs of the group, [default: '1:2']
					required for multiple comparisons
		-alt	 	<s> : KOBAS3.0 GO/Rectome enrichment species abbr. [default: the same as '-spe']
		-mod		<s> : background choose for enrichment, total(T) or expressed(E). [default: E]
		-thread		<s> : thread number for hmmscan to perform Pfam annotation. [default: 1]
		-level		<s> : level of abundance matrix. Choice are prot/pep/frag [default: prot]
					required for peptide-level or fragment-level abundance
		-fast		    :  fast mode to run the pipeline. a part of time-consuming analysis will be skipped.
		-h|?|help	    :  Show this help
=========================================================================
USAGE
}

my ($help, $fa, $expr, $outdir, $group, $groupname, $compare, $fc_cutoff, $fdr_cutoff, $species, $spe_go, $mod, $thread, $level, $fast);
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
		"mod=s" =>\$mod,
		"thread=s" =>\$thread,
		"level=s" =>\$level,
		"fast" =>\$fast,
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
$mod ||= 'E';
$thread ||= 1;
$level ||= 'prot';


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


`mkdir -p $outdir`;

# match input abundance files to fasta sequences
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


my %ab;
open IN, "$expr" or die $!;
open OUT, ">$outdir/.ab.input" or die $!;
my $head=<IN>;
print OUT "$head";
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	if(exists $fa{$tmp[0]}){
		print OUT "$_\n";
		$ab{$tmp[0]}='';
	}
}
close IN;
close OUT;


open IN, "$fa" or die $!;
open OUT, ">$outdir/.seq.input" or die $!;
$/ = '>';
<IN>;
while(<IN>){
        chomp;
		my @tmp = split /\n/;
        my $id = $tmp[0];
        $id =~ s/\>//;
        if(exists $ab{$id}){
                shift(@tmp);
                my $seq = join('',@tmp);
                print OUT ">$id\n$seq\n";
        }
}
$/ = "\n";
$expr = "$outdir/.ab.input";
$fa = "$outdir/.seq.input";

open OUT, ">$outdir/.step0.sh" or die $!;
print OUT "mkdir -p $outdir/Diff $outdir/Anno $outdir/Enrich\n";

#Diff
print OUT "cd $outdir/Diff\n";
print OUT "perl $bin/Diff/g_mapDIAconf.pl -group $group -groupname $groupname -compare $compare >$outdir/Diff/mapDIA.input\n";

print OUT "perl $bin/Diff/diffP.pl -expr $expr -conf $outdir/Diff/mapDIA.input -o $outdir/Diff -group $group -groupname $groupname -compare $compare -fc_cutoff $fc_cutoff -fdr_cutoff $fdr_cutoff -level $level\n";


#Anno
print OUT "cd $outdir/Anno\n";
if(!defined $fast){
	print OUT "sed 's/threads=40/threads=$thread/'  $bin/Anno/annoP.sh >$outdir/Anno/.annoP.sh\n";
}
else{
	print OUT "sed 's/threads=40/threads=$thread/'  $bin/Anno/annoP.sh |grep -v 'pfam'|grep -v '_v3' >$outdir/Anno/.annoP.sh\n";
}

#Enrich
print OUT "cd $outdir/Enrich\n";
print OUT "mkdir -p $outdir/Enrich/All\n";
print OUT "sh $bin/Enrich/fa_len.sh $fa $outdir/Enrich/All/seq.len\n";
if(!defined $fast){
	print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich/All -samplename All -species $species -diffGeneSeq $fa -goanno $outdir/Anno/GO.anno -pfanno $outdir/Anno/PFAM.anno -summary $outdir/Enrich/All/seq.len -ko no -spe_go $spe_go\n";
}

foreach my $each(@compares){
	print OUT "mkdir -p $outdir/Enrich/$each\n";
		print OUT "echo \"perl $bin/Enrich/extractMyChr.v2.pl $fa fa $outdir/Diff/diff_results/$each/$each.diff.id >$outdir/Diff/diff_results/$each/$each.diff.seq\" >>$outdir/Diff/.diffP.sh\n";
		if(!defined $fast){
			print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich/$each -samplename $each -species $species -bg $outdir/Enrich/All/All -mod $mod -diffGeneSeq $outdir/Diff/diff_results/$each/$each.diff.seq -goanno $outdir/Anno/GO.anno -pfanno $outdir/Anno/PFAM.anno -summary $outdir/Enrich/All/seq.len -diffid $outdir/Diff/diff_results/$each/$each.diff.id -updown $outdir/Diff/diff_results/$each/$each.diff.updown -difftb $outdir/Diff/diff_results/$each/$each.diff.xls -ko $ko -spe_go $spe_go\n";
		}else{
			print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich/$each -samplename $each -species $species -mod $mod -diffGeneSeq $outdir/Diff/diff_results/$each/$each.diff.seq -summary $outdir/Enrich/All/seq.len -diffid $outdir/Diff/diff_results/$each/$each.diff.id -difftb $outdir/Diff/diff_results/$each/$each.diff.xls -ko $ko -spe_go $spe_go\n";
		}
}
close OUT;
`sh $outdir/.step0.sh`;


# all in one step pipeline
open OUT, ">$outdir/OneStep.sh" or die $!;
print OUT "\necho Begin Diff Analysis :\ndate\n";
print OUT "cd $outdir/Diff\n";
print OUT "sh .diffP.sh >.diffP.sh.o 2>.diffP.sh.e\n";
print OUT "\necho Begin Annotation:\ndate\n";
print OUT "cd $outdir/Anno\n";
print OUT "sh .annoP.sh $fa >.annoP.sh.o 2>.annoP.sh.e\n";
print OUT "\necho Begin Enrichment:\ndate\n";
if(!defined $fast){
	print OUT "cd $outdir/Enrich/All\n";
	print OUT "sh .All.enrichP.sh >.All.enrichP.sh.o 2>.All.enrichP.sh.e\n";
}
foreach my $each(@compares){
	print OUT "cd $outdir/Enrich/$each\n";
	print OUT "sh .$each.enrichP.sh >.$each.enrichP.sh.o 2>.$each.enrichP.sh.e\n";
}
print OUT "\necho All done!\ndate\n";
close OUT;

`chmod +x $outdir/OneStep.sh`;

