#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description:	oba_module: a modular version of obaDIA
Author:		github.com/yjthu
Copyright:	www.e-omics.com, MDZK, Inc.
Version:	v1.1

Options

	Common:
		-m		<s> : choose the module you want to run [diff/anno/enrich]
		-o		<s> : output directory
		-h|?|help	    :  Show this help

	Diff:	
		-expr 		<s> : Abundance matrix, can be protein-level, peptide-level or fragment-level
		-group 		<s> : group that the samples belong, e.g. sample1/sample2,sample3
		-groupname 	<s> : group names, default group1,group2,...
		-compare 	<s> : how to compare the groups, e.g. 1:2,1:3; Please input in order of trait:control 1vs2
		-fc_cutoff 	<f> : default: 1
		-fdr_cutoff 	<f> : default: 0.1
		-level		<s> : level of abundance matrix. Choice are prot/pep/frag [default: prot]
		
	Anno:
		-fa		<s> : protein squence file, fasta format.
		
	Enrich:
		-name 		<s> : sample name for output files
		-species 	<s> : species for kegg or background annotate file
		-fa		<s> : DE protein sequence fasta

=========================================================================
USAGE
}

my ($help, $m, $outdir, $expr, $group, $groupname, $compare, $fc_cutoff, $fdr_cutoff, $level, $name, $species, $fa);
GetOptions(
		"h|?|help"=>\$help,
		"m=s"=>\$m,
		"o=s"=>\$outdir,
		"expr=s"=>\$expr,
		"group=s"=>\$group,
		"groupname=s"=>\$groupname,
		"compare=s"=>\$compare,
		"fc_cutoff=f"=>\$fc_cutoff,
		"fdr_cutoff=f"=>\$fdr_cutoff,
		"species=s"=>\$species,
		"level=s" =>\$level,
		"name=s" =>\$name,
		"fa=s" =>\$fa,
);

if(!defined($m) || !defined($outdir) || defined($help)){
        &usage;
        exit 0;
}

#=================================
#configure
my $base_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA';
my $bin="$base_dir/bin";
unless (-d $outdir) {
	`mkdir -p $outdir`;
}

open OUT, ">$outdir/.step0.sh" or die $!;
if($m eq 'diff'){
	if(!defined($expr) || !defined($group) || !defined($groupname)){
			&usage;
			exit 0;
	}
	$fc_cutoff ||= 1;
	$fdr_cutoff ||= 0.1;
	$compare ||= "1:2";
	$level ||= 'prot';

	print OUT "mkdir -p $outdir/Diff\n";

	print OUT "cd $outdir/Diff\n";
	print OUT "perl $bin/Diff/g_mapDIAconf.pl -group $group -groupname $groupname -compare $compare >$outdir/Diff/mapDIA.input\n";
	print OUT "perl $bin/Diff/diffP.pl -expr $expr -conf $outdir/Diff/mapDIA.input -o $outdir/Diff -group $group -groupname $groupname -compare $compare -fc_cutoff $fc_cutoff -fdr_cutoff $fdr_cutoff -level $level\n";
}

if($m eq 'anno'){
	if(!defined($fa)){
			&usage;
			exit 0;
	}
	
	print OUT "mkdir -p $outdir/Anno\n";

	print OUT "cd $outdir/Anno\n";
	print OUT "cp $bin/Anno/annoP.sh $outdir/Anno/.annoP.sh\n";
	
}

if($m eq 'enrich'){
		if(!defined($fa) || !defined($name) || !defined($species) ){
			&usage;
			exit 0;
	}
	print OUT "mkdir -p $outdir/Enrich\n";
	print OUT "perl $bin/Enrich/enrichP.pl -directory-kobas-output $outdir/Enrich -samplename $name -species $species -spe_go $species -mod T -diffGeneSeq $fa\n";
}

close OUT;
`sh $outdir/.step0.sh`;


# all in one step pipeline
open OUT, ">$outdir/OneStep.sh" or die $!;
if($m eq 'diff'){
	print OUT "\necho Begin Diff Analysis :\ndate\n";
	print OUT "cd $outdir/Diff\n";
	print OUT "sh .diffP.sh >.diffP.sh.o 2>.diffP.sh.e\n";
}
if($m eq 'anno'){
	print OUT "\necho Begin Annotation:\ndate\n";
	print OUT "cd $outdir/Anno\n";
	print OUT "sh .annoP.sh $fa >.annoP.sh.o 2>.annoP.sh.e\n";
}
if($m eq 'enrich'){
	print OUT "\necho Begin Enrichment:\ndate\n";
	print OUT "cd $outdir/Enrich\n";
	print OUT "sh .$name.enrichP.sh >.$name.enrichP.sh.o 2>.$name.enrichP.sh.e\n";
}
print OUT "\necho All done!\ndate\n";
close OUT;

`chmod +x $outdir/OneStep.sh`;

