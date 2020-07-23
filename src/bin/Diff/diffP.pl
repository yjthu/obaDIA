#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description    Proteomics Differential expression analysis

Options

		-expr 		<s> : Abundance matrix, can be protein-level, peptide-level or fragment-level
		-conf 		<s> : mapDIA config file
		-o 		<s> : output dir
		-group 		<s> : group that the samples belong, e.g. sample1/sample2,sample3
		-groupname 	<s> : group names, default group1,group2,...
		-compare 	<f> : how to compare the groups, e.g. 1:2,1:3; Please input in order of trait:control 1vs2
		-fc_cutoff 	<f> : default: 1
		-fdr_cutoff	<f> : default: 0.1
		-level		<s> : level of abundance matrix. Choice are prot/pep/frag [default: prot]
		-h|?|help   	    :  Show this help
=========================================================================
USAGE
}

my ($help, $expr, $conf, $outdir, $group, $groupname, $compare, $fc_cutoff, $fdr_cutoff, $level);
GetOptions(
		"h|?|help"=>\$help,
		"expr=s"=>\$expr,
		"conf=s"=>\$conf,
		"o=s"=>\$outdir,
		"group=s"=>\$group,
		"groupname=s"=>\$groupname,
		"compare=s"=>\$compare,
		"fc_cutoff=f"=>\$fc_cutoff,
		"fdr_cutoff=f"=>\$fdr_cutoff,
		"level=s" =>\$level,
);

if(!defined($expr) || !defined($conf) || !defined($group) || defined($help)){
        &usage;
        exit 0;
}

#=================================
#configure
my $base_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA';
my $db_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db';
my $mapDIA="mapDIA";
my $Rscript="Rscript";
my $bin="$base_dir/bin/Diff";
#=================================

$fc_cutoff ||= 1;
$fdr_cutoff ||= 0.1;
$compare ||= "1:2";
$level ||= 'prot';

my @groupnames;
my @array=split /,/,$group;
my $sn=$group;
$sn=~ s/\//,/g;
my @sname=split /,/,$sn;
my %hash=();
foreach(@sname){
        $hash{$_}=1;
}
my @samplename=sort keys %hash;
$sn=join(",",@samplename);


my @com=split /,/,$compare;
my $groupnumber=@array;

$compare=join(",",@com);
if (!defined($groupname)) {
	foreach (my $i=1;$i<=$groupnumber;$i++){
		push @groupnames,"group$i";
	}
}else{
	@groupnames=split /,/,$groupname;
}

`mkdir -p $outdir`;
open OUT, ">$outdir/group.txt" or die $!;
print OUT "sample\tgroup\n";
for(my $i=0; $i<=$#groupnames; $i++){
	my @arr = split '/', $array[$i];
	foreach my $e(@arr){
		print OUT "$e\t$groupnames[$i]\n";
	}
}
close OUT;

open OUT, ">$outdir/compare.txt" or die $!;
my @comp = split /,/, $compare;
for(my $i=0; $i<=$#comp; $i++){
	my $comp1 = int((split /:/, $comp[$i])[0]) - 1;
	my $comp2 = int((split /:/, $comp[$i])[1]) - 1;
	print OUT "$groupnames[$comp1]vs$groupnames[$comp2]\t$groupnames[$comp1]\t$groupnames[$comp2]\n";
}
close OUT;

my $gn=join(',',@groupnames);

###########################	 
my $dir;
unless (-d $outdir) {
	`mkdir -p $outdir`;
}

open OUT, ">$outdir/.diffP.sh" or die $!;

print OUT "echo Start Time:\ndate\n";
print OUT "echo ================== Process Tables ==================\n";
print OUT "mkdir -p $outdir/tableProcess $outdir/mapDIA $outdir/diff_results\n";
#print OUT "perl $bin/process_expr_impute.pl $expr $outdir/tableProcess/DIA $group\n";
#print OUT "perl $bin/process_expr.pl $outdir/tableProcess/DIA.impute.xls $outdir/tableProcess/DIA $group\n";

if($level eq 'prot'){
	print OUT "perl $bin/nonZero.pl $expr >$outdir/tableProcess/DIA.count.xls\n";
	print OUT "perl $bin/process_expr_reorder.pl $expr $outdir/tableProcess/DIA $group\n";
	print OUT "sed s'#^FILE=.*#FILE=$outdir/tableProcess/DIA.order.xls#' $conf >$outdir/mapDIA/mapDIA.conf\n";
	print OUT "cd $outdir/mapDIA\n";
	print OUT "$mapDIA $outdir/mapDIA/mapDIA.conf\n";
	print OUT "perl $bin/process_expr_log2tr.pl $outdir/mapDIA/log2_data.txt >$outdir/tableProcess/DIA.normalized.xls\n";
}

if($level eq 'pep'){
	print OUT "sed 's/\t/___/' $expr | perl $bin/nonZero.pl - >$outdir/tableProcess/frag_pep.count.xls\n";
	print OUT "sed 's/\t/___/' $expr | perl $bin/process_expr_reorder.pl - $outdir/tableProcess/frag_pep $group\n";
	print OUT "sed 's/___/\t/' $outdir/tableProcess/frag_pep.order.xls | sed 's/^ID\t/ID\tpep\t/' >$outdir/tableProcess/frag_pep.ab.xls \n";
	print OUT "sed s'#^FILE=.*#FILE=$outdir/tableProcess/frag_pep.ab.xls#' $conf | sed s'#LEVEL=1#LEVEL=2#' | cat - $bin/fragPepConf >$outdir/mapDIA/mapDIA.conf\n";
	print OUT "cd $outdir/mapDIA\n";
	print OUT "$mapDIA $outdir/mapDIA/mapDIA.conf\n";
	print OUT "perl -lane 'print join(\"\\t\",\@F[0..\$#F-1])' $outdir/mapDIA/protein_level.txt >$outdir/tableProcess/DIA.normalized.xls\n";
}
if($level eq 'frag'){
	print OUT "sed 's/\t/___/' $expr | sed 's/\t/___/' | perl $bin/nonZero.pl - >$outdir/tableProcess/frag_pep.count.xls\n";
	print OUT "sed 's/\t/___/' $expr | sed 's/\t/___/' | perl $bin/process_expr_reorder.pl - $outdir/tableProcess/frag_pep $group\n";
	print OUT "sed 's/___/\t/' $outdir/tableProcess/frag_pep.order.xls | sed 's/___/\t/' | sed 's/^ID\t/ID\tpep\tfrag\t/' >$outdir/tableProcess/frag_pep.ab.xls \n";
	print OUT "sed s'#^FILE=.*#FILE=$outdir/tableProcess/frag_pep.ab.xls#' $conf | sed s'#LEVEL=1#LEVEL=3#' | cat - $bin/fragPepConf >$outdir/mapDIA/mapDIA.conf\n";
	print OUT "cd $outdir/mapDIA\n";
	print OUT "$mapDIA $outdir/mapDIA/mapDIA.conf\n";
	print OUT "perl -lane 'print join(\"\\t\",\@F[0..\$#F-2])' $outdir/mapDIA/protein_level.txt >$outdir/tableProcess/DIA.normalized.xls\n";
}


print OUT "\necho ================== diff ==================\n";
print OUT "perl $bin/calrowmeans_v1.pl -rpkm $outdir/tableProcess/DIA.normalized.xls -group $group -groupname $gn -out-rowmeans $outdir/tableProcess/DIA.means.xls\n";
print OUT "perl $bin/gDIAdiffTablesPlus.pl $outdir/mapDIA/analysis_output_wide_format.txt $outdir/tableProcess/DIA.means.xls $groupname $compare $outdir/diff_results $fc_cutoff $fdr_cutoff\n";


foreach my $gp (@com){
	my @tmp = split /:/,$gp;
	my $vs = $groupnames[$tmp[0]-1]."vs".$groupnames[$tmp[1]-1];
	$dir=$outdir."/diff_results/".$vs;
	print OUT "awk '{if(\$7==\"TRUE\"){print \$1}}' $dir/$vs.diff.xls  >$dir/$vs.diff.id\n";
	print OUT "awk '{if(\$7==\"TRUE\"){if(\$4>0){print \$1\"\\tup\"}else{print \$1\"\\tdown\"}}}' $dir/$vs.diff.xls >$dir/$vs.diff.updown\n";
	print OUT "$Rscript $bin/plotVolcano.R $dir/$vs.diff.xls $fdr_cutoff $fc_cutoff\n";
}

print OUT "\necho ================== Plots ==================\n";
print OUT "$Rscript $bin/den_boxplot.R $outdir/tableProcess/DIA.normalized.xls\n";
print OUT "$Rscript $bin/den_boxplot.R $outdir/tableProcess/DIA.means.xls\n";
print OUT "$Rscript $bin/plotCorScatmat.R $outdir/tableProcess/DIA.means.xls\n";
print OUT "$Rscript $bin/plotCorScatmat.R $outdir/tableProcess/DIA.normalized.xls\n";
close OUT;
`chmod +x $outdir/.diffP.sh`;


