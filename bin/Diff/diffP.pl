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

		-expr 		<s> : protein Abundance table
		-conf 		<s> : mapDIA config file
		-o 			<s> : output dir
		-group 		<s> : group that the samples belong, e.g. sample1/sample2,sample3
		-groupname 	<s> : group names, default group1,group2,...
		-compare 	<s> : how to compare the groups, e.g. 1:2,1:3; Please input in order of trait:control 1vs2
		-fc_cutoff 	<f> : default: 1
		-fdr_cutoff <f> : default: 0.1
		-h|?|help   :  Show this help
=========================================================================
USAGE
}

my ($help, $expr, $conf, $outdir, $group, $groupname, $compare, $fc_cutoff, $fdr_cutoff);
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
print OUT "mkdir $outdir/tableProcess\n";
print OUT "perl $bin/process_expr.pl $expr $outdir/tableProcess/DIA $group\n";
print OUT "perl $bin/nonZero.pl $outdir/tableProcess/DIA.normalized.xls >$outdir/tableProcess/DIA.count.xls\n";
print OUT "perl $bin/calrowmeans_v1.pl -rpkm $outdir/tableProcess/DIA.normalized.xls -group $group -groupname $gn -out-rowmeans $outdir/tableProcess/DIA.means.xls\n";

print OUT "\necho ================== mapDIA ==================\n";
print OUT "mkdir $outdir/mapDIA $outdir/diff_results\n";
print OUT "sed s'#^FILE=.*#FILE=$outdir/tableProcess/DIA.normalized.xls#' $conf >$outdir/mapDIA/mapDIA.conf\n";
print OUT "cd $outdir/mapDIA\n";
print OUT "$mapDIA $outdir/mapDIA/mapDIA.conf\n";
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


