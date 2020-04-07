#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description     Generate a shell script to run Proteomics Diff analysis

Options

                -group <s> : group that the samples belong, e.g. sample1/sample2,sample3/sample4
                -groupname <s> : group names, default group1,group2,...
                -compare <s> : how to compare the groups, e.g. 1:2,1:3; Please input in order of trait:control 1vs2
                -h|?|help   :  Show this help
=========================================================================
USAGE
}

my ($help, $group, $groupname, $compare);
GetOptions(
                "h|?|help"=>\$help,
                "group=s"=>\$group,
                "groupname=s"=>\$groupname,
                "compare=s"=>\$compare,
);

if(!defined($group) || defined($help)){
        &usage;
        exit 0;
}

$compare ||= "1:2";

my @group=split /,/, $group;
my @groupname=split /,/, $groupname;
my @compare=split /,/, $compare;

my $group_num=@groupname;
my $rep_num=(split /\//, $group[0]);

my @array;
for(my $i=0; $i<$group_num; $i++){
	for(my $j=0; $j<$group_num; $j++){
		$array[$i][$j]=0;
		if($i==$j){
			$array[$i][$j]='-';
		}
	}
}	
foreach my $com(@compare){
	my $com_row=(split /:/, $com)[0]-1;
	my $com_col=(split /:/, $com)[1]-1;
	$array[$com_row][$com_col]=1;
}

my $com_matrix='';
for(my $i=0; $i<$group_num; $i++){
	my $line='';
        for(my $j=0; $j<$group_num; $j++){
		$line .= "$array[$i][$j] ";
	}
	$line=~s/ $/\n/g;
	$com_matrix .= $line;
}

my $level='1';
my $design='IndependentDesign';
if($rep_num>1){
	$design='replicatedesign';
}	

print <<TXT;

###input file
FILE=DIA.normalized.xls
LEVEL=$level


### Experimental design
EXPERIMENTAL_DESIGN=$design


### Filter
MIN_OBS = 1


### Sample information
LABELS=@groupname
SIZE=$rep_num


### min. max. DE
MIN_DE = .5
MAX_DE =.99


### Contrast matrix for group comparison
CONTRAST=
$com_matrix
TXT

