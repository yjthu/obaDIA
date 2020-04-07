#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
die "Usage: perl $0 <annotation from kobas> <output dir>\n" unless @ARGV==2;

my $db_dir="/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db";
my $data="$db_dir/ko00001.keg";

my $annotation=shift;
my $out=shift;
my $dir=getcwd();

unless(-d $out){
	!system("mkdir $out") or warn "Something goes wrong: $!\n";
}

###	extract ko id, name and definition from ko flat file
my %ko=();
my $annotated_genes=0;
open ANN,$annotation;
<ANN>;
while(<ANN>){
	chomp;
	my @tmp=split;
	my $gid=$tmp[0];
	$ko{$tmp[1]}{$gid}=1;
	$annotated_genes++;
	
}
close ANN;

### 
my ($f,$s,$t,$pid,$k);
my %path2gene1=();
my %path2gene2=();
open DATA,$data;
while(<DATA>){
	if(/^A\<b\>Human Diseases\<\/b\>/){
		last;
	}
	if(/^A\<b\>(.+)\<\/b\>/){
		$f=$1;
		$s="";
		$t="";
	}elsif(/^B\s+\<b\>(.+)\<\/b\>/){
		$s=$1;
		$t="";
	}elsif(/^C\s+\d+\s(.+)/){
		if($1=~/(.*)\s\[PATH:(ko\d+)\]/){
			$t=$1;
			$pid=$2;
		}
	}elsif(/^D\s+(K\d+)\s+(.+)/){
		$k=$1;
#		$an=$2;
		next unless(exists($ko{$k}));
		foreach(keys %{ $ko{$k}}){
			$path2gene1{$f."\t".$s."\t".$t."\t".$pid}{$_}=1;
			$path2gene2{$f."\t".$s}{$_}=1;
		}
	}
}
open STAT,">$out/KEGG_classification.xls";
print STAT "Pathway Hierarchy1\tPathway Hierarchy2\tKEGG Pathway\tPathway ID\tNumber\tIDs\n";
foreach my $key1(sort keys %path2gene1){
	my @tmp=keys %{ $path2gene1{$key1}};
	my $c=@tmp;
	print STAT $key1."\t".$c."\t".join(",",@tmp)."\n";
}
close STAT;

open COUNT,">$out/KEGG_classification_count.txt";
print COUNT "###\tAnnotated Gene Number: $annotated_genes\n";
print COUNT "#Pathway Hierarchy1\tPathway Hierarchy2\tGene Number\n";
foreach my $key2(sort keys %path2gene2){
	my @tmp=keys %{ $path2gene2{$key2}};
	my $c=@tmp;
	print COUNT $key2."\t".$c."\n";
}
close COUNT;

my $R= <<__R__;
##=======================================================================================================================
setwd("$out");
dat<-read.table("KEGG_classification_count.txt",skip=1,sep="\\t")
dat<-dat[order(dat[,1]),]
dat<-subset(dat,dat[,1] !="Human Diseases")
percent=round(dat[,3]/$annotated_genes*100)

labels=as.factor(as.character(dat[,1]))
class<-levels(labels)
n<-floor(max(percent)*1.5/5)*5
a<-sapply(class,function(x) grep(x,labels)[[1]])
b<-c(a[-1],length(percent))
A<-a*1.75-0.5
B<-b*1.75-1.25
pos=(a+b)/2

cols<-rainbow(15)[(1:length(class))*3]
levels(labels)<-cols

png("KEGG_classification.png",width=12,height=9,res=72*4,units="in",type="cairo-png")
par(mar=c(5,20,4,4))
barplot(percent,beside=TRUE,col=as.character(labels),horiz=TRUE,xlim=c(0,max(percent)*1.5),space=0.75,main="KEGG Classification")
text(percent+0.5,(1:length(percent))*1.75-0.5,labels=dat[,3],cex=0.75)
text(-0.1,(1:length(percent))*1.75-0.5,labels=dat[,2],pos=2,xpd=TRUE)
mtext(side=1,line=2,"Percentage")
for(i in 1:length(class)){
	segments(n,A[i],n,B[i],lwd=4,col=cols[i],xpd=TRUE)
}
text(n+0.3,pos*1.75,labels=LETTERS[1:length(class)],xpd=TRUE,col="black",pos=4)
dev.off()

pdf("KEGG_classification.pdf",width=12,height=9)
par(mar=c(5,20,4,4))
barplot(percent,beside=TRUE,col=as.character(labels),horiz=TRUE,xlim=c(0,max(percent)*1.5),space=0.75,main="KEGG Classification")
text(percent+0.5,(1:length(percent))*1.75-0.5,labels=dat[,3],cex=0.75)
text(-0.1,(1:length(percent))*1.75-0.5,labels=dat[,2],pos=2,xpd=TRUE)
mtext(side=1,line=2,"Percentage")
for(i in 1:length(class)){
	segments(n,A[i],n,B[i],lwd=4,col=cols[i],xpd=TRUE)
}
text(n+0.3,pos*1.75,labels=LETTERS[1:length(class)],xpd=TRUE,col="black",pos=4)
dev.off()


##=======================================================================================================================
__R__

open R,"|R --slave --vanilla" or die $!;
print R $R;
close R;

open R,">$out/KEGG_classification.R";
print R $R;
close R;
