#!/usr/bin/env perl
use strict;
use warnings;
use GO::Parser;
die "Usage: Classify genes in the term of gene ontology function.\nExample: perl $0 <gene2go> <output-dir>.\n" unless @ARGV==2;

my $db_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db';
my $obo="$db_dir/gene_ontology.1_2.obo";       # obo file

my $goann=shift;
my $out=shift;
unless(-d $out){
	!system("mkdir $out") or die $!;
}
my $parser=new GO::Parser({handler=>'obj',use_cache=>1});
$parser->parse($obo);
my $graph=$parser->handler->graph;

my %exists_go=();
my $it=$graph->create_iterator;
while(my $ni=$it->next_node()){
        $exists_go{$ni->acc}=1;
}


my @root_go=qw(GO:0008150 GO:0005575 GO:0003674);       # root GO terms
my %go_lev1 = ( 'GO:0008150' => 'Biological Process',
                'GO:0005575' => 'Cellular Component',
                'GO:0003674' => 'Molecular Function'
);
my %Namespaces = ( 'biological_process' => 'Biological Process',
                   'cellular_component' => 'Cellular Component',
                   'molecular_function' => 'Molecular Function'
);

my %go_lev2=();
my %go_lev3=();
my %go_lev4=();

my %path=();

foreach my $termLev1(@root_go){
	foreach my $termLev2(@{$graph->get_child_terms($termLev1)}){
		
		$go_lev2{$termLev2->acc}=$termLev2->name;
		foreach my $termLev3(@{$graph->get_child_terms($termLev2)}){
			$go_lev3{$termLev3->acc}=$termLev3->name;
			foreach my $termLev4(@{$graph->get_child_terms($termLev3)}){
				$go_lev4{$termLev4->acc}=$termLev4->name;
				$path{$termLev1}{$termLev2->acc}{$termLev3->acc}{$termLev4->acc}=1;
			}
		}
	}
}

my %gene2go=();
open GO,$goann;
my %genelist=();
while(<GO>){
	chomp;
	my ($gene,@GOs)=split/\t/;
	foreach my $go(@GOs){
		next unless exists($exists_go{$go});
		foreach my $tmp(@{$graph->get_recursive_parent_terms($go)}){
			$gene2go{$tmp->acc}{$gene}=1;
		}
	}
	$genelist{$gene}=1 if @GOs>0;
}
close GO;
my $total_gene=keys %genelist;

open LV2,">$out/GO_classification_count.txt";
print LV2 "## Total annotated genes: $total_gene\n";
print LV2 "#GO ID (Lev2)\tGO Term (Lev2)\tGO Term (Lev1)\tNumber\n";
open LV4,">$out/GO_classification.xls";
print LV4 "GO ID (Lev1)\tGO Term (Lev1)\tGO ID (Lev2)\tGO Term (Lev2)\tGO ID (Lev3)\tGO Term (Lev3)\tGO ID (Lev4)\tGO Term (Lev4)\tNumber\tList\n";
foreach my $lv1(keys %path){
	foreach my $lv2(keys %{ $path{$lv1}}){
		if(exists $gene2go{$lv2}){
			my $n2=keys %{$gene2go{$lv2}};
			print LV2 $lv2."\t".$go_lev2{$lv2}."\t",$go_lev1{$lv1}."\t".$n2."\n";
		}
		foreach my $lv3(keys %{ $path{$lv1}{$lv2}}){
			foreach my $lv4(keys %{ $path{$lv1}{$lv2}{$lv3}}){
				if(exists $gene2go{$lv4}){
					my $n4=keys %{$gene2go{$lv4}};
					my $l4=join(",",keys %{$gene2go{$lv4}});
					print LV4 $lv1."\t".$go_lev1{$lv1}."\t".$lv2."\t".$go_lev2{$lv2}."\t".$lv3."\t".$go_lev3{$lv3}."\t".$lv4."\t".$go_lev4{$lv4}."\t".$n4."\t".$l4."\n";
				}
			}
		}
	}

}
close LV2;
close LV4;

my $R= <<"__RMAIN__";
#====================================================================================================
setwd("$out")
count<-read.table("GO_classification_count.txt",sep="\\t",skip=1)
colnames(count)<-c("go","name","class","number")
count<-subset(count,number>$total_gene*0.001)
count<-count[order(count\$name),]
count<-count[order(count\$class),]
labels<-as.character(count\$name)
charLn<-max(nchar(labels))

bar_dat<-count\$number*100/$total_gene
ylab<-round(c(0.1,1,10,100)*$total_gene/100)

lns<-nchar(labels)#labels的字节数
lns[lns>45]<-45#超过45的都变成45
x1<-1.75*(1:length(bar_dat))-11
y1<-rep(0.00002,length(bar_dat))
x0<-sapply(1:length(bar_dat),function(x) 1.75*x-lns[x]/5-0.5)
y0<-sapply(1:length(bar_dat),function(x) 10^(log10(0.07)+(log10(0.00002)-log10(0.07))*lns[x]/50))

n1<-grep("Biological",count\$class)[1]
n2<-grep("Cellular",count\$class)[1]
n3<-grep("Molecular",count\$class)[1]
nodes<-c(n1,n2-1,n2,n3-1,n3,length(bar_dat))

pdf("GO_classification_bar.pdf",width=(length(labels)*1.75+14)*0.2,height=3+0.2*(charLn/2+4))
par(mar=c(charLn/2,10,4,5))
barplot(bar_dat,log="y",col="OliveDrab",yaxt="n",xaxt="n",space=0.75,ylim=c(0.1,120),main="Function Classification (GO)")
axis(side=2,at=c(0.1,1,10,100),labels=c(0.1,1,10,100))
axis(side=4,at=c(0.1,1,10,100),labels=ylab)
axis(side=1,at=(0:(length(labels)-1))*1.75+0.5,labels=rep("",length(labels)))
text(1:length(labels)*1.75+0.25,0.07,labels=labels,xpd=TRUE,srt=60,pos=2)

mtext(side=2,line=3,"Percent of proteins")
mtext(side=4,line=3,"Number of proteins")
for(a in nodes){
        segments(x0[a],y0[a],x1[a],y1[a],xpd=TRUE,lwd=1.5)
}
for(b in c(1,3,5)){
        segments(x1[nodes[b]],0.00002,x1[nodes[b+1]],0.00002,xpd=TRUE,lwd=1.5)
        text(x1[floor((nodes[b]+nodes[b+1])/2)],0.000018,count[floor((nodes[b]+nodes[b+1])/2),3],xpd=T,pos=1)
}
dev.off()

png("GO_classification_bar.png",width=(length(labels)*1.75+14)*0.2,height=3+0.2*(charLn/2+4),type="cairo-png",res=72*4,units="in")
par(mar=c(charLn/2,10,4,5))
barplot(bar_dat,log="y",col="OliveDrab",yaxt="n",xaxt="n",space=0.75,ylim=c(0.1,120),main="Function Classification (GO)")
axis(side=2,at=c(0.1,1,10,100),labels=c(0.1,1,10,100))
axis(side=4,at=c(0.1,1,10,100),labels=ylab)
axis(side=1,at=(0:(length(labels)-1))*1.75+0.5,labels=rep("",length(labels)))
text(1:length(labels)*1.75+0.25,0.07,labels=labels,xpd=TRUE,srt=60,pos=2)

mtext(side=2,line=3,"Percent of proteins")
mtext(side=4,line=3,"Number of proteins")
for(a in nodes){
        segments(x0[a],y0[a],x1[a],y1[a],xpd=TRUE,lwd=1.5)
}
for(b in c(1,3,5)){
        segments(x1[nodes[b]],0.00002,x1[nodes[b+1]],0.00002,xpd=TRUE,lwd=1.5)
        text(x1[floor((nodes[b]+nodes[b+1])/2)],0.000018,count[floor((nodes[b]+nodes[b+1])/2),3],xpd=T,pos=1)
}
dev.off()


#====================================================================================================
__RMAIN__

open R,"|R --vanilla --slave" or die $!;
print R $R;
close R;

open R,">$out/GO_classification.R";
print R $R;
close R;

