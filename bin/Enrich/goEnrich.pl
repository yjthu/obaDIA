#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
sub usage
{
        print STDERR <<USAGE;
==================================================================================
Description run GO enrichment test
Version:2020-03-01

Options
<Required>
	-pfam		<s> :    gene2pfamfile GO.anno
	-info		<s> :	 gene length file
	-s			<s> :    names of the gene sets,sep=","
	-deg		<s>	:    file containing a list of genes of interest, i.e. DEGs,
				corresponding to the names of the gene sets
<Optional>
    -o 			<s> :    the output directory, default pwd
	-des 		<s> :    pfam description file
===============================================================================
USAGE
}

my($pfam, $info, $des, $s, $deg, $dir, $help);
GetOptions(
	"h|?|help"=>\$help,
	"pfam=s" =>\$pfam,
	"info=s" =>\$info,
	"des=s" =>\$des,
	"s=s"=>\$s,
	"deg=s"=>\$deg,
	"o=s"=>\$dir,
);
my $db_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db';
$des ||="$db_dir/go.type";

my %hash=('BP'=>'biological_process',
'MF'=>'molecular_function',
'CC'=>'cellular_component');


if (!defined($pfam)||!defined($info)||!defined($s)||!defined($deg)||defined($help)) {
	&usage;
	exit 0;
}
my $pwd = `pwd`;
$pwd =~ s/\n//;
if (defined($dir)){
	unless (-e $dir){
		`mkdir $dir`;
	}
}else{
	$dir =$pwd;
}


open IN, $pfam;
open OUT, ">$dir/gene2go_formated.txt";
while(my $line=<IN>){
        chomp $line;
        my @array=split /\t/, $line;
        my $id=shift(@array);
        foreach my $term (@array){
                print OUT "$id\t$term\n";
        }
}
close IN;
close OUT;


my $R =<< "END";
#------------------------------------------
library("goseq")

#get the gene2cateroty file
pfam<-read.table(paste("$dir", "gene2go_formated.txt", sep="\/"))
pfam<-data.frame(pfam)
all<-levels(pfam\$V1)

#get the geneLength of annotated genes
length<- read.table("$info",sep="\\t")
gene.len<-data.frame(len=length[,2])
rownames(gene.len)<-length[,1];
all.len<-gene.len[all,1]
names(all.len)<-all

groups<-unlist(strsplit("$s", ","))
degs<-unlist(strsplit("$deg", ","))
for (i in 1:length(groups)){
#get the de/non-de of the annotated genes
genelist <- scan(degs[i], what="character",quiet= TRUE)
genes<-as.integer(all%in%genelist)

#run the enrichment test
pwf<- nullp(genes, bias.data=all.len, plot.fit=TRUE)
rownames(pwf)<-names(all.len)
res<-goseq(pwf, gene2cat=pfam, method = "Wallenius")

#adjust for multiple test 
res\$over_represented_pvalue[res\$over_represented_pvalue==0]=min(res\$over_represented_pvalue[res\$over_represented_pvalue>0])/10000
res\$padj<-p.adjust(res\$over_represented_pvalue,method="BH");
res\$over_represented_pvalue<- signif(res\$over_represented_pvalue,5)
colnames(res)[1]<-"Pfam_id"
res<-res[,-3]

write.table(res,file=paste("$dir",paste(groups[i],"identify.GO", sep="."),sep="\/"),row.names=FALSE,quote=F, sep="\t")
}
END
open R,"|R --vanilla --slave" or die $!;
print R $R;
close R;

my @groups=split /,/,$s;
my @degs=split /,/, $deg;
foreach my $i (0 .. (@groups)-1){
	open IN, $degs[$i];
	my %hash1;
	while(<IN>){
		chomp;
		$hash1{$_}=1;
	}
	close IN;

	my $path = "$dir/gene2go_formated.txt";
	open IN, $path;

	my %hash2=();
	my %hash3=();
	my %hashtemp=();
	my $n;
	my $m;
	my %genename=();
	while(<IN>){
		chomp;
		my @tmp=split /\t/, $_;
		$hash2{$tmp[1]}++;
		$n++ unless exists $hashtemp{$tmp[0]};
		if (exists $hash1{$tmp[0]}){
			$hash3{$tmp[1]}++; 
			if (!exists $genename{$tmp[1]}){
				$genename{$tmp[1]}=$tmp[0];
			}else{
				$genename{$tmp[1]}=$genename{$tmp[1]}.",".$tmp[0];
			}
			$m++ unless exists $hashtemp{$tmp[0]};
		}
		$hashtemp{$tmp[0]}++;
	}
	close IN;

	open IN, $des;
	my %hash4=();
	while(<IN>){
		chomp;
		my @tmp1=split /\t/, $_;
		$hash4{$tmp1[0]}=$tmp1[1];
	}

	open IN, "$dir/$groups[$i].identify.GO";
	<IN>;
	open OUT, ">$dir/$groups[$i].identify.GO.xls";
	print OUT "#Term\tDatabase\tID\tTerm_type\tInput_number\tBackground_number\tP_Value\tCorrected_P_Value\tInput\tHyperlink\n";

	
	while(my $q=<IN>){
		chomp $q;
		my @temp=split /\t/, $q;
 		my $first = $hash3{$temp[0]};
		unless (exists $hash3{$temp[0]}){
			$first = 0;
		}
		unless (exists $hash4{$temp[0]}){
			$hash4{$temp[0]}="NA";
		}
		unless (exists $genename{$temp[0]}){
			$genename{$temp[0]}="NA";
		}
		$genename{$temp[0]}=~s/,/\|/g;
		my $link = 'http://amigo.geneontology.org/amigo/search/ontology?q='.$temp[0];
		my $type=$temp[5];
		if(exists $hash{$temp[5]}){
			$type=$hash{$temp[5]};
		}
		if($temp[4] ne 'NA'){
			print OUT "$temp[4]\tGene Ontology\t$temp[0]\t$type\t$temp[2]\t$temp[3]\t$temp[1]\t$temp[6]\t$genename{$temp[0]}\t$link\n";
		}
	}
	close IN;
	close OUT;
}
