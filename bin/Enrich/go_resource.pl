#!/usr/bin/perl
use strict;
use warnings;
use GO::Parser;

my $db_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db';
my $obo_raw="$db_dir/gene_ontology.1_2.obo";
	
my $parser=new GO::Parser({handler=>'obj',use_cache=>1});
$parser->parse($obo_raw);
my $graph=$parser->handler->graph;
my $it=$graph->create_iterator;
my %allGOs=();
while(my $ni=$it->next_node_instance){
		my $acc=$ni->term->acc;
		my $name=$ni->term->name;
		my $onto=$ni->term->namespace();
		$allGOs{$acc}=$name."\t".$onto;
}

open GOANN,$ARGV[0] or die $!;

        my %go_genes=();
        my %genes_go=();
        my $goInfo="gene2GO<-list(";
        my $geneInfo="names(gene2GO)=c(";
                while(<GOANN>){
                        chomp;
                        my @temp=split /\t/;
                        my $geneID=shift @temp;
                        my %gene_GOs=();
                        foreach(@temp){
                                next unless defined($allGOs{$_});
                                my $ref= $graph->get_reflexive_parent_terms($_);
                                foreach my $term_obj(@$ref){
                                        my $acc=$term_obj->acc();
                                        $gene_GOs{$acc}=1;
                                        $go_genes{$acc}{$geneID}=1;
                                        $genes_go{$geneID}{$acc}=1
                                }
                        }
                        my $goID = join ('","',keys %gene_GOs);
                        $goInfo=$goInfo."\"".$geneID."\" = c(\"".$goID."\"),";
                        $geneInfo .="\"$geneID\",";
                }
        close GOANN;

	$goInfo=~ s/,$//;
        $goInfo .=");\n";

        $geneInfo =~ s/,$//;
        $geneInfo .= ");\n";

open OUT2, ">$ARGV[1].resource";
print OUT2 "$goInfo\n$geneInfo\n";
close OUT2;
