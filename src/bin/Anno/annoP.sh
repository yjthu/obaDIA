#!/usr/bin/bash

hmmscan='hmmscan'
signalp='signalp'
diamond='diamond'
Rscript='Rscript'
threads=40

base_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
db_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db

TRINOTATE_HOME=$base_dir/src/Trinotate-v3.1.0-pro
uniprot_sprot=$db_dir/uniprot_sprot.fasta
pfam_db=$db_dir/Pfam-A.hmm
Trinotate_sqlite=$db_dir/Trinotate.sqlite.gz
bin=$base_dir/bin/Anno

fa=$1


$diamond blastp -d $uniprot_sprot.dmnd -q $fa -o blastp.outfmt6 --outfmt 6 --tmpdir . -k 1 --quiet --more-sensitive
$hmmscan --cpu $threads --domtblout TrinotatePFAM.out $pfam_db $fa > pfam.log
$signalp -T . -f short -n signalp.out $fa

grep '>' $fa |sed 's/>//' |awk '{print $1"\t"$1}' >${fa##*/}.gene_trans_map

cp $Trinotate_sqlite Trinotate.sqlite.gz
gunzip -f Trinotate.sqlite.gz
$TRINOTATE_HOME/Trinotate Trinotate.sqlite init --gene_trans_map ${fa##*/}.gene_trans_map --transcript_fasta $fa --transdecoder_pep $fa 
$TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
$TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
$TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_signalp signalp.out
$TRINOTATE_HOME/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls

rm -f Trinotate.sqlite
cut -f 5,8,9,11-14 trinotate_annotation_report.xls >Annotation.txt

cut -f 1,4 Annotation.txt|awk '{if(length($2)>1){print $0}}' >eggNOG.anno
cut -f 1,5 Annotation.txt|awk '{if(length($2)>1){print $0}}' >Kegg.anno
perl $bin/matchEggNOG.pl eggNOG.anno >eggNOG.anno.plotdata
$Rscript $bin/simpleBarplot.R eggNOG.anno.plotdata 6 12
perl $bin/trino2go.pl Annotation.txt >GO.anno
perl $bin/trino2pfam.pl Annotation.txt >PFAM.anno
perl $bin/go_classification_v3.pl GO.anno GO.anno.output
perl $bin/matchKegg.pl Kegg.anno >Kegg.anno.plotdata
perl $bin/kegg_classification_v3.pl Kegg.anno.plotdata Kegg.anno.output
perl $bin/formatAnno.pl Annotation.txt| perl $bin/addAnno.pl eggNOG.anno.plotdata Kegg.anno.plotdata - >Annotation.xls
