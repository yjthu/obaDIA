#!/bin/bash -e

SWISSPROT_SQLITE_DB_URL="https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz";

for file in data/*.gz
do
  if [ ! -e ${file%.gz} ]; then
      gunzip -c $file > ${file%.gz}
  fi
done

if [ ! -d edgeR_trans ]; then
    tar xvf data/edgeR_trans.tar
fi

if [ ! -d edgeR_genes ]; then
    tar xvf data/edgeR_genes.tar
fi

BOILERPLATE="Trinotate.boilerplate.sqlite"

if [ -e $BOILERPLATE ]; then
    echo $BOILERPLATE
    rm $BOILERPLATE
fi

if [ ! -s $BOILERPLATE.gz ]; then
    echo pulling swissprot resource db from ftp site
    wget $SWISSPROT_SQLITE_DB_URL -O $BOILERPLATE.gz
fi

gunzip -c $BOILERPLATE.gz > $BOILERPLATE


sqlite_db="myTrinotate.sqlite"

if [ $* ]; then
    sqlite_db="/tmp/myTrinotate.sqlite"
fi

cp  $BOILERPLATE ${sqlite_db}

echo "###############################"
echo Loading protein set
echo "###############################"

../Trinotate ${sqlite_db} init --gene_trans_map data/Trinity.fasta.gene_to_trans_map --transcript_fasta data/Trinity.fasta --transdecoder_pep data/Trinity.fasta.transdecoder.pep



echo "##############################"
echo Loading swissprot blast results
echo "##############################"

../Trinotate ${sqlite_db} LOAD_swissprot_blastp data/swissprot.blastp.outfmt6
../Trinotate ${sqlite_db} LOAD_swissprot_blastx data/swissprot.blastx.outfmt6


echo "####################################"
echo Loading Custom Database Blast Results
echo "####################################"

# blastP
../Trinotate ${sqlite_db} LOAD_custom_blast --outfmt6 data/custom_pombe.blastp.outfmt6 --prog blastp --dbtype custom_pombe_pep
# blastX, same database
../Trinotate ${sqlite_db} LOAD_custom_blast --outfmt6 data/custom_pombe.blastx.outfmt6 --prog blastx --dbtype custom_pombe_pep
 

echo "#############################"
echo Loading PFAM results
echo "#############################"

../Trinotate ${sqlite_db} LOAD_pfam data/TrinotatePFAM.out 


echo "############################"
echo Loading TMHMM results
echo "############################"

../Trinotate ${sqlite_db} LOAD_tmhmm data/tmhmm.out

echo "###########################"
echo Loading SignalP results
echo "###########################"

../Trinotate ${sqlite_db} LOAD_signalp data/signalp.out


echo "###########################"
echo Loading RNAMMER results
echo "###########################"

../Trinotate ${sqlite_db} LOAD_rnammer data/Trinity.fasta.rnammer.gff


#################################################################
## Load Expression info and DE analysis results for Trinotate-web
#################################################################


# import the expression data (counts, fpkms, and samples)

echo "###################################################"
echo Loading Component Expression Matrix and DE results
echo "###################################################"

# expression data load for genes
../util/transcript_expression/import_expression_and_DE_results.pl --sqlite ${sqlite_db} --gene_mode \
        --samples_file data/samples.txt \
        --count_matrix data/Trinity_genes.counts.matrix \
        --fpkm_matrix data/Trinity_genes.TMM.EXPR.matrix 

# DE results load for genes
../util/transcript_expression/import_expression_and_DE_results.pl --sqlite ${sqlite_db} --gene_mode \
        --samples_file data/samples.txt \
        --DE_dir edgeR_genes


echo "##################################################"
echo Loading Transcript Expression Matrix and DE results
echo "##################################################"

# expression data load for transcripts
../util/transcript_expression/import_expression_and_DE_results.pl --sqlite ${sqlite_db} --transcript_mode \
        --samples_file data/samples.txt \
        --count_matrix data/Trinity_trans.counts.matrix \
        --fpkm_matrix data/Trinity_trans.TMM.EXPR.matrix

# DE results load for transcripts
../util/transcript_expression/import_expression_and_DE_results.pl --sqlite ${sqlite_db} --transcript_mode \
        --samples_file data/samples.txt \
        --DE_dir edgeR_trans


echo "######################################################"
echo Loading transcription profile clusters for transcripts
echo "######################################################"


# import the transcription profile cluster stuff
../util/transcript_expression/import_transcript_clusters.pl --group_name DE_all_vs_all --analysis_name diffExpr.P0.1_C1.matrix.RData.clusters_fixed_P_60 --sqlite ${sqlite_db} edgeR_trans/diffExpr.P0.1_C1.matrix.RData.clusters_fixed_P_60/*matrix


echo "###########################"
echo Generating report table
echo "###########################"

../Trinotate ${sqlite_db} report > Trinotate_report.xls


echo "#########################################"
echo Extracting Gene Ontology Mappings Per Gene
echo "#########################################"

../util/extract_GO_assignments_from_Trinotate_xls.pl  --Trinotate_xls Trinotate_report.xls -G -I > Trinotate_report.xls.gene_ontology

# Load annotations
../util/annotation_importer/import_transcript_names.pl ${sqlite_db} Trinotate_report.xls



# Generate trinotate report summary statistics
../util/report_summary/trinotate_report_summary.pl Trinotate_report.xls Trinotate_report_stats


echo "##########################"
echo done.  See annotation summary file:  Trinotate_report.xls
echo "##########################"

