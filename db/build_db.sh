#!/usr/bin/bash

# download
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz

# unpack
gunzip gene_ontology.1_2.obo.gz
gunzip go.type.gz
gunzip ko00001.keg.gz
gunzip NOG.annotations.tsv.gz

gunzip uniprot_sprot.fasta.gz
gunzip Pfam-A.hmm.gz
gunzip Pfam-A.clans.tsv.gz

# build
../src/diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
../src/hmmpress Pfam-A.hmm
cut -f 1,5 Pfam-A.clans.tsv >Pfam-A.anno

