#!/usr/bin/bash

base_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
db_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db

bin=$base_dir/bin/Enrich


geneSum=$1
pfamAnno=$2
degId=$3
sampleName=$4
outdir=$5

mkdir -p $outdir
perl $bin/pfamEnrich.pl -pfam $pfamAnno -info $geneSum -deg $degId -s $sampleName -o $outdir
perl $bin/extract_top20.pl $outdir/$sampleName.identify.pfam.xls >$outdir/$sampleName.identify.pfam.xls.top20
Rscript $bin/Pathwayscatter2.R $sampleName.PFAM $outdir $outdir/$sampleName.identify.pfam.xls.top20
