#!/usr/bin/bash

base_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
db_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db

bin=$base_dir/bin/Enrich
gotype=$db_dir/go.type

geneSum=$1
goAnno=$2
degId=$3
sampleName=$4
outdir=$5

mkdir -p $outdir
perl $bin/goEnrich.pl -pfam $goAnno -info $geneSum -deg $degId -s $sampleName -o $outdir
perl $bin/extract_top20.pl $outdir/$sampleName.identify.GO.xls >$outdir/$sampleName.identify.GO.xls.top20
Rscript $bin/Pathwayscatter2.R $sampleName.GO $outdir $outdir/$sampleName.identify.GO.xls.top20

perl $bin/go_resource.pl $goAnno $outdir/go
Rscript $bin/goBar.R $outdir/$sampleName.identify.GO.xls $outdir $sampleName
Rscript $bin/topGO.R $outdir $outdir/go.resource $outdir/$sampleName.identify.GO.xls $degId $outdir/$sampleName
