# obaDIA
obaDIA: one-step biological analysis pipeline for quantitative proteomics data from data independent acquisition mass spectrometry

obaDIA uses protein abundance and protein sequence data as inputs and performs differential protein expression analysis, functional annotation and enrichment analysis in a completely automate way. obaDIA was designed for quantitative data independent acquisition proteiomics data initially, it can also be applied to data produced by other quantitative proteomic techniques, such as iTRAQ/TMT, Label-free proteomics. obaDIA is not only easy to use but also runs fast. It usually takes only a few tens of minutes to complete a typical analysis workflow.

## Contents
* [Installation](#installation)
    * [Download](#1-download)
    * [Install](#2-install)
    * [Install dependencies](#3-install-dependencies)
* [Usage](#usage)
    * [Check the input files](#1-check-the-input-files)
    * [Generate a pipeline](#2-generate-a-pipeline)
    * [Run the one-step pipeline](#3-run-the-one-step-pipeline)
* [Parameters](#parameters)
* [Output files](#output-files)

## Installation

### 1. Download
```Bash
git clone https://github.com/yjthu/obaDIA.git
```

### 2. Install
```Bash
# install code
cd obaDIA
bash INSTALL.sh

# download and build database
cd db
bash build_db.sh
```
The step of build database may takes several minites to hours, depend on the network speed.


### 3. Install dependencies

I recommend to use conda for installing dependencies. Of course, you can install them all from the source, but I think that would be a tedious work.

#### 3.1. Install conda
Dowload and install a version of minicoda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

Then, add bioconda channels, creat a new python2 evironment and activate it
```Bash
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n diapp python=2
. activate diapp
```

#### 3.2. Install R packages
```Bash
conda install R==3.6

# then, enter R environment and install required packages
R
```
```R
# Install R packages
# 1. from CRAN
req.pcg <- function(pcg){
    new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
    if (length(new)) install.packages(new, dependencies = T)
    sapply(pcg, require, ch = T)
  }
all.pcg <- c("reshape2", "ggplot2", "corrplot", "pheatmap", "optparse", "PerformanceAnalytics")
req.pcg(all.pcg)

# 2. from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("goseq")
BiocManager::install("goTools")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")

```

#### 3.3. Install kobas 3.0
```Bash
conda install kobas
```

Download the database of kobas 3.0 from [http://kobas.cbi.pku.edu.cn/kobas3/download](http://kobas.cbi.pku.edu.cn/kobas3/download/) and configure the `.kobasrc` file correctly as described.


#### 3.4. Install perl5 and python2 libraries
```Bash
cpanm GO::Parser DBI DBD::SQLite

pip install eventlet bs4 django==1.8.3
```

#### 3.5. Install signalP 4.1

Download and install from [https://services.healthtech.dtu.dk/service.php?SignalP-4.1](https://services.healthtech.dtu.dk/service.php?SignalP-4.1) and export the environment variable in `.bashrc` file.

#### 3.6. Other dependencies

Binary executable program of [`hmmpress`](http://hmmer.org/download.html), [`hmmscan`](http://hmmer.org/download.html), [`diamond`](https://github.com/bbuchfink/diamond), [`sqlite3`](https://www.sqlite.org/download.html), [`mapDIA`](https://sourceforge.net/projects/mapdia) and a slightly modified version of [`Trinotate-v3.1.0`](https://github.com/Trinotate/Trinotate/releases/tag/Trinotate-v3.1.0) software (named `Trinotate-v3.1.0-pro`) have already in the `src` folder, make sure they can work correctly in your system. Otherwise, reinstall them through conda or from source.

Export the environment variable in `.bashrc` file
```Bash
oba_home=
export PATH=$oba_home/src:$oba_home/src/mapDIA:$oba_home/src/Trinotate-v3.1.0-pro:$PATH
```

## Usage

#### 1. Check the input files

If the input fasta file is not consistent with the abundance file, generate a pair of matched input files:
```Bash
# the raw input files are looks like this:
$ wc -l ab
569 ab

$ head -n 5 ab|cut -f 1-3
Uniprot A1      A2
P31946  9905933 7563083
P62258  8696484 6608905
P33176  3294717 2889897
O14654  2167646 1617881

$ grep -c '>' seq
571

$ head -n 2 seq
>sp|Q9ULT8|HECD1_HUMAN E3 ubiquitin-protein ligase HECTD1 OS=Homo sapiens OX=9606 GN=HECTD1 PE=1 SV=3
MADVDPDTLLEWLQMGQGDERDMQLIALEQLCMLLLMSDNVDRCFETCPPRTFLPALCKIFLDESAPDNVLEVTARAITYYLDVSAECTRRIVGVDGAIKA

# the above input files are not consistent with each other in terms of sequence id/number. 
# Use the example script to generate a pair of matched input files:

$ perl matchInputs.pl ab seq

# then, a pair of new files are generated, they looks like this:
$ wc -l ab.tsv
567 ab.tsv

$ grep -c '>' seq.fa
566

$ head -n 2 seq.fa
>Q9ULT8
MADVDPDTLLEWLQMGQGDERDMQLIALEQLCMLLLMSDNVDRCFETCPPRTFLPALCKIFLDESAPDNVLEVTARAITYYLDVSAECTRRIVGVDGAIKA

```

If you use a pair of input files not consistent with each other, an error with be thrown out when generating a pipeline.

#### 2. Generate a pipeline
```Bash
dir=`pwd`

fa=$dir/seq.fa
ab=$dir/ab.tsv
out=$dir/outdir

perl ../bin/oba.pl \
-fas $fa \
-exp $ab \
-out $out \
-group A1/A2/A3,B1/B2/B3,C1/C2/C3,D1/D2/D3,E1/E2/E3,F1/F2/F3 \
-name A,B,C,D,E,F \
-comp 2:1,3:1,4:1,5:1,6:1 \
-fc 1 \
-fdr 0.1 \
-spe hsa

```

#### 3. Run the one-step pipeline

A workflow script will be generated, run it locally or submit it to the cluster. For example:
```Bash
nohup sh $out/OneStep.sh &

```

## Parameters
```Bash

Required:
		-fas            <s> : protein squence file
		-exp            <s> : protein abundance file
		-out            <s> : output directory. [An absolute path is required!]
		-group          <s> : group that the samples belong to: s1/s2/s3,s4/s5/s6
		-name           <s> : group names: g1,g2
		-spe            <s> : kobas3.0 KEGG enrichment species abbr. [refer to "species_abbr.txt" file]

Optional:
		-fc             <f> : default: 1
		-fdr            <f> : default: 0.1
		-comp           <s> : how to compare the groups, just use numbers, 1:2
		-alt            <s> : kobas3.0 GO/Rectome enrichment species abbr. [default: the same as "-spe"]

		-h|?|help           :  Show this help

```


## Output files

```Bash
outdir/                                               
├── Anno                                               # Functional annotation output directory
│   ├── Annotation.txt
│   ├── Annotation.xls                                 # formated annotation summay table
│   ├── blastp.outfmt6
│   ├── eggNOG.anno
│   ├── eggNOG.anno.plotdata
│   ├── eggNOG.anno.plotdata.bar.pdf                   # eggNOG annotation bar plot
│   ├── eggNOG.anno.plotdata.bar.png
│   ├── GO.anno
│   ├── GO.anno.output
│   │   ├── GO_classification_bar.pdf                  # GO annotation bar plot
│   │   ├── GO_classification_bar.png
│   │   ├── GO_classification_count.txt
│   │   ├── GO_classification.R
│   │   └── GO_classification.xls                      # GO annotation detail
│   ├── Kegg.anno
│   ├── Kegg.anno.output
│   │   ├── KEGG_classification_count.txt
│   │   ├── KEGG_classification.pdf                    # KEGG annotation bar chart
│   │   ├── KEGG_classification.png
│   │   ├── KEGG_classification.R
│   │   └── KEGG_classification.xls                    # KEGG annotation detail
│   ├── Kegg.anno.plotdata
│   ├── PFAM.anno
│   ├── pfam.log
│   ├── seq.fa.gene_trans_map
│   ├── signalp.out
│   ├── trinotate_annotation_report.xls
│   ├── TrinotatePFAM.out
│   └── Trinotate.sqlite
├── Diff                                               # Abundance and Differntial expression(DE) analysis directory
│   ├── compare.txt
│   ├── diff_results                                   # DE analysis output directory
│   │   ├── BvsA                                       # One directory for each comparation 
│   │   │   ├── BvsA.diff.id
│   │   │   ├── BvsA.diff.seq
│   │   │   ├── BvsA.diff.updown
│   │   │   ├── BvsA.diff.xls                          # formated differntial expression information table 
│   │   │   ├── BvsA.diff.xls.volcano.pdf              # volcano plot
│   │   │   └── BvsA.diff.xls.volcano.png
│   │   ├── CvsA
│   │   │   ├── ...
│   │   ├── DvsA
│   │   │   ├── ...
│   │   ├── EvsA
│   │   │   ├── ...
│   │   └── FvsA
│   │       ├── ...
│   ├── group.txt
│   ├── mapDIA                                         # raw DE analysis output directory
│   │   ├── analysis_output.txt
│   │   ├── analysis_output_wide_format.txt
│   │   ├── fragment_selection.txt
│   │   ├── log2_data.txt
│   │   ├── mapDIA.conf
│   │   └── param.txt
│   ├── mapDIA.input
│   └── tableProcess                                   # Abundance analysis output directory
│       ├── DIA.count.xls                              # raw abundance table
│       ├── DIA.means.xls                              # mean abundance table
│       ├── DIA.means.xls.abHeatmap.pdf                # mean abundance heatmap
│       ├── DIA.means.xls.abHeatmap.png
│       ├── DIA.means.xls.boxplot.pdf                  # mean abundance box plot
│       ├── DIA.means.xls.boxplot.png
│       ├── DIA.means.xls.corHeatmap.pdf               # mean correlation heatmap
│       ├── DIA.means.xls.corHeatmap.png
│       ├── DIA.means.xls.corMatrix.pdf                # mean correlation-matrix plot
│       ├── DIA.means.xls.corMatrix.png
│       ├── DIA.means.xls.density.pdf                  # mean abundance density plot
│       ├── DIA.means.xls.density.png
│       ├── DIA.means.xls.violon.pdf                   # mean abundance violon plot
│       ├── DIA.means.xls.violon.pdf
│       ├── DIA.normalized.xls                         # normalized abundance table
│       ├── DIA.normalized.xls.abHeatmap.pdf           # normalized abundance heatmap
│       ├── DIA.normalized.xls.abHeatmap.png
│       ├── DIA.normalized.xls.boxplot.pdf             # normalized abundance box plot		
│       ├── DIA.normalized.xls.boxplot.png
│       ├── DIA.normalized.xls.corHeatmap.pdf          # normalized correlation heatmap
│       ├── DIA.normalized.xls.corHeatmap.png
│       ├── DIA.normalized.xls.corMatrix.pdf           # normalized correlation-matrix plot
│       ├── DIA.normalized.xls.corMatrix.png
│       ├── DIA.normalized.xls.density.pdf             # normalized abundance density plot
│       ├── DIA.normalized.xls.density.png
│       ├── DIA.normalized.xls.violon.pdf              # normalized abundance violon plot
│       └── DIA.normalized.xls.violon.png
├── Enrich                                             # Functional enrichment analysis directory
│   ├── All                                            # Expressed protein functional enrichment analysis directory
│   │   ├── All.go_bar.pdf                             # GO enrichment bar plot
│   │   ├── All.go_bar.png
│   │   ├── All.hsa.annotate
│   │   ├── All.hsa.dmout
│   │   ├── All.identify.GO
│   │   ├── All.identify.GO.xls                        # formated GO enrichment result table
│   │   ├── All.identify.Kegg
│   │   ├── All.identify.Rectome
│   │   ├── All.identify.Rectome.top20
│   │   ├── All.identify.Rectome.xls                   # formated Rectome enrichment result table
│   │   ├── All.KeggEnrich
│   │   ├── All.KeggEnrich.html                        # KEGG enrichment html formated report
│   │   ├── All.KeggEnrich.top20
│   │   ├── All.KeggEnrich.top20.scatterplot.pdf       # KEGG enrichment scatter plot
│   │   ├── All.KeggEnrich.top20.scatterplot.png
│   │   ├── All.KeggEnrich.xls                         # formated KEGG enrichment result table
│   │   ├── All.Rectome.Enrich.scatterplot.pdf         # Rectome enrichment scatter plot
│   │   ├── All.Rectome.Enrich.scatterplot.png
│   │   ├── go.resource
│   │   ├── seq.len
│   │   └── src
│   ├── BvsA                                           # DE enrichment analysis directory, one directory for each comparation
│   │   ├── BvsA.annomap
│   │   ├── BvsA.go_bar.pdf                            # GO enrichment bar plot
│   │   ├── BvsA.go_bar.png
│   │   ├── BvsA.GO_BP_DAG.pdf                         # GO enrichment DAG plots (BP)
│   │   ├── BvsA.GO_BP_DAG.png
│   │   ├── BvsA.GO_CC_DAG.pdf                         # GO enrichment DAG plots (CC)
│   │   ├── BvsA.GO_CC_DAG.png
│   │   ├── BvsA.GO_MF_DAG.pdf                         # GO enrichment DAG plots (MF)
│   │   ├── BvsA.GO_MF_DAG.png
│   │   ├── BvsA.hsa.annotate
│   │   ├── BvsA.hsa.dmout
│   │   ├── BvsA.identify.GO
│   │   ├── BvsA.identify.GO.xls                       # formated GO enrichment result table
│   │   ├── BvsA.identify.Kegg
│   │   ├── BvsA.identify.pfam
│   │   ├── BvsA.identify.pfam.xls                     # formated PFAM enrichment result table
│   │   ├── BvsA.identify.pfam.xls.top20
│   │   ├── BvsA.identify.Rectome
│   │   ├── BvsA.identify.Rectome.top20
│   │   ├── BvsA.identify.Rectome.xls                  # formated Rectome enrichment result table
│   │   ├── BvsA.KeggEnrich
│   │   ├── BvsA.KeggEnrich.html
│   │   ├── BvsA.KeggEnrich.top20
│   │   ├── BvsA.KeggEnrich.top20.scatterplot.pdf      # KEGG enrichment scatter plot
│   │   ├── BvsA.KeggEnrich.top20.scatterplot.png
│   │   ├── BvsA.KeggEnrich.xls                        # formated KEGG enrichment result table
│   │   ├── BvsA.PFAM.Enrich.scatterplot.pdf           # PFAM enrichment scatter plot
│   │   ├── BvsA.PFAM.Enrich.scatterplot.png
│   │   ├── BvsA.Rectome.Enrich.scatterplot.pdf        # Rectome enrichment scatter plot
│   │   ├── BvsA.Rectome.Enrich.scatterplot.png
│   │   ├── gene2pfam_formated.txt
│   │   ├── go.resource
│   │   └── src
│   ├── CvsA
│   │   ├── ...
│   ├── DvsA
│   │   ├── ...
│   ├── EvsA
│   │   ├── ...
│   ├── FvsA
│   │   ├── ...
└── OneStep.sh                                         # one-step workflow script

```
