# obaDIA
obaDIA: one-step biological analysis pipeline for data-independent acquisition and other quantitative proteomics data

obaDIA takes a fragment-level, peptide-level or protein-level abundance matrix file from data-independent acquisition (DIA) mass spectrometry experiment,and a FASTA fromat protein sequence file as inputs, performs differential protein expression analysis, functional annotation and enrichment analysis in a completely automate way. obaDIA was designed for data-independent acquisition quantitative proteiomics data initially, it can also be applied to protein-level data produced by other quantitative proteomic techniques, such as iTRAQ/TMT, Label-free proteomics. obaDIA is easy to use and runs fast. All source codes and example data of obaDIA are distributed for academic use only. For any other use, including any commercial use, please contact us first (info@e-omics.com).

## Contents
* [Installation](#installation)
    * [Download](#1-download)
    * [Install](#2-install)
	* [Install database](#3-install-database)
    * [Install dependencies](#4-install-dependencies)
* [Quick start](#quick-start)
	* [From local](#1-from-local)
	* [From docker](#2-from-docker)
* [Usage](#usage)
    * [Input files](#1-input-files)
    * [Edit one of the example scripts](#2-edit-one-of-the-example-scripts)
    * [Run the modified script](#3-run-the-modified-script)
* [Parameters](#parameters)
* [Output files](#output-files)
* [Citation](#citation)

## Installation

### 1. Download
```Bash
git clone https://github.com/yjthu/obaDIA.git
```

### 2. Install

Installation of obaDIA is quite easy, just like this:
```Bash
cd obaDIA
bash INSTALL.sh
```

Alternatively, if you can not get an academic license for [signalP](https://services.healthtech.dtu.dk/service.php?SignalP-4.1), install with `-n` parameter:
```Bash
cd obaDIA
bash INSTALL.sh -n
```

### 3. Install database
Download and build database:
```Bash
cd db
bash build_db.sh
```
This step may take several minites to hours, depend on your network speed.

Then, download species data from Kobas 3.0 database [(ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD)](ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD]), both the files from `seq_pep` and `sqlite3` directories are needed. Abbreviation of sepcies refer to `db/species_abbr.txt`.

Finally, set the home directory of Kobas 3.0 database in `env/.kobasrc`, and copy this file to your home dirctory, like this:
```Bash
cp env/.kobasrc ~
```

### 4. Install dependencies
We have build a docker image with all dependencies installed for obaDIA, if you are a docker user, just get the image through command:
```Bash
docker pull yjcbscau/eomics-base:latest
```

Since an academic license is necessary for download and installation of [signalP 4.1](https://services.healthtech.dtu.dk/service.php?SignalP-4.1), you need to install it independently and set the enviroment variable correctly in `env/.bashrc` file. This step can be skipped when you build obaDIA with `bash INSTALL.sh -n` command.


Alternatively, if you want to install all dependencies from the beginning, please follow these steps: 

#### 4.1. Install conda and activate a new conda environment
Dowload and install a version of minicoda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

Then, add bioconda channels, creat a new python2 evironment and activate it
```Bash
conda create -n obadia python=2
. activate obadia
```

#### 4.2. Install R packages
```Bash
conda install R==3.6 r-xml
```

Then, enter R environment and install required packages
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

#### 4.3. Install Kobas 3.0
```Bash
conda config --add channels bioconda
conda install kobas
```

#### 4.4. Install perl5 and python2 libraries
```Bash
cpanm GO::Parser DBI DBD::SQLite
pip install eventlet bs4 django==1.8.3
```

#### 4.5. Check the build-in dependencies

Binary executable program of [`hmmpress`](http://hmmer.org/download.html), [`hmmscan`](http://hmmer.org/download.html), [`diamond`](https://github.com/bbuchfink/diamond), [`mapDIA`](https://sourceforge.net/projects/mapdia) and a slightly modified version of [`Trinotate-v3.1.0`](https://github.com/Trinotate/Trinotate/releases/tag/Trinotate-v3.1.0) software (named `Trinotate-v3.1.0-pro`) have already in the `src` folder along with the obaDIA release, make sure they can work correctly in your system. Otherwise, reinstall them through conda or from source.

Finally, copy the enviroment variable for obaDIA in `env/.bashrc` file to your home `~/.bashrc` file


## Quick start

#### 1. From local
We provide four example scripts and a set of test data in `example` directory. You can use the example script `example_fast.sh` for a quick test of obaDIA, it shoud be finished within three minutes. The basic command line is like this:

```Bash
cd example
nohup sh example_fast.sh &
```

#### 2. From docker
If you use docker to run obaDIA, you should creat a docker container first. The command line is like this:
```Bash
# 1. set directories
workdir=/storage/data/PROJECT/biouser1/TestDocker
obadir=/storage/data/PROJECT/biouser1/TestDocker/obaDIA
kobasdbdir=/storage/data/PUBLIC/databases/KOBAS_3.0_db
signalpdir=/storage/data/PUBLIC/softwares/SignalP/signalp-4.1

# 2. creat a docker container named obadia
docker run -v $workdir:$workdir -v $obadir:$obadir -v $obadir/env:/root -v $kobasdbdir:$kobasdbdir -v $signalpdir:$signalpdir --name obadia yjcbscau/eomics-base /bin/bash
docker run -it --privileged=true --volumes-from obadia yjcbscau/eomics-base /bin/bash

# 3. run script in the container
cd example
nohup sh example_fast.sh &
```

If you build obaDIA with `bash INSTALL.sh -n` command, `signalpdir` and `-v $signalpdir:$signalpdir` is no more needed. Then, the command line is like this:
```Bash
# 1. set directories
workdir=/storage/data/PROJECT/biouser1/TestDocker
obadir=/storage/data/PROJECT/biouser1/TestDocker/obaDIA
kobasdbdir=/storage/data/PUBLIC/databases/KOBAS_3.0_db

# 2. creat a docker container named obadia
docker run -v $workdir:$workdir -v $obadir:$obadir -v $obadir/env:/root -v $kobasdbdir:$kobasdbdir --name obadia yjcbscau/eomics-base /bin/bash
docker run -it --privileged=true --volumes-from obadia yjcbscau/eomics-base /bin/bash

# 3. run script in the container
cd example
nohup sh example_fast.sh &
```

For your own data, please copy the example script to your work dirctory and edit it refer to the `Usage`.


## Usage

#### 1. Input files
A protein sequence file in FASTA fromat and an abundance matrix file in TSV format are needed. These files can be generated by upstream proteins discovery and abundance estimation softwares for DIA or other quantitative proteomics experiment data. 

For DIA data, three types of abundance matrix file are accepted, including fragment-level, peptide-level and protein-level. For data generated by other techniques, such as DDA, only protein-level abundance matrix is accepted.

There are four example files in `example` directory, they are:
```Bash
seq.fa			# protein sequence file in FASTA format. The protein ID must match the abundance file, uniprot ID is prefered.
ab.tsv			# protein-level abundance matrix file in TSV format. The first column is the protein id, the other columns are the abundance for each sample.
pep.ab.tsv		# peptide-level abundance matrix file in TSV format. The first column is the protein id, the second columns is the peptide sequence, the other columns are the abundance for each sample.
frag.ab.tsv		# fragment-level abundance matrix file in TSV format. The first column is the protein id, the second columns is the peptide sequence, the third columns is the fragment information, the other columns are the abundance for each sample.
```

#### 2. Edit one of the example scripts

We provide four example script files in `example` directory, they are:
```Bash
example_prot.sh		# A standard script to run obaDIA using a protein-level matrix as input.
example_pep.sh		# A standard script to run obaDIA using a peptide-level matrix as input.
example_frag.sh		# A standard script to run obaDIA using a fragment-level matrix as input.
example_fast.sh		# Another script to run obaDIA using a protein-level matrix as input, but use the fast mode to skip several time-consuming steps.
```

What you need to do is copy one of these files to your work dirctory and edit it according to your own needs.
Now, take the `example_prot.sh` file as an example, to explain how the script works:
```Bash

# get the current work directory
dir=`pwd`			

# set the directory of input and output files, must be absolute directory 
fa=$dir/seq.fa
ab=$dir/ab.tsv
out=$dir/outdir

# run the main program of obaDIA to generate a workflow
# the parameters will be explaind in detail in the next section 
perl ../bin/oba.pl \
-fas $fa \
-exp $ab \
-out $out \
-group A1/A2/A3,B1/B2/B3,C1/C2/C3,D1/D2/D3,E1/E2/E3,F1/F2/F3 \
-name A,B,C,D,E,F \
-comp 2:1,3:1,4:1,5:1,6:1 \
-fc 0.5 \
-fdr 0.2 \
-spe hsa \
-thread 40

# then, an all-in-one workflow script named 'OneStep.sh' will be generated in the output dirctory
# run this workflow
sh $out/OneStep.sh >$out/OneStep.sh.o 2>&1
```

#### 3. Run the modified script
Once the modification of `example_prot.sh` is complete, you can run it locally or submit it to the cluster
```Bash
nohup sh example_prot.sh &

```

## Parameters

The parameters and description of obaDIA can be get through the command:
```Bash
perl ../bin/oba.pl
```

All parameters and short descriptions will be shown like this:
```Bash
Required:
		-fas            <s> : protein squence file, fasta format.
		-exp            <s> : protein abundance file, tsv format.
								protein-level, peptide-level or fragment-level file is acceptable
		-out            <s> : output directory. [An absolute path is required!]
		-group          <s> : group that the samples belong to. [eg: s1/s2/s3,s4/s5/s6]
								groups should be seperated by ',',
								samples within a group should be seperated by '/'.
		-name           <s> : the name of each group. [eg: group1,group2]
								should be seperated by ',' and ordered the same as '-group'.
		-spe            <s> : KOBAS3.0 KEGG enrichment species abbr. [refer to "db/species_abbr.txt" file]
								dowload KOBAS database from 'ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD'
Optional:
		-fc             <f> : log2(foldChange) cutoff for DE proteins. [default: 1]
		-fdr            <f> : FDR cutoff for DE proteins by mapDIA. [default: 0.1]
		-comp           <s> : how to compare the groups, use indexs of the group, [default: '1:2']
								required for multiple comparisons
		-alt            <s> : KOBAS3.0 GO/Rectome enrichment species abbr. [default: the same as '-spe']
		-mod            <s> : background choose for enrichment, total(T) or expressed(E). [default: E]
		-thread         <s> : thread number for hmmscan to perform Pfam annotation. [default: 1]
		-level          <s> : level of abundance matrix. Choice are prot/pep/frag [default: prot]
								required for peptide-level or fragment-level abundance
		-fast               :  fast mode to run the pipeline. a part of time-consuming analysis will be skipped.
		-h|?|help           :  Show this help
```

Each parameter is explained in detail as follows:
##### -fas
This parameter is used to set the input protein sequence file. Only standard FASTA format file is accepted, with the first line is the protein identifier start with `>`, the second line is the protein sequence, the protein sequence can also be broken into multiple lines. The protein identifier must be consist of letters and numbers, and no blank or special characters are allowed, except for `_`. For the the protein identifier uniprot ID is prefered because we need to link it to the uniprotKB database, but this is not necessary.

##### -exp
This parameter is used to set the input abundance matrix file. Only standard TSV format file is accepted, which means each column is seperated by tab. For DIA data, protein-level, peptide-level or fragment-level matrix is acceptable, if use protein-level file as input, no need to set the `-level` parameter; if use peptide-level or fragment-level matrix as input, the `-level`  parameter must be set to either `pep` or `frag`.

##### -level
This parameter is used to set the level or type of input abundance matrix file. The default value is `prot`; if use peptide-level or fragment-level matrix as input, it must be set to either `pep` or `frag`.

##### -out
This parameter is used to set the output directory, it must be a absolute path, a reletive path is not allowed.

##### -group
This parameter is used to assign the samples to the groups and set the order of groups. Use sample names but not group names to set it. For samples in the same group, the sample names should be seperated by `/`. For different groups, the seperater between samples should be `,`. The sample names must can be found in the first row of the abundance matrix, and must be consist of letters and numbers, and no blank or special characters are allowed, except for `_`.

##### -name
This parameter is used to set the group names, according to the order of `-group` parameter.

##### -comp
This parameter is used to set the comparison between groups for DE analysis. For only one comparison, the default value is set to `1:2`; for multiple comparisons, it must be set according to the index of the order of `-name` parameter. For example, `-comp 2:1,3:1,4:1,5:1,6:1` means every group shoud be compared to the first group, with the first group as a control.

##### -fc
This parameter is used to set the log2(foldChange) cutoff for DE proteins. The default value is 1.

##### -fdr
This parameter is used to set the FDR(false discovery rate) cutoff generated by mapDIA for DE proteins. The default value is 0.1.

##### -spe
This parameter is used to set the abbreviation of species for KEGG enrichment analysis. At the present release of KOBAS 3.0, more than 4000 species are supported which can refer to the "db/species_abbr.txt" file; for species not in this list, a relative or model species should be selected.

##### -alt
This parameter is used to set the abbreviation of species for GO/Rectome enrichment analysis. At the present release of KOBAS 3.0, dozens of species are supported for Rectome enrichment analysis; default value is the same as `-spe`, but for a species not have the Rectome annotation, a relative or model species must be set by this parameter.

##### -mod
This parameter is used to set the backgroud for enrichment analysis. We proposed a new enrichment strategy in obaDIA, which use expressed protein as backgroud for the enrichment analysis; this is the default choice. If you want perform enrichment analysis use a traditional manner, which use the total proteins in a species as backgroud, set this parameter as `-mod T`.

##### -thread
This parameter is used to set the thread number for hmmscan software in annotation step. Because Pfam annotation is time-consuming, using multiple threads could dramatically accelerate the annotation step.

##### -fast
This parameter is used to skip several time-consuming steps in obaDIA. No value is needed. It can make the total time-consuming of obaDIA pipeline shortened from 1~2 hours to 3~5 minutes.

## Output files

The output files could be slightly different based on the parameters. A typical output files dirctory is as follow, with the main result files marked:
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
│   └── TrinotatePFAM.out
├── Diff                                               # Abundance and Differntial expression(DE) analysis directory
│   ├── compare.txt
│   ├── diff_results                                   # DE analysis output directory
│   │   ├── BvsA                                       # One directory for each comparison 
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
│       ├── DIA.count.xls                              # peak number statistic for each sample
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
│   ├── BvsA                                           # DE enrichment analysis directory, one directory for each comparison
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
└── OneStep.sh                                         # all-in-one workflow script

```


## Citation

If you use obaDIA in your work, please cite the publication as follows:

Yan J, Zhai H, Zhu L, Ding X. obaDIA: one-step biological analysis pipeline for quantitative proteomics data from data independent acquisition mass spectrometry. bioRxiv 2020.05.28.121020; doi: https://doi.org/10.1101/2020.05.28.121020
