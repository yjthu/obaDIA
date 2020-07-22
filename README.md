# obaDIA
obaDIA: one-step biological analysis pipeline for data-independent acquisition and other quantitative proteomics data

obaDIA takes a fragment-level, peptide-level or protein-level abundance matrix file from data-independent acquisition (DIA) mass spectrometry experiment,and a FASTA fromat sequence file as inputs, performs differential protein expression analysis, functional annotation and enrichment analysis in a completely automate way. obaDIA was designed for data-independent acquisition quantitative proteiomics data initially, it can also be applied to protein-level data produced by other quantitative proteomic techniques, such as iTRAQ/TMT, Label-free proteomics. obaDIA is easy to use and runs fast. All source codes and example data of obaDIA are distributed for academic use only. For any other use, including any commercial use, please contact us first (info@e-omics.com).

## Contents
* [Installation](#installation)
    * [Download](#1-download)
    * [Install](#2-install)
	* [Install database](#3-install-database)
    * [Install dependencies](#4-install-dependencies)
* [Quick start](#quick-start)
* [Usage](#usage)
    * [Check the input files](#1-check-the-input-files)
    * [Generate a pipeline](#2-generate-a-pipeline)
    * [Run the one-step pipeline](#3-run-the-one-step-pipeline)
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

Alternatively, if you can not get a academic license for [signalP](https://services.healthtech.dtu.dk/service.php?SignalP-4.1), install with `-n` parameter:
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
This step may takes several minites to hours, depend on your network speed.

Then, download Kobas 3.0 database from [ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD](ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD),  both the files from `seq_pep` and `sqlite3` directories are needed. Abbreviation of sepcies refer to `db/species_abbr.txt`.

Finally, set the home directory of Kobas 3.0 database in `env/.kobasrc`, and copy this file to your home dirctory, like this:
```Bash
cp env/.kobasrc ~
```

### 4. Install dependencies
We have build a docker image with all dependencies installed for obaDIA, if you are a docker user, just get the image through command:
```Bash
docker pull yjcbscau/eomics-base:latest
```

Since [signalP 4.1](https://services.healthtech.dtu.dk/service.php?SignalP-4.1) need an academic license for download and installation, you need to install it independently and set the enviroment variable correctly in `env/.bashrc` file. This step can be skipped when you build obaDIA with `bash INSTALL.sh -n` command.


Alternatively, if you want to install all dependencies from the beginning, please follow these steps: 

#### 4.1. Install conda and activate a new conda environment
Dowload and install a version of minicoda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

Then, add bioconda channels, creat a new python2 evironment and activate it
```Bash
conda create -n oba python=2
. activate oba
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

#### 4.5. Check wrapped dependencies

Binary executable program of [`hmmpress`](http://hmmer.org/download.html), [`hmmscan`](http://hmmer.org/download.html), [`diamond`](https://github.com/bbuchfink/diamond), [`mapDIA`](https://sourceforge.net/projects/mapdia) and a slightly modified version of [`Trinotate-v3.1.0`](https://github.com/Trinotate/Trinotate/releases/tag/Trinotate-v3.1.0) software (named `Trinotate-v3.1.0-pro`) have already in the `src` folder along with the obaDIA release, make sure they can work correctly in your system. Otherwise, reinstall them through conda or from source.

Finally, copy the enviroment variable for obaDIA in `env/.bashrc` file to your home `.bashrc` file


## Quick start

We provide four example scripts and a set of test data in `example` directory. You can use the example script `example_fast.sh` for a quick test of obaDIA, it shoud be finished within three minutes. The basic command line is like this:

```Bash
cd example
nohup sh example_fast.sh &
```

If you use docker to run obaDIA, you should creat a docker container first. The command line is like this:
```Bash
#1. set directories, edit them as your own directory
workdir=/storage/data/PROJECT/biouser1/TestDocker
obadir=/storage/data/PROJECT/biouser1/TestDocker/obaDIA
kobasdbdir=/storage/data/PUBLIC/databases/KOBAS_3.0_db
signalpdir=/storage/data/PUBLIC/softwares/SignalP/signalp-4.1

#2. creat a docker container named obadia
docker run -v $workdir:$workdir -v $obadir:$obadir -v $obadir/env:/root -v $kobasdbdir:$kobasdbdir -v $signalpdir:$signalpdir --name obadia yjcbscau/eomics-base /bin/bash
docker run -it --privileged=true --volumes-from obadia yjcbscau/eomics-base /bin/bash

#3. run script in the container
cd example
nohup sh example_fast.sh &
```

If you build obaDIA with `bash INSTALL.sh -n` command, `signalpdir` and `-v $signalpdir:$signalpdir` is no more needed.
For your own data, pleas copy the example script to your work dirctory and edit it first.

