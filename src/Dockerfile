FROM centos:centos7
MAINTAINER e-Omics Development Team <info@e-omics.com>

RUN yum -y update \
    && yum -y install curl bzip2 \
    && yum -y install vim-enhanced gcc perl perl-App-cpanminus \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ \
    && conda update -n base -c defaults conda \
    && conda install -y R==3.6 \
    && conda install -y kobas \
    && conda install -y r-xml \
    && Rscript -e 'install.packages("reshape2",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("ggplot2",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("corrplot",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("pheatmap",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("optparse",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("PerformanceAnalytics",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'install.packages("BiocManager",dependencies=TRUE,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")' \
    && Rscript -e 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") ; BiocManager::install("GO.db")' \
    && Rscript -e 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") ; BiocManager::install("goseq")' \
    && Rscript -e 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") ; BiocManager::install("goTools")' \
    && Rscript -e 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") ; BiocManager::install("topGO")' \
    && Rscript -e 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") ; BiocManager::install("Rgraphviz")' \
    && pip install eventlet bs4 django==1.8.3 \
    && conda update conda \
    && conda clean --all --yes \
    && rpm -e --nodeps bzip2 \
    && yum -y autoremove \
    && yum clean all \
    && cpanm install GO::Parser DBI DBD::SQLite \
    && rm -fr root/.cpanm

