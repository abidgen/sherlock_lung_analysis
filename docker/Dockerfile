FROM rocker/r-ver:4.3.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    libglpk-dev \
    python3-pip \
    python3-dev \
    wget curl git \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    SigProfilerExtractor==1.1.23 \
    SigProfilerMatrixGenerator==1.2.17 \
    SigProfilerAssignment>=0.1.0

RUN R -e "install.packages(c('remotes','BiocManager'), \
    repos='https://cran.rstudio.com/', quiet=TRUE)"

RUN R -e "BiocManager::install(version='3.17', ask=FALSE, update=FALSE)"

RUN R -e "remotes::install_version('data.table', version='1.14.8', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('dplyr', version='1.1.3', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('tidyr', version='1.3.0', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('ggplot2', version='3.4.3', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('patchwork', version='1.1.3', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('pheatmap', version='1.0.12', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('RColorBrewer', version='1.1-3', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('viridis', version='0.6.4', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('ggpubr', version='0.6.0', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('ggrepel', version='0.9.3', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('cowplot', version='1.1.1', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('Seurat', version='4.3.0', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('harmony', version='0.1.0', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('reticulate', version='1.31', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "BiocManager::install(c( \
    'maftools', \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'MutationalPatterns', \
    'scDblFinder', \
    'SingleCellExperiment', \
    'BiocParallel' \
    ), version='3.17', update=FALSE, ask=FALSE)"

RUN R -e "cat('Seurat:', as.character(packageVersion('Seurat')), '\n'); \
          cat('maftools:', as.character(packageVersion('maftools')), '\n'); \
          cat('ggplot2:', as.character(packageVersion('ggplot2')), '\n')"

ENV OMP_NUM_THREADS=16
WORKDIR /analysis
CMD ["R"]

# Additional R packages for analysis
RUN R -e "BiocManager::install(c( \
    'limma', \
    'edgeR', \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'enrichplot', \
    'DOSE' \
    ), version='3.17', update=FALSE, ask=FALSE)"

RUN R -e "remotes::install_version('survival', version='3.5-7', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "remotes::install_version('survminer', version='0.4.9', \
    repos='https://cran.rstudio.com/', upgrade='never')"

RUN R -e "install.packages('estimate', \
    repos='http://r-forge.r-project.org', \
    dependencies=TRUE)"
