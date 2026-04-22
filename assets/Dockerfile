# Use a base image from rocker with RStudio
FROM rocker/rstudio:4.4.2
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Author:	Austyn Trull

# LABEL adds metadata to an image:
LABEL maintainer="Austyn Trull <atrull@uab.edu>"
LABEL description="Environment which contains tools used by the nf_xpatial pipeline"
LABEL org.opencontainers.image.title="nf_xpatial"
LABEL org.opencontainers.image.description="Environment which contains tools used by the nf_xpatial pipeline"
LABEL org.opencontainers.image.version="0.0.4"
LABEL org.opencontainers.image.authors="Austyn Trull <atrull@uab.edu>"
LABEL org.opencontainers.image.category="Xenium"

# pals, tools

# Installing CRAN with install.packages as install.packages()
RUN R --no-restore --no-save -e "install.packages(c('arrow', 'av', 'data.table', 'dplyr', 'DT', 'geometry', 'ggExtra', 'ggplot2', 'ggplotify'))"
RUN R --no-restore --no-save -e "install.packages(c('jsonlite', 'knitr', 'leidenbase', 'optparse', 'patchwork', 'pdftools', 'Polychrome', 'progressr'))"
RUN R --no-restore --no-save -e "install.packages(c('purrr', 'RColorBrewer', 'remotes', 'R.utils', 'scales', 'SpatialExperiment'))"
RUN R --no-restore --no-save -e "install.packages(c('stringr', 'SummarizedExperiment', 'tidyr', 'tidyverse'))"

# Install github packages
RUN R --no-restore --no-save -e "remotes::install_github('prabhakarlab/Banksy')"

# Install Harmony from immunogenomics repository
RUN R --no-restore --no-save -e "devtools::install_github('immunogenomics/harmony', upgrade='never')"

# Install Seurat and its data
RUN R --no-restore --no-save -e "install.packages(c('Seurat'))"
RUN R --no-restore --no-save -e "remotes::install_github('satijalab/seurat-wrappers')"

# Set automatically when building with --platform
ARG TARGETPLATFORM
ENV TARGETPLATFORM=${TARGETPLATFORM:-linux/amd64}

## Set env variables
ENV PLATFORM=${TARGETPLATFORM}

