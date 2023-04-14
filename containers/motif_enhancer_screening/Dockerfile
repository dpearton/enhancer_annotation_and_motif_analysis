FROM  rocker/tidyverse:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing tidyverse and getopt"


# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      libbz2-dev \
      liblzma-dev

# Install scHelper
RUN   R -e "install.packages('optparse')" && \
      R -e "install.packages('ggplot2')" && \
      R -e "BiocManager::install(c('GenomicRanges', 'Rhtslib', 'Rsamtools', 'GenomicAlignments', 'rtracklayer'))"