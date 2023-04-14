FROM  rocker/tidyverse:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing tidyverse and getopt"

# Install scHelper
RUN   R -e "install.packages('optparse')" && \
      R -e "install.packages('ggplot2')" && \
      R -e "BiocManager::install('GenomicRanges')" && \
      R -e "BiocManager::install('rtracklayer')"

