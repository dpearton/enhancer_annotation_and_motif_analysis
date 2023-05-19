#############################################################################################################################
## Collapse gene name and annotations into single rows after bedtools intersect
#############################################################################################################################
library(tidyverse)
library(optparse)

# Read in command line opts
option_list <- list(
    make_option(c("-p", "--peak_bed"), action = "store", type = "character", help = "Path to peak bed file")
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

dat <- read.table(opt$peak_bed)

dat <- dat[,c(1:4,ncol(dat)-1:ncol(dat))]

colnames(dat) <- c('chrom', 'start', 'end', 'peakid', 'gene_id', 'gene_name')

# Concatenate gene names
dat <- dat %>%
distinct() %>% 
group_by(peakid) %>% 
reframe(chrom=chrom, start=start, end=end, peakid=peakid, gene_id = paste0(gene_id, collapse="|"), gene_name = paste0(gene_name, collapse="|")) %>%
distinct() %>%
relocate(peakid, .after=end)

# Save both annotated peaks and peak bed for motif screening
write.table(dat[,1:4], 'annotated_peaks.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(dat, 'annotated_peaks.tsv', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")