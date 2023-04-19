
#############################################################################################################################
## Annotate peaks within window of genes in provided gtf. Output a filtered annptated peak list and simplified peak bed file. 
#############################################################################################################################

library(GenomicRanges)
library(rtracklayer)
library(tidyverse) 
library(optparse)

# Read in command line opts
option_list <- list(
    make_option(c("-n", "--gene_name_col"), action = "store", type = "character", help = "Column name in gtf containing gene names", default = 'gene_name'),
    make_option(c("-i", "--gene_id_col"), action = "store", type = "character", help = "Column name in gtf containing gene ids", default = 'gene_id'),
    make_option(c("-p", "--peak_bed"), action = "store", type = "character", help = "Path to peak bed file"),
    make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Path to gtf file"),
    make_option(c("-c", "--ctcf_flanking_peaks"), action = "store", type = "character", help = "Path to ctcf_flanking_peaks bed file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

peaks <- read.csv(opt$peak_bed, stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'peakid'), sep = "\t")
# peaks <- read.csv('/temp/test_data/peaks.bed', stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'peakid'), sep = "\t")
# peaks <- read.csv('./results/extend_peaks/extended_peaks.bed', stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'peakid'), sep = "\t")

# keep only unique peaks (incase bed contains duplicated peaks - keep first hit)
peaks <- distinct(peaks, chrom, start, end, .keep_all = TRUE)

gtf <- import(opt$gtf)
# gtf <- import('./results/filter_gtf_gene_list/Gallus_gallus.GRCg6a.97_filtered.gtf')

# Use strand info to get start site of all genes - next we will identify which ctcf_flanking_peaks overlap with this start site
gtf <- resize(gtf, width = 1)

ctcf <- read.csv(opt$ctcf_flanking_peaks, stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'peakid'), sep = "\t")
# ctcf <- read.csv('./results/bedtools_sort_flanking_ctcf/sorted_ctcf_flanking_peaks.bed', stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'peakid'), sep = "\t")

# Create peak ctcf GRanges
ctcf_flanking_peaks_granges <- GRanges(seqnames=ctcf$chrom, ranges=IRanges(start=ctcf$start, end=ctcf$end, names=ctcf$peakid))

# For each ctcf_flanking_peaks, identify overlapping gene start sites
ctcf_flanking_peaks_gtf_hits <- lapply(1:length(ctcf_flanking_peaks_granges), function(x) findOverlaps(gtf, ctcf_flanking_peaks_granges[x,]))
names(ctcf_flanking_peaks_gtf_hits) <- names(ctcf_flanking_peaks_granges)

ctcf_flanking_peaks_gtf_hits <- lapply(ctcf_flanking_peaks_gtf_hits, function(row) {gtf[queryHits(row), c(opt$gene_id_col, opt$gene_name_col)]})

# Filter ctcf_flanking_peaks which have at least one gtf hit
ctcf_flanking_peaks_gtf_hits <- ctcf_flanking_peaks_gtf_hits[unlist(lapply(ctcf_flanking_peaks_gtf_hits, length)) != 0]

# extract gene ids and gene names for each ctcf_flanking_peaks
ctcf_flanking_peaks_gtf_hits <- lapply(ctcf_flanking_peaks_gtf_hits, function(x) {
  out <- as.data.frame(x)[,c(opt$gene_id_col, opt$gene_name_col)]
  gene_id = paste(out[,opt$gene_id_col], collapse = "|")
  gene_name = paste(out[,opt$gene_name_col], collapse = "|")
  return(c(gene_id, gene_name))
  })

# Filter peaks which are annotated to at least one gene in the GTF
peaks <- peaks[match(names(ctcf_flanking_peaks_gtf_hits),peaks$peakid),]

# Make a dataframe of the peak gtf hits
ctcf_flanking_peaks_gtf_hits <- data.frame(ctcf_flanking_peaks_gtf_hits, row.names = c(opt$gene_id_col, opt$gene_name_col)) %>% t() %>% as.data.frame()

# Add annotation columns to peak file
annotated_peaks <- merge(peaks, ctcf_flanking_peaks_gtf_hits, by.x = 'peakid', by.y = 'row.names')

# Rearrange collumns
annotated_peaks <-annotated_peaks[, c('chrom', 'start', 'end', 'peakid', opt$gene_id_col, opt$gene_name_col)]

# Save both annotated peaks and peak bed for motif screening
write.table(annotated_peaks[,1:4], 'annotated_peaks.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(annotated_peaks, 'annotated_peaks.tsv', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

