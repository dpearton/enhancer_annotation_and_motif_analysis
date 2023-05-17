
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
    make_option(c("-w", "--window"), action = "store", type = "integer", help = "Window surrounding TSS for which to screen for peaks", default = 100000))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

peaks <- read.csv(opt$peak_bed, stringsAsFactors = FALSE, sep = "\t", header = FALSE)[,1:4]
# peaks <- read.csv('/home/rstudio/test_data/test_lamprey_atac_peaks.bed', stringsAsFactors = FALSE, sep = "\t", header = FALSE)[,1:4]

colnames(peaks) <- c('chrom', 'start', 'end', 'peakid')

# keep only unique peaks (incase bed contains duplicated peaks - keep first hit)
peaks <- distinct(peaks, chrom, start, end, .keep_all = TRUE)

gtf_window <- import(opt$gtf)
# gtf_window <- import('/home/rstudio/results/filter_gtf_gene_list/pmz.genes_transcript_extract_filtered.gtf')


# Add +- window (default 100kb) to GTF coordinates - make sure that strand is considered when extracting start of gene - this is why we use resize
start(gtf_window) <- start(resize(gtf_window, width = 1))-opt$window
end(gtf_window) <- start(resize(gtf_window, width = 1))+opt$window

# Set negative coordinates to 0
start(gtf_window)[start(gtf_window) < 0] <- 0

# Create peak GRanges
peaks_granges <- GRanges(seqnames=peaks$chrom, ranges=IRanges(start=peaks$start, end=peaks$end, names=peaks$peakid))

# For each peak, identify overlapping peaks. No need to use strand info as we are looking for any overlap from any direction
peak_gtf_hits <- lapply(1:length(peaks_granges), function(x) findOverlaps(query = gtf_window, subject = peaks_granges[x,], ignore.strand=TRUE))
names(peak_gtf_hits) <- names(peaks_granges)

peak_gtf_hits <- lapply(peak_gtf_hits, function(row) {gtf_window[queryHits(row), c(opt$gene_id_col, opt$gene_name_col)]})

# Filter peaks which have at least one gtf hit
peak_gtf_hits <- peak_gtf_hits[unlist(lapply(peak_gtf_hits, length)) != 0]

# extract gene ids and gene names for each peak
peak_gtf_hits <-lapply(peak_gtf_hits, function(x) {
  out <- as.data.frame(x)[,c(opt$gene_id_col, opt$gene_name_col)]
  gene_id = paste(out[,opt$gene_id_col], collapse = "|")
  gene_name = paste(out[,opt$gene_name_col], collapse = "|")
  return(c(gene_id, gene_name))
  })

# Filter peaks which are annotated to at least one gene in the GTF
peaks <- peaks[match(names(peak_gtf_hits),peaks$peakid),]

# Make a dataframe of the peak gtf hits
peak_gtf_hits <- data.frame(peak_gtf_hits, row.names = c(opt$gene_id_col, opt$gene_name_col)) %>% t() %>% as.data.frame()

# Add annotation columns to peak file
annotated_peaks <- merge(peaks, peak_gtf_hits, by.x = 'peakid', by.y = 'row.names')

# Rearrange collumns
annotated_peaks <- annotated_peaks[, c('chrom', 'start', 'end', 'peakid', opt$gene_id_col, opt$gene_name_col)]

# Save both annotated peaks and peak bed for motif screening
write.table(annotated_peaks[,1:4], 'annotated_peaks.bed', row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
write.table(annotated_peaks, 'annotated_peaks.tsv', row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")