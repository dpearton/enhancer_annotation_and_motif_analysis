
###########################################################################
## Annotate FIMO motif hits using annotated enhancer file
###########################################################################

install.packages('getopt')

library(tidyverse)
library(stringr)
library(getopt)

spec = matrix(c(
  'gene_name_col', 'l', 2, "character",
  'gene_id_col', 'c', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

motif_matches <- read.csv(list.files('./', full.names = TRUE, pattern = '*.tsv'), sep='\t', comment.char = "#", stringsAsFactors=FALSE)

motif_matches <- motif_matches %>% separate(sequence_name, into=c("PeakID", "coordinates"), sep="::", convert = TRUE)

# Load putative enhancers from Buzzi 2022
putative_enhancers <- read.delim(list.files('./', full.names = TRUE, pattern = '*.txt'), stringsAsFactors=FALSE)
colnames(putative_enhancers)[1] <- 'PeakID'

# Extract gene id for peaks with motif matches
motif_matches[['gene_id']] <- putative_enhancers[match(motif_matches$PeakID, putative_enhancers$PeakID), opt$gene_id_col]
motif_matches[['gene_name']] <- putative_enhancers[match(motif_matches$PeakID, putative_enhancers$PeakID), opt$gene_name_col]

motif_matches <- motif_matches %>% arrange(motif_alt_id)

write.table(motif_matches, file = './annotated_motif_hits.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')