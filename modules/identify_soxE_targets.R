
###########################################################################
## Compare DE data with DE TFs from Chen et al. (2017) Development
###########################################################################

setwd('/home/rstudio/')

library(tidyverse)
library(stringr)
motif_matches <- read.csv('./fimo_output/fimo.tsv', sep='\t', comment.char = "#", stringsAsFactors=FALSE)

motif_matches <- motif_matches %>% separate(sequence_name, into=c("PeakID", "coordinates"), sep="::", convert = TRUE)

# Load putative enhancers from Buzzi 2022
putative_enhancers <- read.delim('./putative_enhancers_annotated.txt', stringsAsFactors=FALSE)
colnames(putative_enhancers)[1] <- 'PeakID'


# Extract gene id for peaks with motif matches
motif_matches[['gene_id']] <- putative_enhancers[match(motif_matches$PeakID, putative_enhancers$PeakID), 'Entrez.ID']
motif_matches[['gene_name']] <- putative_enhancers[match(motif_matches$PeakID, putative_enhancers$PeakID), 'Gene.Name']

motif_matches <- motif_matches %>% arrange(motif_alt_id)

write.table(motif_matches, file = './otic_genes_with_motif_match.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')