
###########################################################################
## Compare DE data with DE TFs from Chen et al. (2017) Development
###########################################################################

library(openxlsx)
library(tidyverse)

# Download supplementary file 4 from Chen et al. (2017)
# PPR vs otic 5/6ss, PPR vs otic 8/9ss, PPR vs otic 11/12ss. Dataset is already filtered for transcription factors

temp <- tempfile()
download.file("http://www.biologists.com/DEV_Movies/DEV148494/TableS4.xlsx", temp)

# read xlsx file
otic_enr <- read.xlsx(temp, startRow = 15) 
unlink(temp)

# assign column names
colnames(otic_enr)[1:13] <- paste0(colnames(otic_enr)[1:13], c(rep('_normalised_count', 4), rep('_foldChange', 3),
                                                               rep('_pval', 3), rep('_padj', 3)))

# assign row names
rownames(otic_enr) <- otic_enr[,17]
otic_enr[,17] <- NULL

# remove genes from Chen dataset which are not DE in at least one of the stages (not sure why these are in the supplementary file)
otic_enr <- otic_enr[apply(otic_enr, 1, function(x) any(!is.na(x[c("5-6ss_foldChange", "8-9ss_foldChange", "11-12ss_foldChange")]))),]

gene_ids <- otic_enr %>% pull(ensembl_gene_id)

writeLines(gene_ids, './de_gene_ids.txt')