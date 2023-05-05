###########################################################################
## Extract DEGs from Chen et al. (2017) Development
###########################################################################

# Download supplementary file 3 from Chen et al. (2017)
# PPR vs otic 5/6ss, PPR vs otic 8/9ss, PPR vs otic 11/12ss.
# install.packages('openxlsx')
library(openxlsx)
library(biomaRt)

# read xlsx file
otic_enr <- read.xlsx('./TableS2.xlsx', startRow = 15)

# DE gene ids
gene_ids <- otic_enr$ensembl_gene_id
gene_ids <- gene_ids[!gene_ids == "-"]


# Get biomart GO annotations for TFs
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'ggallus_gene_ensembl',
                      version = 80)

biomart_out <- getBM(attributes=c('chromosome_name','start_position','end_position', "external_gene_name", 'ensembl_gene_id'),
        filters = 'ensembl_gene_id',
        values = gene_ids,
        mart = ensembl)

# replace missing gene names with gene ids
biomart_out$external_gene_name[biomart_out$external_gene_name == ""] <- biomart_out$ensembl_gene_id[biomart_out$external_gene_name == ""]

biomart_out$chromosome_name <- paste0('chr', biomart_out$chromosome_name)

write.table(biomart_out[,c('chromosome_name','start_position','end_position', "external_gene_name")], './temp.bed', quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)


# read in output from ucsc lift over
lift_over <- read.csv('~/Downloads/hglft_genome_2b32d_51c3d0.bed', sep ="\t", col.names = c('chromosome_name','start_position','end_position', "external_gene_name"))


# znf385c = ENSGALG00000003403

goi <- c('PICK1','SPRY2','ENSGALG00000010722','HEY1','TSPAN13','MAP7','Blimp-1','RNF150','SPRY1','ZNF330','PHLDA2','STK39','HS6ST1','PLSCR1','MYO1E','MMD','AUTS2','EYA2','ZBTB16','ENSGALG00000003403','ZNF462')

# Check which of our key genes are not in the lift over - spry2 and plscr1 are missing because they are not in Chen et al. 2017 initial DE list
goi[!goi %in% lift_over[[4]]]


lift_over[[1]] <- sub("chr", "", lift_over[[1]])

lift_over <- dplyr::arrange(lift_over, chromosome_name, start_position)

# Get biomart ids for new genome version
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'ggallus_gene_ensembl',
                      version = 97)

biomart_out <- getBM(attributes=c('chromosome_name','start_position','end_position', "external_gene_name", 'ensembl_gene_id'),
                     filters = c("chromosome_name", "start","end"),
                     values = list("chromosome_name" = lift_over$chromosome_name, "start" = lift_over$start_position, "end" = lift_over$end_position),
                     mart = ensembl)


temp <- lapply(1:nrow(lift_over), function(x) {
  getBM(attributes=c('chromosome_name','start_position','end_position', "external_gene_name", 'ensembl_gene_id'),
        filters = c("chromosome_name", "start","end"),
        values = list("chromosome_name" = lift_over[x,1], "start" = lift_over[x,2], "end" = lift_over[x,3]),
        mart = ensembl)
})


write.csv(temp, file = './galgal6_chen_degs.csv', quote = FALSE, row.names = FALSE, col.names = TRUE)









# 
# 
# biomart_out[biomart_out$external_gene_name %in% goi,]
# 
# 
# biomart_out[biomart_out$external_gene_name %in% 'HEY1',]
# 
# temp[temp$external_gene_id %in% goi,]
# 
# 
# missing_genes <- goi[!goi %in% biomart_out$external_gene_name]
# 
# 
# dplyr::filter(temp, external_gene_id %in% missing_genes)
# 
# temp[temp$external_gene_id %in% missing_genes,]
