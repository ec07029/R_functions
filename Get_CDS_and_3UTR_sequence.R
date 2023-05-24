# Obtain CDS and 3'UTR sequence


library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# CDS
results.cds <- getBM(attributes = c('entrezgene_id',
                                    'coding'),
                     filters = 'entrezgene_id',
                     values = '11215',
                     mart = ensembl)

# 3'UTR
results <- getBM(attributes = c('entrezgene_id',
                                '3utr'),
                 filters = 'entrezgene_id',
                 values = unique(molm.filter$entrez.id),
                 mart = ensembl)

write.csv(results, "3UTR_sequences.csv", row.names = FALSE)