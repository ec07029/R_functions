# author: K Chu

# Get gene symbol from genomic coordinates

# Create grange object
miclip.with.distance.gr <- miclip.df.with.distance %>% select(c("seqnames", "start", "end", "strand")) %>% makeGRangesFromDataFrame

# Get genes
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotate)

genes.human <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Find coordinate overlaps between grange object and genomic database
overlaps <- findOverlaps(miclip.with.distance.gr, genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

# Obtain entrez id from overlaps. Annotate into gene symbols
miclip.gene.symbol <- lookUp(genes.human[subjectHits(overlaps),]$gene_id, 'org.Hs.eg.db', 'SYMBOL')
miclip.gene.symbol.final <- as.vector(unlist(miclip.gene.symbol))

# Add entrez id to the rows with overlap (aka. queryHits)
miclip.df.with.distance$entrez.id <- NA
miclip.df.with.distance$entrez.id[queryHits(miclip.with.distance.and.gene.gr)] <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)[subjectHits((miclip.with.distance.and.gene.gr))]$gene_id

# Add gene
miclip.df.with.distance$gene <- NA
miclip.df.with.distance$gene[queryHits(overlaps)] <- miclip.gene.symbol.final
