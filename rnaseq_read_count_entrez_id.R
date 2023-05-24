# author: K Chu

library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx) # xlsx package completely unusable

register(MulticoreParam(workers = 6))

bam.files <- Sys.glob("*.dedupped.bam")

ebg1 <- exonsBy( TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene" )
se1 <- summarizeOverlaps(features=ebg1, reads=bam.files,
        mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
saveRDS( se1, file = "rnaseq_read_count_entrez_id.rds" )
