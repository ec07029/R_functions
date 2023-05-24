library(Rsubread)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx) # xlsx package completely unusable

register(MulticoreParam(workers = 6))

bam.files <- Sys.glob("*.dedupped.bam")

fc <- featureCounts( files = bam.files,
               annot.inbuilt="mm10", 
               annot.ext = "/data/leslie/chuk/index_Karen/GRCm38_GencodeM15_with_humanMSI2-linker-ADAR/gencode.vM15.primary_assembly.annotation_MSI2-linker-ADAR_KC.gtf",
               isGTFAnnotationFile = TRUE, 
               useMetaFeatures = TRUE, 
               isPairedEnd = TRUE, 
               ignoreDup = TRUE )

saveRDS(fc, file="rnaseq_read_count_feature_counts.rds")


# ebg1 <- exonsBy( TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene" )
# se1 <- summarizeOverlaps(features=ebg1, reads=bam.files,
#         mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
# saveRDS( se1, file = "rnaseq_read_count_entrez_id.rds" )