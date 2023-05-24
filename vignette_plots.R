library(DESeq2)
library(data.table)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
library(ggplot2)
library("ggrepel")
library(data.table)
library(dplyr)
library(openxlsx)
library(ComplexHeatmap)

# setwd("~/mount/chuk/METTL3_Hanzhi/Project_10991_B_WT_KO_LSK_SON_OE/JAX_0453/align_to_normal_reference_genome/DESeq2/DESeq2_vignette_plots/")

# Import data
reads <- readRDS("../../rnaseq_read_count_entrez_id.rds")
reads.df <- as.data.frame( assay(reads) )

# Remove all rows with only zero value
reads.filtered <- reads.df[apply(reads.df, 1, function(x) !all(x==0)),]

# Create dds
sample.id <- colnames(reads.filtered)
sample.id <- strsplit(sample.id, "_IGO")
sample.id <- sapply(sample.id, "[[", 1)

condition <- ifelse(grepl("K", colnames(reads.filtered)), "METTL3.KO", "METTL3.WT")
condition2 <- ifelse(grepl("WS|KS", colnames(reads.filtered)),
                     "SON.OE", "SON.WT")

coldata <- data.frame(Sample = as.factor(sample.id), 
                      condition = as.factor(condition),
                      condition2 = as.factor(condition2))

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(reads.filtered), 
  colData = coldata, 
  design = ~condition)

dds$condition <- relevel(dds$condition, "METTL3.WT") # intercept that DESeq2 calculates depends on the first factor. 

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

# Effects of transformations on the variance
ntd <- normTransform(dds)
library("vsn")
png("meanSdPlot.png", type="cairo", 800, 800)
meanSdPlot(assay(ntd))
dev.off()

# Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition","Sample")])
png("heatmap_of_count_matrix.png", type="cairo", 800, 800)
p <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
print(p)
dev.off()

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- gsub("_IGO.*", "", rownames(sampleDistMatrix))
colnames(sampleDistMatrix) <- gsub("_IGO.*", "", colnames(sampleDistMatrix))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("heatmap_of_sample-to-sample_distances.png", type="cairo", 800, 800)
p <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
print(p)
dev.off()

# PCA
pcaData <- plotPCA(vsd, intgroup=c("Sample", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=Sample, shape=cell)) +
png("pca_mettl3_condition.png", type="cairo", 800, 800)
p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=5) +
  #geom_label(
    #label=pcaData$Sample, 
    #size = 3,
    #nudge_x = 2, nudge_y = 0 ) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(plot.title = element_text(size=40)) +
  theme(axis.text=element_text(size=30, color="black"), axis.title=element_text(size=30, color="black")) +
  theme(axis.line = element_line(colour = "black", size=1))
print(p)
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("Sample", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=Sample, shape=cell)) +
png("pca_son_condition.png", type="cairo", 800, 800)
p <- ggplot(pcaData, aes(PC1, PC2, color=condition2)) +
  geom_point(size=5) +
  #geom_label(
  #label=pcaData$Sample, 
  #size = 3,
  #nudge_x = 2, nudge_y = 0 ) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(plot.title = element_text(size=40)) +
  theme(axis.text=element_text(size=30, color="black"), axis.title=element_text(size=30, color="black")) +
  theme(axis.line = element_line(colour = "black", size=1))
print(p)
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("Sample", "condition"), returnData=TRUE)
pcaData$Sample <- gsub("Sample_", "", pcaData$Sample)
pcaData$Sample <- gsub("_IGO.*", "", pcaData$Sample)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=Sample, shape=cell)) +
png("pca_son_condition_labels.png", type="cairo", 800, 800)
p <- ggplot(pcaData, aes(PC1, PC2, color=condition2)) +
  geom_point(size=5) +
  geom_label(
  label=pcaData$Sample,
  size = 3,
  nudge_x = 2, nudge_y = 0 ) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(plot.title = element_text(size=40)) +
  theme(axis.text=element_text(size=30, color="black"), axis.title=element_text(size=30, color="black")) +
  theme(axis.line = element_line(colour = "black", size=1))
print(p)
dev.off()

