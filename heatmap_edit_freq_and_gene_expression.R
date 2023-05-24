# Heatmap: Edit frequency and Gene Expression of HyperTRIBE gene targets

# setwd("/data/leslie/chuk/Z138_MSI2_Chiara/Analysis_trimmed_reads_test_auto-script/DESeq2/")

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

# Get list of all possible comparisons for DESeq2
snp.counts <- read.xlsx("/data/leslie/chuk/Z138_MSI2_Chiara/Project_10874_B/JAX_0444/trim_galore/Z138_MSI2_Chiara_snp_counts.xlsx")
sample.names <- gsub("[0-9]_IGO.*", "", colnames(snp.counts))
sample.names.only <- sample.names[ grepl("Sample", sample.names) ]
sample.names.only.dedup <- unique(sample.names.only)
comparison.matrix <- combn(sample.names.only.dedup, 2) # gives all possible pairwise combos

# Import VST reads
reads.vst <- read.csv("reads.vst.csv")
rownames(reads.vst) <- reads.vst$X
reads.vst <- reads.vst %>% select(-c("X"))

# Z-transform
reads.vst.transform <- t( scale(t(reads.vst)) )
reads.vst.transform.df <- as.data.frame(reads.vst.transform)
write.csv(reads.vst.transform.df, "reads.vst_z-transformed.csv")

for (i in 1:ncol(comparison.matrix)) {
# for (i in 1:2) {
  
  print(i)
  
  query.comparison <- comparison.matrix[,i]
  group.one <- query.comparison[1]
  group.two <- query.comparison[2]
  
  # Set comparison name
  comparison.name <- paste0(group.two, "_vs_", group.one)
  output.dir <- paste0(comparison.name, "/")
  hypertribe.dir <- paste0("../beta-binomial_on_atlas/", 
                           group.two, "_vs.over_", group.one,
                           "/output/")
  
  # Get hypertribe data
  atlas.name <- c("Sample_MSI2ADAR-Control-_vs_Sample_MIGR1-Control-", 
                  "Sample_MSI2ADAR-GSK591-_vs_Sample_MIGR1-GSK591-",
                  "Sample_MSI2ADAR-GSK712-_vs_Sample_MIGR1-GSK712-",
                  "Sample_MSI2ADAR-Ro-_vs_Sample_MIGR1-Ro-")
  if ( (comparison.name %in% atlas.name) == TRUE ) {
    
    hypertribe.query <- Sys.glob(paste0("../", group.two, "_vs.over_", group.one, "/output/*snp_counts_significance_fpkm.csv"))
    writeLines(hypertribe.query, paste0(output.dir, "hypertribe_query_filepath_for_heatmap_edit_freq_and_gene_expression.txt"))
    
  } else { 
    hypertribe.query <- Sys.glob(paste0(hypertribe.dir, "*snp_counts_significance_fpkm.csv"))
    writeLines(hypertribe.query, paste0(output.dir, "hypertribe_query_filepath_for_heatmap_edit_freq_and_gene_expression.txt")) 
    }
  
  # Filter HyperTRIBE data
  hypertribe.data <- read.csv(hypertribe.query)
  fpkm.index <- which( grepl("fpkm", colnames(hypertribe.data)) )
  hypertribe.data.filter <- hypertribe.data[ hypertribe.data$p.adj < 0.1 &
                                             hypertribe.data$diff.frequency >= 0.1 &
                                             hypertribe.data[,fpkm.index[1]] >= 5 &
                                             hypertribe.data[,fpkm.index[2]] >= 5, ]
  
  if ( nrow(hypertribe.data.filter) == 0 ) {
    
    writeLines("No significant hypertribe targets", 
               paste0(output.dir, "no_significant_hypertribe_targets_heatmap_edit.freq_and_gene.expression.txt"))
    print("next")
    next
    
  }
  
  # Get maximum edit frequency for each significant HyperTRIBE gene target
  hypertribe.filter.max.edit.freq <- hypertribe.data.filter %>% group_by(entrez.id) %>% slice(which.max(diff.frequency))
  hypertribe.filter.max.edit.freq.df <- as.data.frame(hypertribe.filter.max.edit.freq)
  write.csv(hypertribe.filter.max.edit.freq.df,
            paste0(output.dir, comparison.name, "_filtered_max.diff.frequency.only.csv"),
            row.names = FALSE)
  
  # Subset VST z-transformed reads for significant hypertribe targets
  reads.vst.transform.df.subset <- reads.vst.transform.df[ rownames(reads.vst.transform.df) %in% hypertribe.filter.max.edit.freq.df$entrez.id, ]
  
  # Prepare heatmap input
  reads.vst.transform.df.subset$entrez.id <- rownames(reads.vst.transform.df.subset)
  freq.index <- which( grepl(".frequency", colnames(hypertribe.filter.max.edit.freq.df)) )
  hypertribe.filter.max.edit.freq.values <- hypertribe.filter.max.edit.freq.df %>% select(c("entrez.id", "gene.symbol"))
  hypertribe.filter.max.edit.freq.values <- cbind(hypertribe.filter.max.edit.freq.values,
                                                  hypertribe.filter.max.edit.freq.df[,freq.index])
  
  heatmap.combined <- merge(reads.vst.transform.df.subset, 
                            hypertribe.filter.max.edit.freq.values, 
                            by="entrez.id", all=TRUE)
  rownames(heatmap.combined) <- heatmap.combined$gene.symbol
  heatmap.combined.final <- heatmap.combined %>% select(-c("entrez.id", "gene.symbol"))
  write.csv(heatmap.combined.final, paste0(output.dir, "heatmap_combined_final.csv"), row.names = FALSE)
  
  
  # Get edit freq for heatmap
  heatmap.diff.freq.input <- heatmap.combined.final [ ,grepl(".frequency", colnames(heatmap.combined.final))]
  
  # Get gene expression for heatmap
  freq.index <- which( grepl(".frequency", colnames(heatmap.combined.final)))
  heatmap.gene.expression.input <- heatmap.combined.final[,-freq.index]
  colnames(heatmap.gene.expression.input) <- gsub("_IGO.*", "", colnames(heatmap.gene.expression.input))
  
  # Plot heatmap
  # rows clustered
  png(paste0( output.dir, comparison.name, "_heatmap_edit.freq_and_gene.expression_rows.only.clustered.png"), type="cairo", 1000, 800)
  p <- Heatmap(as.matrix(heatmap.diff.freq.input), use_raster=T, raster_device = "CairoPNG",
               column_title = paste0("Edit Freq"),
               name="Edit Freq",
               show_row_names = TRUE,
               cluster_columns = FALSE) +
    Heatmap(as.matrix(heatmap.gene.expression.input), use_raster=T, raster_device = "CairoPNG",
            column_title = paste0(comparison.name, "\nGene Expression of Significant HyperTRIBE gene targets\np.adj < 0.1, diff.freq >= 0.1, fpkm >= 5\n# of HyperTRIBE targets: ", length(unique(hypertribe.data.filter$entrez.id))),
            name="Z-transformed VST read counts",
            show_row_names = TRUE,
            cluster_columns = FALSE) 
  print(p)
  
  # Add borders
  decorate_heatmap_body("Edit Freq", {
    
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
    
  })
  
  decorate_heatmap_body("Z-transformed VST read counts", {
    
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
    
  })
  
  dev.off()
  
  # Write gene symbol order in heatmap
  write.csv( rownames(heatmap.combined.final)[row_order(p)],
             paste0( output.dir, comparison.name, 
                     "_heatmap_edit.freq_and_gene.expression_gene_symbol_order_rows.only.clustered.csv"), 
             row.names = FALSE )
  
  # rows and columns clustered
  png(paste0( output.dir, comparison.name, "_heatmap_edit.freq_and_gene.expression.png"), type="cairo", 1000, 800)
  p <- Heatmap(as.matrix(heatmap.diff.freq.input), use_raster=T, raster_device = "CairoPNG",
               column_title = paste0("Edit Freq"),
               name="Edit Freq",
               show_row_names = TRUE,
               cluster_columns = TRUE) +
    Heatmap(as.matrix(heatmap.gene.expression.input), use_raster=T, raster_device = "CairoPNG",
            column_title = paste0(comparison.name, "\nGene Expression of Significant HyperTRIBE gene targets\np.adj < 0.1, diff.freq >= 0.1, fpkm >= 5\n# of HyperTRIBE targets: ", length(unique(hypertribe.data.filter$entrez.id))),
            name="Z-transformed VST read counts",
            show_row_names = TRUE,
            cluster_columns = TRUE) 
  print(p)
  
  # Add borders
  decorate_heatmap_body("Edit Freq", {
    
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
    
  })
  
  decorate_heatmap_body("Z-transformed VST read counts", {
    
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
    
  })
  
  dev.off()
  
  # Write gene symbol order in heatmap
  write.csv( rownames(heatmap.combined.final)[row_order(p)],
             paste0( output.dir, comparison.name, 
                     "_heatmap_edit.freq_and_gene.expression_gene_symbol_order.csv"), 
             row.names = FALSE )
  
}









