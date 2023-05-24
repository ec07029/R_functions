# Volcano plots: HyperTRIBE targets colored

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
    writeLines(hypertribe.query, paste0(output.dir, "hypertribe_query_filepath.txt"))
    
  } else { 
    hypertribe.query <- Sys.glob(paste0(hypertribe.dir, "*snp_counts_significance_fpkm.csv"))
    writeLines(hypertribe.query, paste0(output.dir, "hypertribe_query_filepath.txt")) }
  
  hypertribe.data <- read.csv(hypertribe.query)
  fpkm.index <- which( grepl("fpkm", colnames(hypertribe.data)) )
  hypertribe.data.filter <- hypertribe.data[ hypertribe.data$p.adj < 0.1 &
                                               hypertribe.data$diff.frequency >= 0.1 &
                                               hypertribe.data[,fpkm.index[1]] >= 5 &
                                               hypertribe.data[,fpkm.index[2]] >= 5, ]
  
  # Get DESeq2 data
  deseq.query <- Sys.glob(paste0(output.dir, "*DESeq2Results.csv"))
  res <- read.csv(deseq.query)
  
  # Prep data for plotting
  # Set significance column
  num.name <- group.two
  denom.name <- group.one
  padj.thres <- 0.05
  volcano.input <- as.data.frame( res[ order( res$padj, decreasing = F ), ] )
  volcano.input <- mutate(volcano.input, sig=ifelse(volcano.input$padj < padj.thres, "Sig", "Not Sig"))
  
  # Set hypertribe target column
  volcano.input$hypertribe.target <- "No"
  hypertribe.gene.targets <- unique(hypertribe.data.filter$entrez.id)
  
  if ( length(hypertribe.gene.targets) != 0 ) {
    
    volcano.input[ volcano.input$entrez.id %in% hypertribe.gene.targets, ]$hypertribe.target <- "Yes"
    
    # Full volcano plot
    png( paste0( output.dir, comparison.name, "_volcano_with_hypertribe_targets_colored.png"), type="cairo", 800, 800 )
    p <- ggplot(data=volcano.input, aes(x=log2FoldChange, y=-log10(pvalue), 
                                        colour= hypertribe.target, shape=sig)) + 
      geom_point(alpha=0.4, size=4) +
      theme_minimal() +
      #theme(legend.position="none") +
      geom_text_repel(data=volcano.input[1:30,], aes(label=gene), size = 5,box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"), color="black") +
      xlab(paste0("\nlog2(", num.name, "/", denom.name, ")")) + ylab("-log10(p)\n") +
      ggtitle(paste0(comparison.name, ";\npadj < ", padj.thres, ";\n# of DESeq2 sig genes: ", nrow( subset(volcano.input, sig=="Sig") ), "\n# of hypertribe gene targets: ", length(hypertribe.gene.targets), "\nHyperTRIBE filter: p.adj < 0.1, diff.freq >= 0.1, fpkm >= 5")) +
      scale_color_manual( values = c( "Yes"='red3', "No" ='black' ) ) + 
      theme(plot.title = element_text(size=20)) +
      theme(axis.text=element_text(size=20, color="black"), axis.title=element_text(size=20, color="black")) +
      theme(axis.line = element_line(colour = "black", size=1))
    print(p)
    dev.off()
    
    # Volcano plot with just hypertribe targets
    volcano.input.subset <- subset(volcano.input, hypertribe.target=="Yes")
    png( paste0( output.dir, comparison.name, "_volcano_HYPERTRIBE_TARGETS_PLOTTED_ONLY.png"), type="cairo", 800, 800 )
    p <- ggplot(data=volcano.input.subset, aes(x=log2FoldChange, y=-log10(pvalue), 
                                               shape= hypertribe.target, colour=sig)) + 
      geom_point(alpha=1, size=4) +
      theme_minimal() +
      #theme(legend.position="none") +
      geom_text_repel(data=volcano.input.subset[1:30,], aes(label=gene), size = 5,box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"), color="black") +
      xlab(paste0("\nlog2(", num.name, "/", denom.name, ")")) + ylab("-log10(p)\n") +
      ggtitle(paste0(comparison.name, ";\npadj < ", padj.thres, ";\n# of DESeq2 sig genes: ", nrow( subset(volcano.input, sig=="Sig") ), "\n# of hypertribe gene targets: ", length(hypertribe.gene.targets), "\nHyperTRIBE filter: p.adj < 0.1, diff.freq >= 0.1, fpkm >= 5")) +
      scale_color_manual( values = c( "Sig"='red3', "Not Sig" ='black' ) ) + 
      theme(plot.title = element_text(size=20)) +
      theme(axis.text=element_text(size=20, color="black"), axis.title=element_text(size=20, color="black")) +
      theme(axis.line = element_line(colour = "black", size=1))
    print(p)
    dev.off()
    
  } else { writeLines("No significant hypertribe targets", 
                      paste0(output.dir, "no_significant_hypertribe_targets.txt")) }
  
  
}











