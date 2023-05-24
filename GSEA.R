# author: "K Chu"

library(fgsea)
setwd("~/Desktop")
examplePathways <- read.table("lsc_ADAR_DMSO|ADAR_Ro_DESeq2Results_log2FC_lessthan-1.5_padj_lessthan_0.01_decreasing_log2FC_order_gene_list.grp")
exampleRanks <- read.table("41467_2019_10523_MOESM5_ESM_SHEET3_MOLM13_RNAseq_Ro_treatment_gene_list_decreasing_log2FC_order.rnk", header=TRUE)

mypathway <- list(mypathway=examplePathways[,1])

ranks <- exampleRanks
ranks <- setNames(ranks$log2FoldChange, ranks$gene)

fgseaRes <- fgsea(pathways = mypathway, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

plotEnrichment(mypathway[[1]], ranks)

gsea <- function(grp.file, rnk.file) {
  
  examplePathways <- read.table(grp.file)
  exampleRanks <- read.table(rnk.file, header=TRUE)
  
  mypathway <- list(mypathway=examplePathways[,1])
  
  ranks <- exampleRanks
  ranks <- setNames(ranks$log2FoldChange, ranks$gene)

  fgseaRes <- fgsea(pathways = mypathway, 
                    stats = ranks,
                    minSize=15,
                    maxSize=500,
                    nperm=10000)
  
  plotEnrichment(mypathway[[1]], ranks)
  
}

gsea("lsc_MIGR1_DMSO|MIGR1_Ro_DESeq2Results_log2FC_lessthan-2_padj_lessthanorequalto_0.05_gene_list_decreasing_log2FC_order.grp",
     "41467_2019_10523_MOESM5_ESM_SHEET4_shRNAs_in_AML_and_CML-BC_cell_lines_from_Kharas_etal_decreasing_log2FC_order.rnk")


