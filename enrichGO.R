library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(lattice)
library(dplyr)
library(DESeq2)
library(SummarizedExperiment)
library(reshape2)
library(openxlsx)
library(ggplot2)


background <- unique( as.character( rownames(reads.mat) ) ) # All data

run.enrichrgo <- function(entrez.id, background, cell.name) {
  
  targets <- unique( as.character( unique( entrez.id ) ) )
  background <- unique( as.character( background ) )
  
  ncg <- enrichGO(targets,
                  #organism = 'mmu',
                  #keyType = "kegg",
                  OrgDb = "org.Mm.eg.db",
                  pAdjustMethod = "BH",
                  universe = background )
  
  ncg.results <- ncg$Description
  
  ncg.plot <- enrichplot::dotplot(ncg, x="Count", showCategory=length(ncg$Description))
  png(paste0(folder, "enrichr_output/", cell.name, "_enrichGO_fpkm_greaterthanorequalto_1.png"), 1200, 800)
  print(ncg.plot)
  dev.off()
  
  return(ncg.results)
  
}

lt.results <- run.enrichrgo(lt.yuheng.filter$entrez.id, rownames(lt_fpkm.filter),
                            "LT_yuheng_sig_genes_with_with_cell.type_gene.expression_background")