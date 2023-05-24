# cannot work if dplyr is loaded
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

gene.symbol_to_entrez.id <- function(gene.list) {
  
  entrez.id <- select(hs, 
                      keys = gene.list,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
  
  return(entrez.id)
  
}

msi2.lethality.genes_entrez.id <- gene.symbol_to_entrez.id(msi2.lethality.filter$Gene)
