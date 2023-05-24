library(enrichR)
library(openxlsx)
library(dplyr)

# import
list1 <- read.xlsx("/data/kharas/chuk/Z138_MSI2_Chiara/CRISPR_screen/CRISPR_list_1.xlsx")
list3.z138 <- read.xlsx("/data/kharas/chuk/Z138_MSI2_Chiara/CRISPR_screen/CRISPR_list_3.xlsx", sheet=1)
list3.oci <- read.xlsx("/data/kharas/chuk/Z138_MSI2_Chiara/CRISPR_screen/CRISPR_list_3.xlsx", sheet=2)

# clean up
list3.oci <- list3.oci %>% dplyr::select(-("X4")) #remove empty column carried over from excel
colnames(list3.oci) <- colnames(list3.z138)

list3.z138$FDR <- p.adjust( list3.z138$PValue, 'BH' )
list3.oci$FDR <- p.adjust( list3.oci$PValue, 'BH' )

# enriched and depleted
list1.enriched <- list1[ list1$logFC >=2 & list1$FDR < 0.05, ]
list1.depleted <- list1[ list1$logFC <= -2 & list1$FDR < 0.05, ]

list3.z138.enriched <- list3.z138[ list3.z138$LogFC >=1 & list3.z138$FDR < 0.05, ]
list3.z138.depleted <- list3.z138[ list3.z138$LogFC <= -1 & list3.z138$FDR < 0.05, ]

list3.oci.enriched <- list3.oci[ list3.oci$LogFC >=1 & list3.oci$FDR < 0.05, ]
list3.oci.depleted <- list3.oci[ list3.oci$LogFC <= -1 & list3.oci$FDR < 0.05, ]


#EnrichR
enrichr.database <- function(list.df, title) {
  
  dbs <- c("MSigDB_Oncogenic_Signatures", 
           "NCI-60_Cancer_Cell_Lines", 
           "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
           "WikiPathways_2019_Human",
           "KEGG_2019_Human",
           "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
           "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",
           "GO_Biological_Process_2018",
           "GO_Cellular_Component_2018",
           "GO_Molecular_Function_2018")
  
  enriched <- enrichr(c(list.df$NAME), dbs)
  
  for (i in 1:length(dbs)) {
    
    folder.name <- dbs[i]
    dir.create(folder.name)
    
    enriched.df <- enriched[[i]]
    
    write.csv(enriched.df, 
              paste0(folder.name, "/EnrichR_", title, "_", dbs[i], ".csv"),
              row.names = FALSE)
    
    png(paste0(folder.name, "/EnrichR_", title, "_", dbs[i], ".png"), 1200, 800, type="cairo")
    p <- plotEnrich(enriched[[i]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
    print(p)
    dev.off()
    
  }

}

enrichr.database(list1.depleted, "list1.depleted")
enrichr.database(list3.z138.enriched, "list3.z138.enriched")
enrichr.database(list3.z138.depleted, "list3.z138.depleted")
enrichr.database(list3.oci.enriched, "list3.oci.enriched") # Count isn't accurate for row 15
enrichr.database(list3.oci.depleted, "list3.oci.depleted")

# Filtered keyterms
# Oxidative phosphorylation, Mitochondrial function, sirtuin, RNA binding
enrichr.database <- function(list.df, title) {
  
  dbs <- c("MSigDB_Oncogenic_Signatures", 
           "NCI-60_Cancer_Cell_Lines", 
           "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
           "WikiPathways_2019_Human",
           "KEGG_2019_Human",
           "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
           "GTEx_Tissue_Sample_Gene_Expression_Profiles_down")
  
  enriched <- enrichr(c(list.df$NAME), dbs)
  
  for (i in 1:length(dbs)) {
    
    folder.name <- dbs[i]
    dir.create(folder.name)
    
    enriched.df <- enriched[[i]]
    enriched.df <- enriched.df[ enriched.df$Adjusted.P.value < 0.05, ]
    enriched.df <- enriched.df[ grepl("oxidative|phosphor|mito|RNA|sirtuin|bind", enriched.df$Term),]
    
    write.csv(enriched.df, 
              paste0(folder.name, "/EnrichR_Filtered_KeyTerms", title, "_", dbs[i], "_p-adj_lessthan_0.05.csv"),
              row.names = FALSE)
    
    png(paste0(folder.name, "/EnrichR_Filtered_KeyTerms", title, "_", dbs[i], "_p-adj_lessthan_0.05.png"), 1200, 800, type="cairo")
    p <- plotEnrich(enriched.df, showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
    print(p)
    dev.off()
    
  }
  
}

# enrichr.database(list1.depleted, "list1.depleted")
# enrichr.database(list3.z138.enriched, "list3.z138.enriched")
# enrichr.database(list3.z138.depleted, "list3.z138.depleted")
# enrichr.database(list3.oci.enriched, "list3.oci.enriched")
# enrichr.database(list3.oci.depleted, "list3.oci.depleted")













