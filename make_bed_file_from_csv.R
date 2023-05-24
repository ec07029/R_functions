# create bed file from csv

setwd("~/mount/chuk/sc-HyperTRIBE/from_Nathan/10X_genomics_Minimap_aligner/Analysis/")

df <- read.csv("MSI2-ADAR_10X_snp_counts_significance_FILTERED_p.adj_0.1_diff.freq_0.1.csv")

df.subset <- df[,1:2]

df.subset$end <- df.subset$start
df.subset$start <- df.subset$start - 1

write.table(df.subset, "MSI2-ADAR_10X_snp_counts_significance_FILTERED_p.adj_0.1_diff.freq_0.1.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE)
