# author: "K Chu"

library(dplyr)

setwd("/Users/chuk/Documents/Syncrip/distance to nearest motif/data/")

# Read in files
hsc <- read.csv("181214.HSC-results-fpkm.csv")
mpp <- read.csv("181214.MPP-results-fpkm.csv")

# Clean up: remove rownames X column (unnecessary column)
hsc <- hsc %>% select(-X)
mpp <- mpp %>% select(-X)

# Calculate editing frequency
calculate.edit.frequency <- function(df1) {
  
  ref.counts <- df1[,grep('ref.count', colnames(df1))]
  alt.counts <- df1[,grep('alt.count', colnames(df1))]
  
  adar.index <- which( grepl("ADAR", colnames(ref.counts)) )
  mig.index <- which( grepl("MIG", colnames(ref.counts)) ) 
  
  alt.freq <- alt.counts / ( ref.counts + alt.counts )
  df.stats <- data.frame(diff.frequency = rowMeans(alt.freq[, adar.index ], na.rm=T) - rowMeans(alt.freq[, mig.index ],na.rm=T),
                         ADAR.frequency = rowMeans(alt.freq[,adar.index], na.rm=T),
                         MIG.frequency = rowMeans(alt.freq[, mig.index ], na.rm=T) )
  
}

hsc.edit.freq <- calculate.edit.frequency(hsc)
mpp.edit.freq <- calculate.edit.frequency(mpp)

hsc <- cbind( hsc[,1:8], hsc.edit.freq, hsc[,-(1:8)] )
mpp <- cbind( mpp[,1:8], mpp.edit.freq, mpp[,-(1:8)] )

write.csv(hsc, "181214.HSC-results-fpkm_ADDED_DIFF_FREQ-KAREN.csv", row.names = FALSE)
write.csv(mpp, "181214.MPP-results-fpkm_ADDED_DIFF_FREQ-KAREN.csv", row.names = FALSE)
