# author: K Chu

# Get max edited site
library(dplyr)
df.filter.max.edit.freq <- df.filter %>% group_by(gene.symbol) %>% slice(which.max(diff.frequency))
df.filter.max.edit.freq.df <- as.data.frame(df.filter.max.edit.freq)
