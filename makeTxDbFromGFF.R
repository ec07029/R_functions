library("GenomicFeatures")

gff = makeTxDbFromGFF("/data/leslie/chuk/index_Karen/gencode.v19.annotation/gencode.v19.annotation.gtf")
saveRDS(gff, "/data/leslie/chuk/data_resources/gencode.v19.annotation.gtf_makeTxDbFromGFF.RDS")