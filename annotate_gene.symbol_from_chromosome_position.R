library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

exbygene <- exonsBy( TxDb.Hsapiens.UCSC.hg19.knownGene, "gene" )
ol <- findOverlaps( iclip.gr, exbygene )
entrez.id <- names( exbygene )[ subjectHits(ol) ]
gene.symbol <- mget( entrez.id, org.Hs.egSYMBOL, ifnotfound = NA )
gene.symbol <- unlist( gene.symbol )
anno <- DataFrame( entrez.id, gene.symbol )

iclip$entrez.id <- NA
iclip$gene.symbol <- NA
iclip$entrez.id[ queryHits(ol) ] <- anno$entrez.id
iclip$gene.symbol[ queryHits(ol) ] <- anno$gene.symbol