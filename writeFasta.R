# author: K Chu

# write fasta file using dataframe that has column names "ensembl_gene_id" and "3utr".
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"ID"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

setwd("/Users/chuk/Documents/MSI2-SYNCRIP/output_files/")

# Write all sequences
writeFasta(as.data.frame(random.genomic.sequences.gr.in.HSC.no.overlap.with.background),
           paste0("HOMER_BACKGROUND_mouse_genome_", nucleotide.distance.to.expand, "bp_leftright_window_in_HSC.fa"))
write.csv(random.genomic.sequences.gr.in.HSC.no.overlap.with.background$sequence, "background_sequences_in_HSC.csv")
