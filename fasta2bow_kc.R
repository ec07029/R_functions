# author: K Chu

inputFile <- "/Users/karen/Downloads/hy395-selex_embed-aee607680965/data/train_seqs.fasta"
labelFile <- "/Users/karen/Downloads/hy395-selex_embed-aee607680965/data/train_labels.txt"

#fasta2bow
# convert a fasta file to bag of words representation
library(Biostrings)
library(Matrix)
library(parallel)
load("~/Downloads/hy395-selex_embed-aee607680965/data/dict_8mer2wc8mer.rdt")
vocab_8mer <- colnames(oligonucleotideFrequency(DNAStringSet("AAAAAAAA"), 8))
pairwise.kmers <- pairwise.kmers[vocab_8mer,]

# functions
# Creates matrix (row = sequence in "input", columns = kmer) with frequency of kmer
makeFeats <- function(seqs, k, rc=F) {
  # Obtains the frequency of kmers of Oligo in sliding window (window length = k)
  a <- oligonucleotideFrequency(seqs, k)
  feats <- a
  if (rc) {
    seqs_rc <- reverseComplement(seqs)
    b <- oligonucleotideFrequency(seqs_rc, k)
    feats <- a+b
  }
  return(feats)
}

# Creating a string with labels and kmers that occur in input sequence
Feats2bow <- function(feats, labels, conn) {
  res <- mclapply(1:nrow(feats), function(j) {
    # Creates string of all kmers with value > 1
    tmp <- paste0(colnames(feats)[feats[j,]>0], collapse=" ")
    label <- labels[j]
    # Attach hashtag symbol to each element in label (elements separated by space)
    label_formated <- paste0(paste0("#",strsplit(label, " ")[[1]]), collapse=" ")
    # Paste label before kmer sequence
    tmp <- paste0(label_formated," ",tmp)
    return(tmp)
  }, mc.cores = 4)
  res <- unlist(res)
  write(res, conn, append = T)
}

get.time <- function () {
  return(as.numeric(format(Sys.time(), "%s")))
}

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly=TRUE)
inputFile <- args[1]
labelFile <- args[2]
outputFile <- args[3]
wc <- as.logical(args[4])

input <- readDNAStringSet(inputFile, "fasta")
labels <- as.character(read.table(labelFile, sep="\t", stringsAsFactors=F)$V1)
N <- length(input) # input sequences size
batch_size <- 5000
batches <- ceiling(N/batch_size)
print(paste0("total ",N, " sequences."))

time.start <- get.time()
file.create(outputFile)
fileCon <- file(outputFile, "a")
for (i in 1:batches) {
  # generate batch index
  start <- (i-1)*batch_size + 1
  end <- min(i*batch_size, N)
  # generate feature matrix and labels
  tmp <- Matrix(makeFeats(input[start:end], 8), sparse=T) # rows = each sequence in "input; col = kmer
  feats_tmp <- tmp%*%pairwise.kmers
  feats_nowc_tmp <- feats_tmp[, -grep("N",colnames(feats_tmp))]
  labels_tmp <- labels[start:end]
  if (wc == T) {
    Feats2bow(feats_tmp, labels_tmp, fileCon)
  } else {
    Feats2bow(feats_nowc_tmp, labels_tmp, fileCon)
  }
  time.end <- get.time()
  print(paste0("Process ",min(i*batch_size, N), " sequences takes ", round( (time.end-time.start)/60, 3), " min."))
}
close(fileCon)

