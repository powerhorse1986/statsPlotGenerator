if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("tximport")
BiocManager::install("RCurl")
BiocManager::install("DESeq2")
BiocManager::install("readr")
BiocManager::install("rhdf5")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(tximport)
library(DESeq2)
library(rhdf5)
library(EnsDb.Hsapiens.v75)

samples <- read.table("/home/lima/Project/simulation/Comparision/samples.txt", header = TRUE)
samples

###########################################################
# kallisto
###########################################################
kallPath <- "/home/lima/Project/simulation/Comparision/Kall"
txdb <- EnsDb.Hsapiens.v75
df <- transcripts(txdb, return.type = "DataFrame")
tx2gene <- df[, 8 : 7]  # tx ID, then gene ID

###########################################################
# kallisto TSV files for real data
###########################################################
files <- file.path(kallPath, "Real", samples$GeneID, "abundance.tsv")
names(files) <- samples$GeneID
# some problems with rhdf5 packages
#txi.kall <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
head(txi.kall.tsv$counts)
dds <- DESeqDataSetFromTximport(txi.kall.tsv, samples, ~1)

# performe a rlog/varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds, blind = FALSE)

# convert counts to a dataframe and save it into csv
rld_df <- data.frame(assay(rld))
write.csv(rld_df, file.path(kallPath, "Real", "KallistoReadCounts.csv"))

##############################################################
# RSEM counts
##############################################################
RSEMPath <- "/home/lima/Project/simulation/Comparision/RSEM"
files <- file.path(RSEMPath, "Real", samples$GeneID, "genes.results")
names(files) <- samples$GeneID

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem)
txi.rsem$length <- txi.rsem$length + 0.5
# create DESeqDataSet object from tximport, ~1 indicating we are simply exploring the data
dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~1)

# rld <- varianceStabilizingTransformation(dds, blind = FALSE) # varianceStabilizingTransformation is not required due to it normalizes the coutns
# rld_df <- data.frame(assay(rld))
mcols(dds) <- DataFrame(mcols(dds))
#rld <- rld[rowSums(counts(rld)) > 2.5]
cutOff = 0

while(cutOff <= 10) {
  # filtered <- rld[rowSums(assay(rld)) >= cutOff, ]
  # keep <- rowSums(counts(rld)) >= cutOff
  # output <- rld_df[keep, ]
  output <- dds[rowSums(counts(dds)) >= cutOff]
  output <- data.frame(assay(output))
  outPutName = paste("CountsCutOffValue", cutOff, ".csv", sep = "")
  print(outPutName)
  outPutPath = file.path(RSEMPath, "Real", outPutName)
  write.csv(output, outPutPath)
  cutOff = cutOff + 2.5
}