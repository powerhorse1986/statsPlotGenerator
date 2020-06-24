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

library(tximport)
library(DESeq2)
library(rhdf5)
library(EnsDb.Hsapiens.v75)
library(ensembldb)

pt <- read.table("/home/lima/Project/simulation/Comparision/pt.txt", header = TRUE)
pt

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
files <- file.path(kallPath, "Real", pt$GeneID, "abundance.tsv")
names(files) <- pt$GeneID
# some problems with rhdf5 packages
#txi.kall <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txi.kall.tsv$length <- txi.kall.tsv$length + 0.5
# txi.kall.tsv$counts <- log2(txi.kall.tsv$counts + 1)
head(txi.kall.tsv$counts)
dds <- DESeqDataSetFromTximport(txi.kall.tsv, pt, ~1)

# performe a rlog/varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds, blind = FALSE)
# type(rld)
# convert counts to a dataframe and save it into csv
# rld_df <- data.frame(assay(rld))
df <- data.frame(assay(dds))
write.csv(df, file.path("/home/lima/Project/simulation/Comparision/pt", "KallistoReadCounts.csv"))

##############################################################
# RSEM counts
##############################################################
RSEMPath <- "/home/lima/Project/simulation/Comparision/RSEM"
files <- file.path(RSEMPath, "Real", pt$GeneID, "genes.results")
names(files) <- pt$GeneID

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

head(txi.rsem)
txi.rsem$length <- txi.rsem$length + 0.5
# txi.rsem$counts <- log2(txi.rsem$counts + 1)
# create DESeqDataSet object from tximport, ~1 indicating we are simply exploring the data
# countData <- round(txi.rsem[['counts']])
# colData(txi.rsem)
dds <- DESeqDataSetFromTximport(txi.rsem, pt, ~1)

nrow(dds)
# dds[rowSums(counts(dds)) > 2.5, 0]
# rld <- varianceStabilizingTransformation(dds, blind = FALSE) # varianceStabilizingTransformation is not required due to it normalizes the coutns
# rld_df <- data.frame(assay(rld))
# mcols(dds) <- DataFrame(mcols(dds))
#rld <- rld[rowSums(counts(rld)) > 2.5]
cutOff = 0

while(cutOff <= 10) {
  keep <- rowSums(counts(dds)) >= cutOff
  output <- dds[keep, ]
  output <- data.frame(assay(output))
  outPutName = paste("CountsCutOffValue", cutOff, ".csv", sep = "")
  print(outPutName)
  outPutPath = file.path("/home/lima/Project/simulation/Comparision/pt", outPutName)
  write.csv(output, outPutPath)
  cutOff = cutOff + 2.5
}


########################################################################################################################################################
# not that raw
samples <- read.table("/home/lima/Project/simulation/Comparision/samples.txt", header = TRUE)
samples

###########################################################
# kallisto
###########################################################
kallPath <- "/home/lima/Project/simulation/Comparision/Kall"
txdb <- EnsDb.Hsapiens.v75
df <- transcripts(txdb, return.type = "DataFrame")
tx2gene <- df[, 8 : 7]  # tx ID, then gene ID
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# k <- keys(txdb, keytype = "GENEID")
# df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
# tx2gene <- df[, 2:1]
# tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# seqlevelsStyle(txdb) <- "Ensembl"
# k <- keys(txdb, keytype = "GENEID")
# df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
# tx2gene <- df[, 2:1]
# tx2gene <- select(txdb, k, "GENEID")
###########################################################
# kallisto TSV files for real data
###########################################################
files <- file.path(kallPath, "Real", samples$GeneID, "abundance.tsv")
names(files) <- samples$GeneID
# some problems with rhdf5 packages
#txi.kall <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txi.kall.tsv$length <- txi.kall.tsv$length + 0.5
# txi.kall.tsv$counts <- log2(txi.kall.tsv$counts + 1)
head(txi.kall.tsv$counts)
condition = factor(samples$Condition, levels = c("Veh", "High", "Low", "Pt"))
samples$Condition <- condition
dds <- DESeqDataSetFromTximport(txi.kall.tsv, samples, ~ Condition)
nrow(dds)
# performe a rlog/varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds, blind = FALSE)
nrow(rld)
# type(rld)
# convert counts to a dataframe and save it into csv
rld_df <- data.frame(assay(rld))

write.csv(rld_df, file.path("/home/lima/Project/simulation/Comparision/all/normalized", "NormalizedKallistoReadCounts.csv"))

##############################################################
# RSEM counts
##############################################################
RSEMPath <- "/home/lima/Project/simulation/Comparision/RSEM"
files <- file.path(RSEMPath, "Real", samples$GeneID, "genes.results")
names(files) <- samples$GeneID

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

head(txi.rsem)
txi.rsem$length <- txi.rsem$length + 0.5
condition = factor(samples$Condition, levels = c("Veh", "High", "Low", "Pt"))
samples$Condition <- condition
# txi.rsem$counts <- log2(txi.rsem$counts + 1)
# create DESeqDataSet object from tximport, ~1 indicating we are simply exploring the data
# countData <- round(txi.rsem[['counts']])
# colData(txi.rsem)
dds <- DESeqDataSetFromTximport(txi.rsem, samples, design = ~ Condition)

nrow(dds)
# dds[rowSums(counts(dds)) > 2.5, 0]
rld <- varianceStabilizingTransformation(dds, blind = FALSE) # varianceStabilizingTransformation is not required due to it normalizes the coutns
nrow(rld)
rld_df <- data.frame(assay(rld))
# mcols(dds) <- DataFrame(mcols(dds))
# rld <- rld[rowSums(counts(rld)) > 2.5]
cutOff = 0

while(cutOff <= 10) {
  keep <- rowSums(counts(dds)) >= cutOff
  output <- rld_df[keep, ]
  # output <- data.frame(assay(output))
  # output = data.frame(assay(output))
  outPutName = paste("NormalizedCountsCutOffValue", cutOff, ".csv", sep = "")
  print(outPutName)
  outPutPath = file.path("/home/lima/Project/simulation/Comparision/all/normalized", outPutName)
  write.csv(output, outPutPath)
  cutOff = cutOff + 2.5
}
