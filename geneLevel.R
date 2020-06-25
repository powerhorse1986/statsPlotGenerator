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
BiocManager::install("RMariaDB")
BiocManager::install("tidyverse")

library(tximport)
library(DESeq2)
library(rhdf5)
library(EnsDb.Hsapiens.v75)
library(ensembldb)
library(GenomicFeatures)
library(RMariaDB)
library(tidyverse)

samples <- read.table("/home/lima/Project/simulation/Comparision/pt.txt", header = TRUE)
samples

###########################################################
# kallisto
###########################################################
kallPath <- "/home/lima/Project/simulation/Comparision/Kall/Tmp"
gtffile <- file.path("/home/lima/Project/simulation/GRCh37", "Homo_sapiens.GRCh37.87.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
GR <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx2gene <- mcols(GR)
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene[,  ] <- lapply(tx2gene[, ], as.character)

###########################################################
# kallisto TSV files for real data
###########################################################
files <- file.path(kallPath, samples$GeneID, "abundance.tsv")
names(files) <- samples$GeneID
# some problems with rhdf5 packages
#txi.kall <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txi.kall.tsv$length <- txi.kall.tsv$length + 0.5
# txi.kall.tsv$counts <- log2(txi.kall.tsv$counts + 1)
head(txi.kall.tsv$counts)
dds <- DESeqDataSetFromTximport(txi.kall.tsv, samples, ~1)
nrow(dds)
# performe a rlog/varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds, blind = FALSE)
# type(rld)
# convert counts to a dataframe and save it into csv
# rld_df <- data.frame(assay(rld))
df <- data.frame(assay(rld))
write.csv(df, file.path("/home/lima/Project/simulation/Comparision/pt/Normalized", "NormalizedKallCounts.csv"))

##############################################################
# RSEM counts
##############################################################
RSEMPath <- "/home/lima/Project/simulation/Comparision/RSEM/Tmp"
# files <- file.path(RSEMPath, pt$GeneID, "genes.results")
files <- file.path(RSEMPath, "PT", "genes.results")
#names(files) <- pt$GeneID
names(files) <- "PT"
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

files1 <- file.path(RSEMPath, "PTSimulated", "genes.results")
#names(files) <- pt$GeneID
names(files1) <- "PTSimulated"
txi.rsem1 <- tximport(files1, type = "rsem", txIn = FALSE, txOut = FALSE)

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
# txdb <- EnsDb.Hsapiens.v75
# df <- transcripts(txdb, return.type = "DataFrame")
# tx2gene <- df[, 8 : 7]  # tx ID, then gene ID

#############################################################
# txdb1 <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 75, circ_seqs = character(), 
                             # server="ensembldb.ensembl.org", username="anonymous", password=NULL, port=0L, tx_attrib=NULL)
# GR1 <- transcripts(txdb1, columns = c("tx_name", "gene_id"))
# df1 <- mcols(GR1)
# tx2gene1 <- df1[, 1 : 2]
# colnames(tx2gene1) <- c("TXNAME", "GENEID")
# tx2gene1[, ] <- lapply(tx2gene1[, ], as.character)
# tx2gene1

gtffile <- file.path("/home/lima/Project/simulation/GRCh37", "Homo_sapiens.GRCh37.87.gtf")
txdb2 <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
GR2 <- transcripts(txdb2, columns = c("tx_name", "gene_id"))
tx2gene2 <- mcols(GR2)
colnames(tx2gene2) <- c("TXNAME", "GENEID")
tx2gene2[,  ] <- lapply(tx2gene2[, ], as.character)
#############################################################


# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# k <- keys(txdb, keytype = "GENEID")
# df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
# tx2gene <- df[, 2:1]
# class(tx2gene)

###########################################################
# kallisto TSV files for real data
###########################################################
files <- file.path(kallPath, "Real", samples$GeneID, "abundance.tsv")
names(files) <- samples$GeneID
# some problems with rhdf5 packages
#txi.kall <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene2, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
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


#############################################################
# Normalize samples by treatment group vs tumor
#############################################################
samples <- read.table("/home/lima/Project/simulation/Comparision/veh.txt", header = TRUE)
samples
kallPathThree <- "/home/lima/Project/simulation/Comparision/Kall"
gtffile <- file.path("/home/lima/Project/simulation/GRCh37", "Homo_sapiens.GRCh37.87.gtf")
txdb3 <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
GR3 <- transcripts(txdb3, columns = c("tx_name", "gene_id"))
tx2gene3 <- mcols(GR3)
colnames(tx2gene3) <- c("TXNAME", "GENEID")
tx2gene3[,  ] <- lapply(tx2gene3[, ], as.character)

files <- file.path(kallPathThree, "Real", samples$GeneID, "abundance.tsv")
names(files) <- samples$GeneID

txi.kall.tsv <- tximport(files, type = "kallisto", dropInfReps=TRUE, tx2gene = tx2gene2, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
txi.kall.tsv$length <- txi.kall.tsv$length + 0.5

# condition = factor(samples$Condition, levels = c("Veh", "High", "Low", "Pt"))
# condition = factor(samples$Condition, levels = c("Low", "Pt"))
# samples$Condition <- condition
dds <- DESeqDataSetFromTximport(txi.kall.tsv, samples, ~ 1)
nrow(dds)
rld <- varianceStabilizingTransformation(dds, blind = FALSE)
nrow(rld)
rld_df <- data.frame(assay(rld))
write.csv(rld_df, file.path("/home/lima/Project/simulation/Comparision/veh/Normalized", "kallistoVeh.csv"))