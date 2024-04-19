#Script DESEq2 2024
#Install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages(apeglm)

install.packages("pheatmap")

install.packages("xlsx")
library("xlsx")
library(DESeq2)
library(apeglm)
library(pheatmap)

#Set directory
setwd("C:/Users/Ana Luiza/Desktop/Deseq2_2024")
getwd()

#load count matrix
dat <- read.table("clipboard", header = T, row.names = 1)
dat
info <- read.table("colData_leaf.txt", header = T, sep = '\t')

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

#if want remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#main DESEq
ddsDE <- DESeq(dds)

#export normalized read counts
normCounts <- counts(ddsDE, normalized = T)
normCounts
write.csv(normCounts, "normalized_leaf.csv")
write.table(normCounts, sep = "\t", "Normalized_Leaf.txt")


#DESEq results
res <- results(ddsDE, alpha = 0.05)
res

resOrdered <- res[order(res$padj),]
resOrdered
write.csv(resOrdered, "deSeq_order_leaves.csv")
write.xlsx(resOrdered, "deSeq_order_leaves.xlsx")

head(res[order(res$padj),], 20)


#summary results
summary(res)

#get significant genes
resSig <- subset(res, padj < 0.05)
resSig
res_Sig <- as.data.frame(resSig)
write.table(res_Sig, sep = "\t", "resSigLeaf.txt")

#get up regulated genes
upReg <- subset(resSig, log2FoldChange > 0)
upReg
up_reg <- as.data.frame(upReg)
write.table(up_reg, sep = "\t", "Up_reg_Leaf.txt")

#get down regulated genes
downReg <- subset(resSig, log2FoldChange < 0)
downReg
down_Reg <- as.data.frame(downReg)
write.table(downReg, sep = "\t", "Down_reg_Leaf.txt")

plotMA(ddsDE, ylim=c(-10,10))

plotCounts(dds, gene="TE_00001097_Chr8_0__SIRE", intgroup="condition")
hist( res$pvalue, breaks=20, col="grey" )

#plots
normCounts <- read.csv("normalized_leaf.csv", row.names = 1)
deSeqRes <- read.csv("deSeq_order_leaves.csv", row.names = 1)
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")
deSeqRes

#pheatmap
signi <- subset(deSeqRes,padj <= 0.05)
allsig <- merge(normCounts, signi, by=0)

sigCounts <- allsig[,2:7]
row.names(sigCounts) <- allsig$Row.names

pheatmap(sigCounts)
