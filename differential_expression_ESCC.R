#Author: K.T.Shreya Parthasarathi
#Date: 23/03/2022
#Purpose: Estimate differentially expressed ion channels, lipid metabolism genes and EMT-related genes from rna-seq datasets of patients with ESCC (DESeq2)

#Install libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

setwd("path\\to\\current\\working\\directory")

#Import libraries
library("DESeq2")
library("ggplot2")

#Data import and preprocessing
file = read.table("ESCC_EMT_GSE32424.tsv", sep = '\t',header = TRUE)
dim(file)
file = file[,2:ncol(file)]

rownames(file) = file$hgnc_symbol
dim(file)
class(file)


#DESeq2
countdata = as.matrix(file[-1])
dim(countdata)
colnames(countdata)
rownames(countdata)

condition = factor(c(rep("tumor",7), rep("normal",5)))
condition
coldata = data.frame(row.names = colnames(countdata), condition)
coldata
ddsFull = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =~condition)
dds = DESeq(ddsFull)
res <- results( dds )
summary(res)
res_ordered = res[order(res$padj),]


#Upregulated and downregulated identification
res_d = as.data.frame(res_ordered)

res_d$diffexpressed <- "No change"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "Upregulated" 
res_d$diffexpressed[res_d$log2FoldChange > 0.6 & res_d$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "Downregulated"
res_d$diffexpressed[res_d$log2FoldChange < -0.6 & res_d$pvalue < 0.05] <- "Downregulated"


res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
colnames(res_d)[1] <- "genes"

#Export results to csv files 
a = res_d[which(res_d$diffexpressed=='Upregulated'),]
b = res_d[which(res_d$diffexpressed=='Downregulated'),]
write.table(a, file = 'gse32424_ESCC_EMT_up.tsv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(b, file = 'gse32424_ESCC_EMT_down.tsv', sep = '\t', row.names=FALSE, col.names=TRUE)


