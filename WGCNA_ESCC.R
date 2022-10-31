#Author: K.T.Shreya Parthasarathi
#Date: 25/03/2022
#Purpose: Estimate co-expressed ion channels, lipid metabolism genes and EMT-related genes from rna-seq datasets of patients with ESCC (WGCNA)

setwd("path\\to\\current\\working\\directory")

#Import libraries
library("DESeq2")
library("ggplot2")

#Data import and preprocessing
file = read.table("non-redundant_genes_GSE32424.tsv", sep = '\t',header = TRUE)
dim(file)

rownames(file) = file$hgnc_symbol
dim(file)
class(file)


#DESeq2 for data normalization
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

vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)

expr_normalized <- wpn_vsd

expr_normalized[1:2,1:12]

dim(expr_normalized)


input_mat = t(expr_normalized)
input_mat[1:5,1:10]


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#Also downloaded "impute" and "preprocessor"
#BiocManager::install("WGCNA")
library(WGCNA)
library(flashClust)

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence - ESCC"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - ESCC"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

softPower1 = 18
adjacency1 = adjacency(input_mat, power = softPower1, type = "signed")
dissTOM1   = 1-TOMsimilarity(adjacency1, TOMType="signed")
geneTree1  = flashClust(as.dist(dissTOM1), method="average")

par(mfrow=c(1,1))
minModuleSize = 50
dynamicMods1 = cutreeDynamic(dendro = geneTree1, distM = dissTOM1,
                             deepSplit = 0, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicMods1
table(dynamicMods1)
dynamicColors1 = labels2colors(dynamicMods1)
colors1 = table(dynamicColors1)
colors1

plotDendroAndColors(geneTree1, dynamicColors1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "ESCC - KMIO - No.of Samples = 12 pairs")

par(mar=c(5.1, 4.1, 4.1, 2.1))
traitdata1 = read.csv("GSE32424_sample_traits.txt", sep = '\t',header = TRUE)
traitdata1
patient1 = rownames(input_mat)
patient1
traitRows1 = match(patient1, traitdata1$Sample_ID)
traitRows1
datTraits1 = traitdata1[traitRows1, -1]
names(datTraits1)
nGenes1 = ncol(input_mat)
nSamples1 = nrow(input_mat)
MEs0 = moduleEigengenes(input_mat, dynamicColors1)$eigengenes
MEs1 = orderMEs(MEs0)
MEs1
moduleTraitCor1 = cor(MEs1, datTraits1, use = "p")
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples1)
textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(",signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1) = dim(moduleTraitCor1)
labeledHeatmap(Matrix = moduleTraitCor1,
               xLabels = names(datTraits1),
               yLabels = names(MEs1),
               ySymbols = NULL,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix1,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1), 
               main = paste("Module- Trait relationships - KMIO"), ylab = "Gene Expression-Based Modules")

par(mfrow = c(1,1))
tumor = as.data.frame(datTraits1$Tumor)
names(tumor) = "tumor"
modNames1 = substring(names(MEs1), 3)
geneModuleMembership1 = as.data.frame(cor(input_mat, MEs1, use = "p"))
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples1))
names(geneModuleMembership1) = paste("MM", modNames1, sep="")
names(MMPvalue1) = paste("p.MM", modNames1, sep="")
geneTraitSignificance1 = as.data.frame(cor(input_mat, tumor, use = "p"))
GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples1))
names(geneTraitSignificance1) = paste("GS.", names(tumor), sep="")
names(GSPvalue1) = paste("p.GS.", names(tumor), sep="")
module = "blue"
column = match(module, modNames1)
column
moduleGenes = dynamicColors1==module
moduleGenes
verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),
                   abs(geneTraitSignificance1[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Tumor",
                   main = paste("GSE32424\nModule membership vs. gene significance\n"),
                   cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, col = module)


#Import to cytoscape
TOM1 = TOMsimilarityFromExpr(input_mat, power = 18)
modules = c("blue")
genes = colnames(input_mat)
genes
inModule = is.finite(match(dynamicColors1, modules));
modGenes = genes[inModule];
modTOM = TOM1[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("GSE32424_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("GSE32424_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.0,
                               nodeNames = modGenes,
                              nodeAttr = dynamicColors1[inModule])









