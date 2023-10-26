#Walnut RNASeq analysis
#Set working directory
setwd("/Users/mikemartinez/Desktop/Walnut_RNASeq/Proximal_vs_distal/")

suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

#Differential gene expression
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(DEGreport))

#Graphics and visualizations
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(patchwork))

#GSEA analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
#------------------------------------------------------------------------------#





#-----Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/")

#-----Read in the raw counts for the proximal colon
rawProx <- read.csv("Walnut_Proximal_Counts.csv", header = TRUE, sep = ",")

#-----Map ENSEMBL Ids to gene symbols, remove duplicated genes, set rownames, remove X and Ensembl columns
rawProx$Ensembl <- mapIds(org.Hs.eg.db, key = rawProx$X, column = "ENSEMBL",
                      keytype = "SYMBOL", multiVals = "first")
rawProx <- rawProx[!duplicated(rawProx$Ensembl),]
rownames(rawProx) <- paste(rawProx$Ensembl, rawProx$X, sep = " - ")
rawProx$X <- NULL
rawProx$Ensembl <- NULL

#-----Design file: we only want to take the rows that are metadata for the proximal colon
design <- read.csv("New_Walnut_Design.csv", header = TRUE, sep = ",")
designProx <- design[design$Location == "Proximal",]
designProx$URoAClass <- factor(designProx$URoAClass, levels = c("Low", "Med", "High"))

#-----Order the design file by Urolithin status
designProx.ordered <- designProx %>%
  arrange(URoAClass)
rownames(designProx.ordered) <- designProx.ordered$Sample

#-----We need to ensure that the samples in the counts file are in the same order as the design file
proxSampleOrder <- designProx.ordered$Sample
rawProx <- rawProx[,proxSampleOrder]

#-----Sanity Check: Verify colnames(rawProx) and the same order as design
all(colnames(raw) %in% rownames(design))
colnames(raw)[!(colnames(raw) %in% rownames(design))]
all(colnames(raw) == rownames(design))

#-----Low URolithinA is the first 13 columns of the counts
#---Let's set a 25$ detection rate for the genes of 10 counts to pass filtering

#Filter
minimumCountpergene <- 10
MinSampleWithminimumgeneCounts <- 5

#Filter out low read counts 
rawProx <- rawProx[rowSums(data.frame(rawProx>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

# #-----Set thresholds
# minimumCountpergeneLowUroA <- 10
# MinSampleWithminimumgeneCountsLowUroA <- 3
# 
# #-----Filter out low read counts for each Urolithin group. 
# lowUroAFilt <- rowSums(rawProx[,1:13] > minimumCountpergeneLowUroA) > MinSampleWithminimumgeneCountsLowUroA
# rawProx <- rawProx[lowUroAFilt,]
# medUroAFilt <- rowSums(rawProx[,14:27] > minimumCountpergeneLowUroA) > MinSampleWithminimumgeneCountsLowUroA
# rawProx <- rawProx[medUroAFilt,]
# highUroAFilt <- rowSums(rawProx[,27:ncol(rawProx)] > minimumCountpergeneLowUroA) > MinSampleWithminimumgeneCountsLowUroA
# rawProx <- rawProx[highUroAFilt,]

#-----Set up DESeq2 object
proxdds <- DESeqDataSetFromMatrix(countData = round(rawProx),
                              colData = designProx.ordered,
                              design = ~ URoAClass)

#-----Relevel
proxdds$URoAClass <- relevel(proxdds$URoAClass, ref = "Low")

#-----Run DESeq2
proxdds <- DESeq(proxdds)
write_rds(proxdds, "Proximal_dds.rds")

#-----Let's check out the dispersion plot for proximal
plotDispEsts(proxdds, main = "Proximal Dispersion Estimates")

#-----Variance-stabilized transformation
vsd <- vst(proxdds)

#-----Plot PCA
PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(designProx.ordered)), size = 3, max.overlaps = Inf) +
  ggtitle("Proximal Colon") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("ProximalPCA.pdf", PCA, width = 12, height = 8)

#-----Let's look at the results for high vs low
phl <- as.data.frame(results(proxdds, name = "URoAClass_High_vs_Low"))
phl <- phl[order(phl$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(phl))
phl$Symbols <- labels
write.csv(phl, file = "Proximal_High_vs_Low_all_DEGs.csv")
phl.sig <- phl[phl$padj < 0.05,]
phl.sigGenes <- rownames(phl.sig)

#-----Gene counts plot for the significant genes
proxCountsPlotHL <- degPlotWide(proxdds, rownames(proxdds[phl.sigGenes,]), group="URoAClass")
proxCountsPlotHL <- proxCountsPlot + ggtitle("Proximal High vs Low: Significant Genes")
ggsave("Proximal_High_Vs_Low_SignificantGene_Plots.pdf", proxCountsPlotHL)

#-----Volcano plot
phl.volcano <- EnhancedVolcano(phl,
                           lab = phl$Symbols,
                           title = "Proximal High vs Low",
                           subtitle = "",
                           legendPosition = "bottom",
                           x = 'log2FoldChange',
                           y = 'pvalue') 
ggsave("Proximal_High_vs_Low_VolcanoPlot.pdf", phl.volcano, width = 12, height = 8)

#-----Let's look at the results for medium vs low
pml <- as.data.frame(results(proxdds, name = "URoAClass_Med_vs_Low"))
pml <- pml[order(pml$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(pml))
pml$Symbols <- labels
write.csv(pml, file = "Proximal_Med_vs_Low_all_DEGs.csv")
pml.sig <- pml[pml$padj < 0.05,]

#-----Let's look at the results for high vs med
pmh <- as.data.frame(results(proxdds, contrast = c("URoAClass", "Med", "High")))
pmh <- pmh[order(pmh$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(pmh))
pmh$Symbols <- labels
write.csv(pmh, file = "Proximal_Med_vs_High_all_DEGs.csv")
pmh.sig <- pmh[pmh$padj < 0.05,]
pmh.sigGenes <- rownames(pmh.sig)

#-----Gene counts plot for the significant genes
proxCountsPlotMH <- degPlotWide(proxdds, rownames(proxdds[pmh.sigGenes,]), group="URoAClass")
proxCountsPlotMH <- proxCountsPlotMH + ggtitle("Proximal Med vs High: Significant Genes")
ggsave("Proximal_Med_Vs_High_SignificantGene_Plots.pdf", proxCountsPlotMH)

#-----Volcano plot
pmh.volcano <- EnhancedVolcano(pmh,
                               lab = pmh$Symbols,
                               title = "Proximal Med vs High",
                               subtitle = "",
                               legendPosition = "bottom",
                               x = 'log2FoldChange',
                               y = 'pvalue') 
ggsave("Proximal_Med_vs_High_VolcanoPlot.pdf", pmh.volcano, width = 12, height = 8)










