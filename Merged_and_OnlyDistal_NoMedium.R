#Running only the high and low samples from proximal and distal colon together to see how it affects the results

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
suppressPackageStartupMessages(library(msigdbr))

#Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_ProxDist_merged/")

#-----Read in the raw counts
proxCounts <- read.csv("../Counts/Walnut_Proximal_Counts.csv", header = TRUE, sep = ",")
distCounts <- read.csv("../Counts/Walnut_Distal_Counts.csv", header = TRUE, sep = ",")
counts <- merge(proxCounts, distCounts, by = "X", all = TRUE)

#-----Read in the metadata
meta <- read.csv("../Metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
meta <- meta[meta$Sample != "WS47.1",]
metanomed <- meta[meta$URoAClass != "Med",]

#-----Find the samples that are medium producers
medProducers <- meta[meta$URoAClass == "Med",]$Sample
matchingColumns <- which(colnames(counts) %in% medProducers)


#-----Filter these samples out of the counts
countsnomed <- counts[,-matchingColumns]

#-----Filter out the outlier
outlier <- "WS47.1"
outlierCol <- which(colnames(countsnomed) %in% outlier)
countsnomed <- countsnomed[,-outlierCol]

#-----Set rownames of counts and organize the metadata in decreasing order, match in counts
countsnomed$Ensembl <- mapIds(org.Hs.eg.db, key = countsnomed$X, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first")
countsnomed <- countsnomed[!duplicated(countsnomed$Ensembl),]
countsnomed <- countsnomed[!grepl("^IG", countsnomed$X),]
countsnomed <- countsnomed[!grepl("^RPL", countsnomed$X),]
countsnomed <- countsnomed[!grepl("^RPS", countsnomed$X),]
rownames(countsnomed) <- paste(countsnomed$Ensembl, countsnomed$X, sep = " - ")

#-----Order the design file by Urolithin status
meta.ordered <- metanomed %>%
  arrange(URoAClass)
rownames(meta.ordered) <- meta.ordered$Sample
write.csv(meta.ordered, file = "ProximalAndDistal_NoMedium_NoOutlierWS47.1_Design.csv")

#-----We need to ensure that the samples in the counts file are in the same order as the design file
sampleOrder <- meta.ordered$Sample
countsnomed <- countsnomed[,sampleOrder]
write.csv(countsnomed, file = "ProximalAndDistal_NoMedium_NoOutlierWS47.1_Counts.csv")

#-----Sanity Check: Verify colnames(rawProx) and the same order as design
all(colnames(countsnomed) %in% rownames(meta.ordered))
colnames(countsnomed)[!(colnames(countsnomed) %in% rownames(meta.ordered))]
all(colnames(countsnomed) == rownames(meta.ordered))

#-----Filter
minimumCountpergene <- 10
MinSampleWithminimumgeneCounts <- 5

#-----Filter out low read counts 
countsnomed <- countsnomed[rowSums(data.frame(countsnomed>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#-----Set up DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(countsnomed),
                                  colData = meta.ordered,
                                  design = ~ URoAClass)

#-----Relevel
dds$URoAClass <- relevel(dds$URoAClass, ref = "Low")

#-----Run DESeq2
dds <- DESeq(dds)
write_rds(dds, file = "ProximalAndDistalCombined_NoMedium_dds.rds")

#-----Let's check out the dispersion plot for proximal
plotDispEsts(dds, main = "Distal Dispersion Estimates")

#-----Variance-stabilized transformation
vsd <- vst(dds)

#-----Plot PCA
PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(meta.ordered)), size = 3, max.overlaps = Inf) +
  ggtitle("Proximal and Distal Colon Merged, No Medium Producers") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("DistalPCA.pdf", PCA, width = 12, height = 8)

#-----Let's look at the comparison between high and low filtered
hl <- as.data.frame(results(dds, name = "URoAClass_High_vs_Low"))
hl <- hl[order(hl$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(hl))
hl$Symbols <- labels
write.csv(hl, file = "ProxDist_High_vs_Low_all_DEGs.csv")
hl.sig <- hl[hl$padj < 0.05,]
hl.sig <- na.omit(hl.sig)
hl.sigGenes <- rownames(hl.sig)
hl.fc <- hl[abs(hl$log2FoldChange) > 2,]

#-----Function to prepare DEGs for GSEA
getEntrez <- function(degs) {
  res.ordered <- degs
  res.ordered$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(res.ordered))
  res.ordered$Entrez <- mapIds(org.Hs.eg.db, key = res.ordered$Ensembl,
                               column = "ENTREZID", keytype = "ENSEMBL",
                               multiVals = "first")
  res.ordered.genes <- res.ordered$log2FoldChange
  
  #Assign Entrez IDs as names for the genes
  names(res.ordered.genes) <- res.ordered$Entrez
  
  #Remove duplicated Entrez IDs and their corresponding values
  unique_entrez_genes <- names(res.ordered.genes[!duplicated(names(res.ordered.genes))])
  unique_genes <- res.ordered.genes[unique_entrez_genes]
  unique_genes <- sort(unique_genes, decreasing = TRUE)
  
  return(unique_genes)
}

#-----Load in the genesets we want to look at 
m_t2gC2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene, gs_subcat) %>%
  filter(gs_subcat %in% c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"))

m_t2gC5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene,  gs_subcat) %>%
  filter(gs_subcat %in% c("GO", "GO:BP", "GO:CC", "GO:MF"))

#-----Combine into one large dataframe of gene sets
pathways <- rbind(m_t2gC2, m_t2gC5)

#-----Run the function on the three comparisons
hlGenes <- getEntrez(hl)


#-----Set a random seed 
set.seed(03061999)

#-----GSEA on the distal High vs Low comparison
hlGsea <- GSEA(hlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
hlGseaOb <- setReadable(hlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(hlGseaOb, "ProxDist_High_vs_Low_GSEA_results.rds")
hlGsea.res <- as.data.frame(setReadable(hlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(hlGsea.res, file = "ProximalDistal_High_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in High", suppressed = "Enriched in Low")
)
hlGsea.dotplot <- dotplot(hlGseaOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 40, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Proximal and Distal High vs Low")




#-----Read in the distal counts
distal <- read.csv("../Counts/Walnut_Distal_Counts.csv", header = TRUE, sep = ",")
meta <- read.csv("../Metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
meta <- meta[meta$Location == "Distal",]
medProducers <- meta[meta$URoAClass == "Med",]$Sample
matchingColumns <- which(colnames(distal) %in% medProducers)


#-----Filter these samples out of the counts
counts <- distal[,-matchingColumns]
meta <- meta[meta$URoAClass != "Med",]

#-----Set rownames of counts and organize the metadata in decreasing order, match in counts
counts$Ensembl <- mapIds(org.Hs.eg.db, key = counts$X, column = "ENSEMBL",
                              keytype = "SYMBOL", multiVals = "first")
counts <- counts[!duplicated(counts$Ensembl),]
counts <- counts[!grepl("^IG", counts$X),]
counts <- counts[!grepl("^RPL", counts$X),]
counts <- counts[!grepl("^RPS", counts$X),]
rownames(counts) <- paste(counts$Ensembl, counts$X, sep = " - ")
write.csv(counts, file = "Distal_Counts_NoMed.csv")

#-----Order the design file by Urolithin status
meta.ordered <- meta %>%
  arrange(URoAClass)
rownames(meta.ordered) <- meta.ordered$Sample
write.csv(meta.ordered, file = "Distal_NoMedium_NoOutlier_Design.csv")

#-----We need to ensure that the samples in the counts file are in the same order as the design file
sampleOrder <- meta.ordered$Sample
counts <- counts[,sampleOrder]

#-----Sanity Check: Verify colnames(rawProx) and the same order as design
all(colnames(count) %in% rownames(meta.ordered))
colnames(counts)[!(colnames(counts) %in% rownames(meta.ordered))]
all(colnames(counts) == rownames(meta.ordered))

#-----Set up DESeq2 object
distdds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = meta.ordered,
                                  design = ~ URoAClass)

#-----Relevel
distdds$URoAClass <- relevel(distdds$URoAClass, ref = "Low")

#-----Run DESeq2
distdds <- DESeq(distdds)
write_rds(distdds, "Distal_NoMed_dds.rds")
distdds <- readRDS("Distal_NoMed_dds.rds")

#-----Let's check out the dispersion plot for proximal
plotDispEsts(distdds, main = "Distal Dispersion Estimates")

#-----Variance-stabilized transformation
vsd <- vst(distdds)

#-----Plot PCA
PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(meta.ordered)), size = 3, max.overlaps = Inf) +
  ggtitle("Distal Colon No Medium Producers") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("Distal_NoMedium_High_vs_Low_PCA.tiff", PCA, dpi = 800)


#-----Let's look at the comparison between high and low filtered
dhl <- as.data.frame(results(distdds, name = "URoAClass_High_vs_Low"))
dhl <- dhl[order(dhl$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(dhl))
dhl$Symbols <- labels
write.csv(dhl, file = "Distal_NoMedium_High_vs_Low_all_DEGs.csv")
dhl.sig <- dhl[dhl$padj < 0.05,]
dhl.sig <- na.omit(dhl.sig)
dhl.sigGenes <- rownames(dhl.sig)
dhl.fc <- dhl[abs(dhl$log2FoldChange) > 2,]

showGenes <- c("ENSG00000099866 - MADCAM1", "ENSG00000171316 - CHD7")
MADCAM1counts <- counts[counts$Symbols == " MADCAM1",]
write.csv(MADCAM1counts, file = "Distal_High_vs_Low_MADCAM1_Counts.csv")

MADCAM1counts <- plotCounts(distdds, gene = "ENSG00000099866 - MADCAM1", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "MADCAM1",
           xlab = "URoA Group", returnData = TRUE)
MADCAM1plot <- ggplot(MADCAM1counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("MADCAM1") +
  theme(plot.title = element_text(size = 26))
MADCAM1plot
ggsave("Distal_High_vs_Low_MADCAM1plot.tiff", MADCAM1plot, dpi = 800)

MAML1counts <- plotCounts(distdds, gene = "ENSG00000161021 - MAML1", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "MAML1",
                            xlab = "URoA Group", returnData = TRUE)
MAML1plot <- ggplot(MAML1counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("MAML1") +
  theme(plot.title = element_text(size = 26))
MAML1plot
ggsave("Distal_High_vs_Low_MAML1plot.tiff", MAML1plot, dpi = 800)

CHD7counts <- plotCounts(distdds, gene = "ENSG00000171316 - CHD7", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "CHD7",
                          xlab = "URoA Group", returnData = TRUE)
CHD7plot <- ggplot(CHD7counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("CHD7") +
  theme(plot.title = element_text(size = 26))
CHD7plot
ggsave("Distal_High_vs_Low_CHD7plot.tiff", CHD7plot, dpi = 800)



distCountsPlotHL <- degPlotWide(distdds, rownames(distdds[showGenes,]), group="URoAClass")
distCountsPlotHL <- distCountsPlotHL + ggtitle("Distal High vs Low: Significant Genes")

counts$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(counts))

#-----Run the function on the three comparisons
dhlGenes <- getEntrez(dhl)


#-----Set a random seed 
set.seed(03061999)

#-----GSEA on the distal High vs Low comparison
dhlGsea <- GSEA(dhlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
dhlGseaOb <- setReadable(dhlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(dhlGseaOb, "Dist_NoMedium_High_vs_Low_GSEA_results.rds")
dhlGsea.res <- as.data.frame(setReadable(dhlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(dhlGsea.res, file = "Distal_NoMedium_High_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in High", suppressed = "Enriched in Low")
)
dhlGsea.dotplot <- dotplot(dhlGseaOb, x = "GeneRatio", color = "p.adjust", 
                          showCategory = 30, 
                          label_format = 50, 
                          split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Distal High vs Low No Medium")
dhl.dotplot <- dhlGsea.dotplot + theme(strip.text = element_text(size = 12))
ggsave("/users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/Distal_High_vs_Low_NoMedium.tiff", dhl.dotplot, dpi = 800, width = 24, height = 12)



dhlGSEA <- read_rds("/users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/Dist_NoMedium_High_vs_Low_GSEA_results.rds")


heatplot(dhlGseaOb, showCategory = c("GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_OXIDASE",
                                     "REACTOME_PEROXISOMAL_LIPID_METABOLISM",
                                     "WP_ESTROGEN_METABOLISM",
                                     "GOBP_ESTROGEN_METABOLIC_PROCESS",
                                     "WP_CODEINE_AND_MORPHINE_METABOLISM"),
                                     foldChange = dhlGenes)

heatplot(dhlGseaOb, showCategory = 5,
         foldChange = dhlGenes)



UTG1A3counts <- plotCounts(distdds, gene = "ENSG00000288702 - UGT1A3", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "UTG1A3",
                         xlab = "URoA Group", returnData = TRUE)
UTG1A3plot <- ggplot(UTG1A3counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("UTG1A3") +
  theme(plot.title = element_text(size = 26))
UTG1A3plot
ggsave("Distal_High_vs_Low_CHD7plot.tiff", CHD7plot, dpi = 800)


estrogenGenes <- c("ENSG00000137869 - CYP19A1",
                   "ENSG00000140465 - CYP1A1",
                   "ENSG00000100197 - CYP2D6",
                   "ENSG00000160868 - CYP3A4",
                   "ENSG00000106258 - CYP3A5",
                   "ENSG00000160870 - CYP3A7",
                   "ENSG00000278535 - DHRS11",
                   "ENSG00000198189 - HSD17B11",
                   "ENSG00000086696 - HSD17B2",
                   "ENSG00000203857 - HSD3B1",
                   "ENSG00000080511 - RDH8",
                   "ENSG00000196502 - SULT1A1",
                   "ENSG00000109193 - SULT1E1",
                   "ENSG00000241635 - UGT1A1",
                   "ENSG00000288702 - UGT1A3",
                   "ENSG00000196620 - UGT2B15",
                   "ENSG00000156096 - UGT2B4")

#Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/")
#Make a counts plot for each of the genes in the estrogen metabolism geneSet
for(i in estrogenGenes) {
  counts <- plotCounts(distdds, gene = i, intgroup = "URoAClass", normalized = TRUE, transform = FALSE,
                             xlab = "URoA Group", returnData = TRUE)
  countsPlot <- ggplot(counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
    geom_boxplot(width = 0.2) +
    geom_jitter(width = 0.09) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20)) +
    guides(fill = FALSE) +
    ggtitle(i) +
    theme(plot.title = element_text(size = 26))
  ggsave(paste(i,"Distal_High_vs_Low.tiff",sep = "_"), countsPlot, dpi = 600)
}

counts <- as.data.frame(counts(distdds, normalized = TRUE))
counts$Gene <- rownames(counts)
countsEstrogenSubset <- counts[counts$Gene %in% estrogenGenes,]
countsEstrogenSubset$Gene <- NULL

UrolithinAClass <- meta.ordered$URoAClass
names(UrolithinAClass) <- meta.ordered$Sample 
UrolithinAClass <- as.data.frame(UrolithinAClass)

pheatmap(as.matrix(countsEstrogenSubset),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 12,
         annotation_col = UrolithinAClass)


estrogenGenes2 <- c("ENSG00000100197 - CYP2D6",
                   "ENSG00000160868 - CYP3A4",
                   "ENSG00000106258 - CYP3A5",
                   "ENSG00000160870 - CYP3A7",
                   "ENSG00000278535 - DHRS11",
                   "ENSG00000198189 - HSD17B11",
                   "ENSG00000086696 - HSD17B2",
                   "ENSG00000203857 - HSD3B1",
                   "ENSG00000080511 - RDH8",
                   "ENSG00000196502 - SULT1A1",
                   "ENSG00000109193 - SULT1E1",
                   "ENSG00000241635 - UGT1A1",
                   "ENSG00000288702 - UGT1A3",
                   "ENSG00000196620 - UGT2B15")

counts <- as.data.frame(counts(distdds, normalized = TRUE))
counts$Gene <- rownames(counts)
countsEstrogenSubset <- counts[counts$Gene %in% estrogenGenes2,]
countsEstrogenSubset$Gene <- NULL

UrolithinAClass <- meta.ordered$URoAClass
names(UrolithinAClass) <- meta.ordered$Sample 
UrolithinAClass <- as.data.frame(UrolithinAClass)

pheatmap(as.matrix(countsEstrogenSubset),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = 12,
         annotation_col = UrolithinAClass)
















