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
suppressPackageStartupMessages(library(msigdbr))
#------------------------------------------------------------------------------#

#-----Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_19th_2023_Distal_No_IGV.MIT.RIB/")

#-----Read in the raw counts for the proximal colon
rawDist <- read.csv("counts/Walnut_Distal_Counts.csv", header = TRUE, sep = ",")

#-----Map ENSEMBL Ids to gene symbols, remove duplicated genes, set rownames, remove X and Ensembl columns
rawDist$Ensembl <- mapIds(org.Hs.eg.db, key = rawDist$X, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first")
rawDist <- rawDist[!duplicated(rawDist$Ensembl),]
rawDist <- rawDist[!grepl("^IG", rawDist$X),]
rawDist <- rawDist[!grepl("^RPL", rawDist$X),]
rawDist <- rawDist[!grepl("^RPS", rawDist$X),]
rownames(rawDist) <- paste(rawDist$Ensembl, rawDist$X, sep = " - ")
rawDist$X <- NULL
rawDist$Ensembl <- NULL

#-----Design file: we only want to take the rows that are metadata for the proximal colon
design <- read.csv("metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
designDist <- design[design$Location == "Distal",]
designDist$URoAClass <- factor(designDist$URoAClass, levels = c("Low", "Med", "High"))

#-----Order the design file by Urolithin status
designDist.ordered <- designDist %>%
  arrange(URoAClass)
rownames(designDist.ordered) <- designDist.ordered$Sample

#-----We need to ensure that the samples in the counts file are in the same order as the design file
distSampleOrder <- designDist.ordered$Sample
rawDist <- rawDist[,distSampleOrder]

#-----Sanity Check: Verify colnames(rawProx) and the same order as design
all(colnames(rawDist) %in% rownames(designDist.ordered))
colnames(rawDist)[!(colnames(raw) %in% rownames(designDist.ordered))]
all(colnames(rawDist) == rownames(designDist.ordered))

#-----Filter
minimumCountpergene <- 10
MinSampleWithminimumgeneCounts <- 5

#-----Filter out low read counts 
rawDist <- rawDist[rowSums(data.frame(rawDist>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]

#-----Set up DESeq2 object
distdds <- DESeqDataSetFromMatrix(countData = round(rawDist),
                                  colData = designDist.ordered,
                                  design = ~ URoAClass)

#-----Relevel
distdds$URoAClass <- relevel(distdds$URoAClass, ref = "Low")

#-----Run DESeq2
distdds <- DESeq(distdds)
write_rds(distdds, "Distal_dds.rds")
distdds <- readRDS("Distal_dds.rds")

#-----Let's check out the dispersion plot for proximal
plotDispEsts(distdds, main = "Distal Dispersion Estimates")

#-----Variance-stabilized transformation
vsd <- vst(distdds)

#-----Plot PCA
PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(designDist.ordered)), size = 3, max.overlaps = Inf) +
  ggtitle("Distal Colon") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("DistalPCA.pdf", PCA, width = 12, height = 8)

#-----Let's look at the comparison between high and low filtered
dhl <- as.data.frame(results(distdds, name = "URoAClass_High_vs_Low"))
dhl <- dhl[order(dhl$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(dhl))
dhl$Symbols <- labels
write.csv(dhl, file = "Distal_High_vs_Low_all_DEGs.csv")
dhl.sig <- dhl[dhl$padj < 0.05,]
dhl.sigGenes <- rownames(dhl.sig)
dhl.fc <- dhl[abs(dhl$log2FoldChange) > 2,]

#-----Gene counts plot for the significant genes
distCountsPlotHL <- degPlotWide(distdds, rownames(distdds[dhl.sigGenes,]), group="URoAClass")
distCountsPlotHL <- distCountsPlotHL + ggtitle("Distal High vs Low: Significant Genes")
ggsave("Distal_High_Vs_Low_SignificantGene_Plots.pdf", distCountsPlotHL)

#-----Volcano plot
dhl.volcano <- EnhancedVolcano(dhl,
                               lab = dhl$Symbols,
                               title = "Distal High vs Low",
                               subtitle = "",
                               legendPosition = "bottom",
                               x = 'log2FoldChange',
                               y = 'padj') 
ggsave("Distal_High_vs_Low_VolcanoPlot.pdf", dhl.volcano, width = 12, height = 8)

#-----Let's look at the comparison between med and low filtered
dml <- as.data.frame(results(distdds, name = "URoAClass_Med_vs_Low"))
dml <- dml[order(dml$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(dml))
dml$Symbols <- labels
write.csv(dhl, file = "Distal_Med_vs_Low_all_DEGs.csv")
dml.sig <- dml[dml$padj < 0.05,]
dml.sigGenes <- rownames(dml.sig)
dml.fc <- dml[abs(dml$log2FoldChange) > 3,]
dml.fc.genes <- rownames(dml.fc)

#-----Volcano plot
dml.volcano <- EnhancedVolcano(dml,
                               pCutoff = 0.05,
                               lab = dml$Symbols,
                               title = "Distal Med vs Low",
                               subtitle = "",
                               legendPosition = "bottom",
                               x = 'log2FoldChange',
                               y = 'pvalue') 
ggsave("Distal_Med_vs_Low_VolcanoPlot.pdf", dml.volcano, width = 20, height = 8)

#-----Gene plots for MOST significant genes with fold changes > abs(3)
distCountsPlotML <- degPlotWide(distdds, rownames(distdds[dml.fc.genes,]), group="URoAClass")
distCountsPlotML <- distCountsPlotML + ggtitle("Distal Med vs Low: Top Significant Genes: abs(l2FC) > 3")
ggsave("Distal_Med_Vs_Low_SignificantGene_Plots.pdf", distCountsPlotML)

#-----Prepare a heatmap
dmlCounts <- as.data.frame(counts(distdds, normalized = TRUE))
dmlTopCounts <- dmlCounts[dml.fc.genes,]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(dmlTopCounts))
rownames(dmlTopCounts) <- labels 
URoAClass <- designDist.ordered$URoAClass
names(URoAClass) <- designDist.ordered$Sample
URoAClass <- as.data.frame(URoAClass)
dmlFC <- dml.fc$log2FoldChange
names(dmlFC) <- rownames(dmlTopCounts)
dmlFC <- as.data.frame(dmlFC)

#------Print heatmap
dml.heatmap <- pheatmap(as.matrix(dmlTopCounts),
                        cluster_rows = TRUE, 
                        cluster_cols = TRUE, 
                        scale = "row",
                        annotation_col = URoAClass,
                        annotation_row = dmlFC)

#-----Let's do the comparison between Med and High
dmh <- as.data.frame(results(distdds, contrast = c("URoAClass", "Med", "High")))
dmh <- dmh[order(dmh$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(dmh))
dmh$Symbols <- labels
write.csv(dmh, file = "Distal_Med_vs_High_all_DEGs.csv")
dmh.sig <- dmh[dmh$padj < 0.05,]
dmh.sigGenes <- rownames(dmh.sig)


#-----Gene plots for MOST significant genes with fold changes > abs(3)
distCountsPlotMH <- degPlotWide(distdds, rownames(distdds[dmh.sigGenes,]), group="URoAClass")
distCountsPlotMH <- distCountsPlotMH + ggtitle("Distal Med vs High Significant Genes")
ggsave("Distal_Med_Vs_High_SignificantGene_Plots.pdf", distCountsPlotMH)

#-----Volcano plot
dmh.volcano <- EnhancedVolcano(dmh,
                               lab = dmh$Symbols,
                               title = "Distal High vs Low",
                               subtitle = "",
                               legendPosition = "bottom",
                               x = 'log2FoldChange',
                               y = 'padj') 
ggsave("Distal_Med_vs_High_VolcanoPlot.pdf", dmh.volcano, width = 12, height = 8)

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
dhlGenes <- getEntrez(dhl)
dmlGenes <- getEntrez(dml)
dmhGenes <- getEntrez(dmh)

#-----Set a random seed 
set.seed(03061999)

#-----GSEA on the distal High vs Low comparison
dhlGsea <- GSEA(dhlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
dhlGseaOb <- setReadable(dhlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(dhlGseaOb, "Distal_High_vs_Low_GSEA_results.rds")
dhlGsea.res <- as.data.frame(setReadable(dhlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(dhlGsea.res, file = "Distal_High_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in High", suppressed = "Enriched in Low")
)
dhlGsea.dotplot <- dotplot(dhlGseaOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Distal High vs Low")

#-----GSEA on the distal Med vs Low comparison
dmlGsea <- GSEA(dmlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
dmlGseaOb <- setReadable(dmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(dmlGseaOb, "Distal_Med_vs_Low_GSEA_results.rds")
dmlGsea.res <- as.data.frame(setReadable(dmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(dmlGsea.res, file = "Distal_Med_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in Low")
)
dmlGsea.dotplot <- dotplot(dmlGseaOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Distal Med vs Low")

#-----GSEA on the distal Med vs High comparison
dmhGsea <- GSEA(dmhGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
dmhGseaOb <- setReadable(dmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(dmhGseaOb, "Distal_Med_vs_High_GSEA_results.rds")
dmhGsea.res <- as.data.frame(setReadable(dmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(dmhGsea.res, file = "Distal_Med_vs_High_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in High")
)
dmhGsea.dotplot <- dotplot(dmhGseaOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Distal Med vs High")








