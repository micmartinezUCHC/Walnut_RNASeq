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
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/")

#-----Read in the raw counts for the proximal colon
rawProx <- read.csv("counts/Walnut_Proximal_Counts.csv", header = TRUE, sep = ",")

#-----Map ENSEMBL Ids to gene symbols, remove duplicated genes, set rownames, remove X and Ensembl columns
rawProx$Ensembl <- mapIds(org.Hs.eg.db, key = rawProx$X, column = "ENSEMBL",
                      keytype = "SYMBOL", multiVals = "first")
rawProx <- rawProx[!duplicated(rawProx$Ensembl),]
rownames(rawProx) <- paste(rawProx$Ensembl, rawProx$X, sep = " - ")
rawProx$X <- NULL
rawProx$Ensembl <- NULL

#-----Design file: we only want to take the rows that are metadata for the proximal colon
design <- read.csv("metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
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
all(colnames(rawProx) %in% rownames(designProx.ordered))
colnames(raw)[!(colnames(raw) %in% rownames(design))]
all(colnames(rawProx) == rownames(designProx.ordered))

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


#---------------------------#
#---------------------------#
#---------------------------#
#-----October 19th 2023-----#

#Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_18th_2023/")

#-----Read in the results files
phl <- read.csv("../October_18th_2023_Proximal_Analysis/Proximal_High_vs_Low/Proximal_High_vs_Low_all_DEGs.csv", header = TRUE, sep = ",", row.names = 1)
phl <- phl[order(phl$log2FoldChange, decreasing = TRUE),]
pml <- read.csv("Proximal_Med_vs_Low/Proximal_Med_vs_Low_all_DEGs.csv", header = TRUE, sep = ",", row.names = 1)
pml <- pml[order(phl$log2FoldChange, decreasing = TRUE),]
pmh <- read.csv("Proximal_Med_vs_High/Proximal_Med_vs_High_all_DEGs.csv", header = TRUE, sep = ",", row.names = 1)
pmh <- pmh[order(phl$log2FoldChange, decreasing = TRUE),]

#-----Let's write a function to get the entrez IDs and assembl the lists for GSEA
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

#-----Run the function on the three comparisons
phlGenes <- getEntrez(phl)
pmlGenes <- getEntrez(pml)
pmhGenes <- getEntrez(pmh)

#-----Load in the genesets we want to look at 
m_t2gC2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene, gs_subcat) %>%
  filter(gs_subcat %in% c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"))

m_t2gC5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene,  gs_subcat) %>%
  filter(gs_subcat %in% c("GO", "GO:BP", "GO:CC", "GO:MF"))

#-----Combine into one large dataframe of gene sets
pathways <- rbind(m_t2gC2, m_t2gC5)

#-----Set a random seed 
set.seed(03061999)

#-----Run GSEA
phlGsea <- GSEA(phlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
phlGseaOb <- setReadable(phlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(phlGseaOb, "Proximal_High_vs_Low_GSEA_results.rds")
phlGsea.res <- as.data.frame(setReadable(phlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(phlGsea.res, file = "Proximal_High_vs_Low_GSEA_results.csv")

phlGSEAOb <- readRDS("../October_18th_2023_Proximal_Analysis/Proximal_High_vs_Low/GSEA/Proximal_High_vs_Low_GSEA_results.rds")
custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in High", suppressed = "Enriched in Low")
)
phlGsea.dotplot <- dotplot(phlGSEAOb, x = "GeneRatio", color = "p.adjust", 
                                showCategory = 20, 
                                label_format = 50, 
                                split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Proximal High vs Low")

phlCNET <- cnetplot(phlGSEAOb, foldChange = phlGenes, colorEdge = TRUE,
                    node_label = "gene",
                    showCategory = c("GOMF_BIOACTIVE_LIPID_RECEPTOR_ACTIVITY",
                                     "BIOCARTA_INFLAM_PATHWAY",
                                     "BIOCARTA_IL17_PATHWAY"))

pmlGsea <- GSEA(pmlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
pmlGseaOb <- setReadable(pmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(pmlGseaOb, "Proximal_Med_vs_Low_GSEA_results.rds")
pmlGsea.res <- as.data.frame(setReadable(pmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(pmlGsea.res, file = "Proximal_Med_vs_Low_GSEA_results.csv")

pmlGSEAOb <- readRDS("../October_18th_2023_Proximal_Analysis/Proximal_Med_vs_Low/GSEA/Proximal_Med_vs_Low_GSEA_results.rds")
custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in Low")
)
pmlGsea.dotplot <- dotplot(pmlGSEAOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Proximal Med vs Low")

pmhGsea <- GSEA(pmhGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
pmhGseaOb <- setReadable(pmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(pmhGseaOb, "Proximal_Med_vs_High_GSEA_results.rds")
pmhGsea.res <- as.data.frame(setReadable(pmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(pmhGsea.res, file = "Proximal_Med_vs_High_GSEA_results.csv")

pmhGSEAOb <- readRDS("../October_18th_2023_Proximal_Analysis/Proximal_Med_vs_High/GSEA/Proximal_Med_vs_High_GSEA_results.rds")
custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in High")
)
pmhGsea.dotplot <- dotplot(pmhGSEAOb, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Proximal Med vs High")

overlap <- intersect(phlGsea.res$Description, pmlGsea.res$Description)
overlap <- intersect(overlap, pmhGsea.res$Description)

#-----Create a new dataframe containing just these terms and the enrichment scores
phlOverlap <- phlGsea.res[phlGsea.res$Description %in% overlap,]
pmlOverlap <- pmlGsea.res[pmlGsea.res$Description %in% overlap,]
pmhOverlap <- pmhGsea.res[pmhGsea.res$Description %in% overlap,]

#-----Collate information
sharedGeneSets <- data.frame(matrix(nrow = length(overlap)))
rownames(sharedGeneSets) <- overlap
sharedGeneSets$Prox.HighVLow.ES <- phlOverlap$enrichmentScore
sharedGeneSets$Prox.HighVLow.Padj <- phlOverlap$p.adjust
sharedGeneSets$Prox.HighVLow.qval <- phlOverlap$qvalue
sharedGeneSets$Prox.MedVLow.ES <- pmlOverlap$enrichmentScore
sharedGeneSets$Prox.MedVLow.Padj <- pmlOverlap$p.adjust
sharedGeneSets$Prox.MedVLow.qval <- pmlOverlap$qvalue
sharedGeneSets$Prox.MedVHigh.ES <- pmhOverlap$enrichmentScore
sharedGeneSets$Prox.MedVHigh.Padj <- pmhOverlap$p.adjust
sharedGeneSets$Prox.MedVHigh.qval <- pmhOverlap$qvalue
sharedGeneSets$matrix.nrow...length.overlap.. <- NULL

#-----Add a sign variable
sharedGeneSets$Prox.HighVLow.Sign <- ifelse(sharedGeneSets$Prox.HighVLow.ES > 0, "Enriched in High", "Enriched in Low")
sharedGeneSets$Prox.MedVLow.Sign <- ifelse(sharedGeneSets$Prox.MedVLow.ES > 0, "Enriched in Med", "Enriched in Low")
sharedGeneSets$Prox.MedVHigh.Sign <- ifelse(sharedGeneSets$Prox.MedVHigh.ES > 0, "Enriched in Med", "Enriched in High")

#-----Write as a csv
write.csv(sharedGeneSets, file = "SharedGeneSets_Between_3ComparisonGroups.csv")

#
 #
  #
   #
    #   Now let's try the same analysis, but removing IGV genes, mitochondrial genes, and ribosomal genes
   #
  #
 #
#

#-----Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_Proximal_No_Medium/")

#-----Read in the raw counts for the proximal colon
rawProx <- read.csv("../Counts/Walnut_Proximal_Counts.csv", header = TRUE, sep = ",")

#-----Map ENSEMBL Ids to gene symbols, remove duplicated genes, set rownames, remove X and Ensembl columns
rawProx$Ensembl <- mapIds(org.Hs.eg.db, key = rawProx$X, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first")
rawProx <- rawProx[!duplicated(rawProx$Ensembl),]
rawProx <- rawProx[!grepl("^IG", rawProx$X),]
rawProx <- rawProx[!grepl("^RPL", rawProx$X),]
rawProx <- rawProx[!grepl("^RPS", rawProx$X),]
rownames(rawProx) <- paste(rawProx$Ensembl, rawProx$X, sep = " - ")
rawProx$X <- NULL
rawProx$Ensembl <- NULL

#-----Design file: we only want to take the rows that are metadata for the proximal colon
design <- read.csv("../Metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
designProx <- design[design$Location == "Proximal",]
designProx$URoAClass <- factor(designProx$URoAClass, levels = c("Low", "Med", "High"))

#-----Let's remove that one outlier from the design matrix
designProx <- designProx[designProx$Sample != "WS47.1",]

#-----Order the design file by Urolithin status
designProx.ordered <- designProx %>%
  arrange(URoAClass)
rownames(designProx.ordered) <- designProx.ordered$Sample
write.csv(designProx.ordered, file = "ProximalWalnut_Design.csv")

#-----We need to ensure that the samples in the counts file are in the same order as the design file
proxSampleOrder <- designProx.ordered$Sample
rawProx <- rawProx[,proxSampleOrder]
write.csv(rawProx, file = "Ordered_rawProximal_Counts.csv")

#-----Let's try removing that one outlier samples
rawProx$WS47.1 <- NULL

#-----Sanity Check: Verify colnames(rawProx) and the same order as design
all(colnames(rawProx) %in% rownames(designProx.ordered))
colnames(raw)[!(colnames(raw) %in% rownames(design))]
all(colnames(rawProx) == rownames(designProx.ordered))

#Filter
minimumCountpergene <- 10
MinSampleWithminimumgeneCounts <- 5

#Filter out low read counts 
rawProx <- rawProx[rowSums(data.frame(rawProx>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]


#-----Set up DESeq2 object
proxdds <- DESeqDataSetFromMatrix(countData = round(rawProx),
                                  colData = designProx.ordered,
                                  design = ~ URoAClass)

#-----Relevel
proxdds$URoAClass <- relevel(proxdds$URoAClass, ref = "Low")

#-----Run DESeq2
proxddsFilt <- DESeq(proxdds)
write_rds(proxddsFilt, "Proximal_ddsFilt.rds")
proxddsFilt <- read_rds("Proximal_ddsFilt.rds")

#-----Visualize dispersion estimates
plotDispEsts(proxddsFilt, main = "Proximal (No IGVs, Ribosomal Genes) Dispersion Estimates")

#-----Visualize PCA
vsd <- vst(proxddsFilt)

#-----Plot PCA
PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(designProx.ordered)), size = 3, max.overlaps = Inf) +
  ggtitle("Proximal Colon: No IGVs, Ribosome Proteins, Outlier Removed") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("Filtered_ProximalPCA.pdf", PCA, width = 12, height = 8)


#-----Let's look at the comparison between high and low filtered
phl.filt <- as.data.frame(results(proxddsFilt, name = "URoAClass_High_vs_Low"))
phl.filt$Gene <- rownames(phl.filt)
counts <- as.data.frame(counts(proxddsFilt))
counts$Gene <- rownames(counts)
phl.filt <- merge(phl.filt, counts, by = "Gene", all = TRUE)
rownames(phl.filt) <- phl.filt$Gene
phl.filt <- phl.filt[order(phl.filt$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(phl.filt))
phl.filt$Symbols <- labels
phl.filt$Gene <- NULL
write.csv(phl.filt, file = "Filtered_Proximal_High_vs_Low_all_DEGs.csv")
phl.filt.sig <- phl.filt[phl.filt$padj < 0.05,]
phl.filt.sigGenes <- rownames(phl.filt.sig)
phl.fc <- phl.filt[abs(phl.filt$log2FoldChange) > 2,]

#-----Gene counts plot for the significant genes
filtproxCountsPlotHL <- degPlotWide(proxddsFilt, rownames(proxddsFilt[phl.filt.sigGenes,]), group="URoAClass")
filtproxCountsPlotHL <- filtproxCountsPlotHL + ggtitle("Filtered Proximal High vs Low: Significant Genes")
ggsave("Filtered Proximal_High_Vs_Low_SignificantGene_Plots.pdf", filtproxCountsPlotHL)

#-----Volcano plot
filtphl.volcano <- EnhancedVolcano(phl.filt,
                               lab = phl.filt$Symbols,
                               title = "Filtered Proximal High vs Low",
                               subtitle = "",
                               legendPosition = "bottom",
                               x = 'log2FoldChange',
                               y = 'pvalue') 
ggsave("Filtered_Proximal_High_vs_Low_VolcanoPlot.pdf", filtphl.volcano, width = 12, height = 8)


#-----Let's look at the comparison between med and low filtered
pml.filt <- as.data.frame(results(proxddsFilt, name = "URoAClass_Med_vs_Low"))
pml.filt$Gene <- rownames(pml.filt)
pml.filt <- merge(pml.filt, counts, by ="Gene", all = TRUE)
pml.filt <- pml.filt[order(pml.filt$log2FoldChange, decreasing = TRUE),]
rownames(pml.filt) <- pml.filt$Gene
pml.filt$Gene <- NULL
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(pml.filt))
pml.filt$Symbols <- labels
write.csv(pml.filt, file = "Filtered_Proximal_Med_vs_Low_all_DEGs.csv")
pml.filt.sig <- pml.filt[pml.filt$padj < 0.05,]
pml.filt.sigGenes <- rownames(pml.filt.sig)
pml.fc <- phl.filt[abs(pml.filt$log2FoldChange) > 2,]

#-----Gene counts plot for the significant genes
filtproxCountsPlotML <- degPlotWide(proxddsFilt, rownames(proxddsFilt[pml.filt.sigGenes,]), group="URoAClass")
filtproxCountsPlotML <- filtproxCountsPlotML + ggtitle("Filtered Proximal Med vs Low: Significant Genes")
ggsave("Filtered Proximal_Med_Vs_Low_SignificantGene_Plots.pdf", filtproxCountsPlotML)

#-----Volcano plot
filtpml.volcano <- EnhancedVolcano(pml.filt,
                                   lab = pml.filt$Symbols,
                                   title = "Filtered Proximal Med vs Low",
                                   subtitle = "",
                                   legendPosition = "bottom",
                                   x = 'log2FoldChange',
                                   y = 'pvalue') 
ggsave("Filtered_Proximal_Med_vs_Low_VolcanoPlot.pdf", filtpml.volcano, width = 12, height = 8)

#-----Let's look at the comparison between med and low filtered
pmh.filt <- as.data.frame(results(proxddsFilt, contrast = c("URoAClass", "Med", "High")))
pmh.filt$Gene <- rownames(pmh.filt)
pmh.filt <- merge(pmh.filt, counts, by = "Gene", all = TRUE)
rownames(pmh.filt) <- pmh.filt$Gene
pmh.filt$Gene <- NULL
pmh.filt <- pmh.filt[order(pmh.filt$log2FoldChange, decreasing = TRUE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(pmh.filt))
pmh.filt$Symbols <- labels
write.csv(pmh.filt, file = "Filtered_Proximal_Med_vs_High_all_DEGs.csv")
pmh.filt.sig <- pmh.filt[pmh.filt$padj < 0.05,]
pmh.filt.sigGenes <- rownames(pmh.filt.sig)
pmh.fc <- pmh.filt[abs(pmh.filt$log2FoldChange) > 2,]

#-----Gene counts plot for the significant genes
filtproxCountsPlotMH <- degPlotWide(proxddsFilt, rownames(proxddsFilt[pmh.filt.sigGenes,]), group="URoAClass")
filtproxCountsPlotMH <- filtproxCountsPlotMH + ggtitle("Filtered Proximal Med vs High: Significant Genes")
ggsave("Filtered Proximal_Med_Vs_High_SignificantGene_Plots.pdf", filtproxCountsPlotMH)

#-----Volcano plot
filtpmh.volcano <- EnhancedVolcano(pmh.filt,
                                   lab = pml.filt$Symbols,
                                   title = "Filtered Proximal Med vs High",
                                   subtitle = "",
                                   legendPosition = "bottom",
                                   x = 'log2FoldChange',
                                   y = 'padj') 
ggsave("Filtered_Proximal_Med_vs_High_VolcanoPlot.pdf", filtpmh.volcano, width = 12, height = 8)

#-----Run GSEA
#-----Run the function on the three comparisons
filt.phlGenes <- getEntrez(phl.filt)
filt.pmlGenes <- getEntrez(pml.filt)
filt.pmhGenes <- getEntrez(pmh.filt)

#-----Set a random seed 
set.seed(03061999)

#-----Run GSEA
filt.phlGsea <- GSEA(filt.phlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
filt.phlGseaOb <- setReadable(filt.phlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(filt.phlGseaOb, "Filtered_Proximal_High_vs_Low_GSEA_results.rds")
filt.phlGsea.res <- as.data.frame(setReadable(filt.phlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(filt.phlGsea.res, file = "Filtered_Proximal_High_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in High", suppressed = "Enriched in Low")
)
filt.phlGsea.dotplot <- dotplot(filt.phlGseaOb, x = "GeneRatio", color = "p.adjust", 
        showCategory = 20, 
        label_format = 50, 
        split = ".sign") +
      facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Filtered Proximal High vs Low")

filt.phlCNET <- cnetplot(filt.phlGseaOb, foldChange = filt.phlGenes, colorEdge = TRUE,
                    node_label = "gene",
                    showCategory = "GOBP_NCRNA_METABOLIC_PROCESS")



filt.pmlGsea <- GSEA(filt.pmlGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
filt.pmlGseaOb <- setReadable(filt.pmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(filt.pmlGseaOb, "Filtered_Proximal_Med_vs_Low_GSEA_results.rds")
filt.pmlGsea.res <- as.data.frame(setReadable(filt.pmlGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(filt.pmlGsea.res, file = "Filtered_Proximal_Med_vs_Low_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in Low")
)
filt.pmlGsea.dotplot <- dotplot(filt.pmlGseaOb, x = "GeneRatio", color = "p.adjust", 
                                showCategory = 20, 
                                label_format = 50, 
                                split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Filtered Proximal Med vs Low")

filt.pmhGsea <- GSEA(filt.pmhGenes, minGSSize = 10, maxGSSize = 500, eps = 1e-30, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE, by = "fgsea", TERM2GENE = pathways)
filt.pmhGseaOb <- setReadable(filt.pmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_rds(filt.pmhGseaOb, "Filtered_Proximal_Med_vs_High_GSEA_results.rds")
filt.pmhGsea.res <- as.data.frame(setReadable(filt.pmhGsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
write.csv(filt.pmhGsea.res, file = "Filtered_Proximal_Med_vs_High_GSEA_results.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Enriched in Med", suppressed = "Enriched in High")
)
filt.pmhGsea.dotplot <- dotplot(filt.pmhGseaOb, x = "GeneRatio", color = "p.adjust", 
                                showCategory = 20, 
                                label_format = 50, 
                                split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  ggtitle("Filtered Proximal Med vs High")


filt.overlap <- intersect(filt.phlGsea.res$Description, filt.pmlGsea.res$Description)
filt.overlap <- intersect(filt.overlap, filt.pmhGsea.res$Description)

#-----Create a new dataframe containing just these terms and the enrichment scores
filt.phlOverlap <- filt.phlGsea.res[filt.phlGsea.res$Description %in% filt.overlap,]
filt.pmlOverlap <- filt.pmlGsea.res[filt.pmlGsea.res$Description %in% filt.overlap,]
filt.pmhOverlap <- filt.pmhGsea.res[filt.pmhGsea.res$Description %in% filt.overlap,]

#-----Collate information
filt.sharedGeneSets <- data.frame(matrix(nrow = length(filt.overlap)))
rownames(filt.sharedGeneSets) <- filt.overlap
filt.sharedGeneSets$Prox.HighVLow.ES <- filt.phlOverlap$enrichmentScore
filt.sharedGeneSets$Prox.MedVLow.ES <- filt.pmlOverlap$enrichmentScore
filt.sharedGeneSets$Prox.MedVHigh.ES <- filt.pmhOverlap$enrichmentScore
filt.sharedGeneSets$matrix.nrow...length.filt.overlap.. <- NULL
filt.sharedGeneSets$Prox.HighVLow.Padj <- filt.phlOverlap$p.adjust
filt.sharedGeneSets$Prox.HighVLow.qval <- filt.phlOverlap$qvalue
filt.sharedGeneSets$Prox.MedVLow.Padj <- filt.pmlOverlap$p.adjust
filt.sharedGeneSets$Prox.MedVLow.qval <- filt.pmlOverlap$qvalue
filt.sharedGeneSets$Prox.MedVHigh.Padj <- filt.pmhOverlap$p.adjust
filt.sharedGeneSets$Prox.MedVHigh.qval <- filt.pmhOverlap$qvalue

#-----Add a sign variable
filt.sharedGeneSets$Prox.HighVLow.Sign <- ifelse(filt.sharedGeneSets$Prox.HighVLow.ES > 0, "Enriched in High", "Enriched in Low")
filt.sharedGeneSets$Prox.MedVLow.Sign <- ifelse(filt.sharedGeneSets$Prox.MedVLow.ES > 0, "Enriched in Med", "Enriched in Low")
filt.sharedGeneSets$Prox.MedVHigh.Sign <- ifelse(filt.sharedGeneSets$Prox.MedVHigh.ES > 0, "Enriched in Med", "Enriched in High")

#-----Write as a csv
write.csv(filt.sharedGeneSets, file = "Filtered.SharedGeneSets_Between_3ComparisonGroups.csv")



#Read the files
raw <- read.csv("/Users/mikemartinez/Desktop/RNA_Sequencing/Counts/Walnut_Distal_Counts.csv", header = TRUE, sep = ",")
dds <- readRDS("/Users/mikemartinez/Desktop/RNA_Sequencing/October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/Distal_NoMed_dds.rds")

plot_total_counts(dds)
plot_library_complexity(dds)
plot_gene_detection(dds)
plot_biotypes(dds)

ddsdf <- as.data.frame(results(dds))
dds <- filter_genes(dds, min_count = 10, min_rep = 5)
ddsdf <- as.data.frame(results(dds))
vsd <- vst(dds)
mean_sd_plot(vsd)

plot_sample_clustering(vsd, anno_vars = "URoAClass", distance = "spearman")
pca_res <- plot_pca(vsd, show_plot = FALSE)
plot_loadings(pca_res, color_by = "gc_content")
plot_pca_scatters(vsd, n_PCs = 5, color_by = "URoAClass")
ggsave("PCA_scatters_Distal_High_vs_Low.pdf")
ds_res <- lfcShrink(dds, coef = "URoAClass_High_vs_Low", lfchreshold = log2(0.5), type = "normal", parallel = TRUE)
plot_ma(ds_res, dds)












