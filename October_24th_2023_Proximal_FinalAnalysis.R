#This is to finalize the analysis on the proximal walnut RNA seq between high and low

#Load libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(ggh4x)
library(pheatmap)
library(EnhancedVolcano)

#Set working directory
setwd("/Users/mikemartinez/Desktop/RNA_Sequencing/October_24th_Final_Proximal_NoMed_Analysis/")

#Read in the proximal dds object with no medium samples included
dds <- readRDS("../October_20th_2023_Proximal_No_Medium/Proximal_dds_NoMedium.rds")
res <- as.data.frame(results(dds))
res$Gene <- rownames(res)
counts <- as.data.frame(counts(dds, normalized = TRUE))
counts$Gene <- rownames(counts)
results <- merge(res, counts, by = "Gene", all = TRUE)
counts$Gene <- NULL

write.csv(results, file = "Proximal_Colon_High_vs_Low_all_DEGs_NoMedium.csv")

results.filt <- results[results$padj < 0.05,]
results.fc <- results[abs(results$log2FoldChange) > 0.5,]
results.fc <- results.fc[order(results.fc$log2FoldChange, decreasing = TRUE),]

#-----Plot the 3 significant genes that came up
DNASE2counts <- plotCounts(dds, gene = "ENSG00000105612 - DNASE2", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "DNASE2",
                            xlab = "URoA Group", returnData = TRUE)
DNASE2plot <- ggplot(DNASE2counts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("DNASE2") +
  theme(plot.title = element_text(size = 26))
DNASE2plot
ggsave("Proximal_High_vs_Low_DNASE2plot.tiff", DNASE2plot, dpi = 800)

VIPcounts <- plotCounts(dds, gene = "ENSG00000146469 - VIP", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "VIP",
                           xlab = "URoA Group", returnData = TRUE)
VIPplot <- ggplot(VIPcounts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("VIP") +
  theme(plot.title = element_text(size = 26))
VIPplot
ggsave("Proximal_High_vs_Low_VIPplot.tiff", VIPplot, dpi = 800)

HLAcounts <- plotCounts(dds, gene = "ENSG00000204287 - HLA-DRA", intgroup = "URoAClass", normalized = TRUE, transform = FALSE, main = "HLA-DRA",
                        xlab = "URoA Group", returnData = TRUE)
HLAplot <- ggplot(HLAcounts, aes(x = URoAClass, y = count, fill = URoAClass)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.09) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  guides(fill = FALSE) +
  ggtitle("HLA-DRA") +
  theme(plot.title = element_text(size = 26))
HLAplot
ggsave("Proximal_High_vs_Low_HLA-DRAplot.tiff", HLAplot, dpi = 800)


design <- read.csv("/Users/mikemartinez/Desktop/RNA_Sequencing/Metadata/ProximalWalnut_Design.csv", header = TRUE, sep = ",")
design <- design[design$URoAClass != "Med",]
Group <- design$URoAClass
names(Group) <- design$Sample
Group <- as.data.frame(Group)

ProximalHM <- pheatmap(results.fc[,8:ncol(results.fc)],
                       scale = "row",
                       color = viridis(500),
                       cluster_cols = TRUE,
                       cluster_rows = TRUE,
                       show_rownames = FALSE,
                       annotation_col = Group,
                       main = "Proximal High vs Low, all Genes with abs(log2FC) > 0.5",
                       fontsize = 20)


#Let's do the same thing for distal really quick too 
distdds <- readRDS("../October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/Distal_NoMed_dds.rds")
distresults <- as.data.frame(results(distdds))
distresults$Gene <- rownames(distresults)
distCounts <- as.data.frame(counts(distdds))
distCounts$Gene <- rownames(distCounts)
distalResults <- merge(distresults, distCounts, by = "Gene", all = TRUE)
distCounts$Gene <- NULL
distResults.filt <- distalResults[distalResults$padj < 0.05,]
distResults.filt <- na.omit(distResults.filt)
distResults.fc <- distalResults[abs(distalResults$log2FoldChange) > 0.5,]
distResults.fc <- na.omit(distResults.fc)


distDesign <- read.csv("../Metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
distDesign <- distDesign[distDesign$Location == "Distal" & distDesign$URoAClass != "Med",]
Group <- distDesign$URoAClass
names(Group) <- distDesign$Sample
Group <- as.data.frame(Group)

DistalHM <- pheatmap(distResults.fc[,8:ncol(results.fc)],
                       scale = "row",
                       color = viridis(500),
                       cluster_cols = TRUE,
                       cluster_rows = TRUE,
                       show_rownames = FALSE,
                       annotation_col = Group,
                      clustering_method = "complete",
                       main = "Distal High vs Low, all Genes with abs(log2FC) > 0.5",
                       fontsize = 20)


#-----Quickly remake the GSEA plot for distal
distgsea <- readRDS("../October_20th_2023_ProxDist_merged/DisatlOnly_NoMedium/Dist_NoMedium_High_vs_Low_GSEA_results.rds")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Activated", suppressed = "Suppressed")
)
dhlGsea.dotplot <- dotplot(distgsea, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  theme(strip.text = element_text(size = 20),
        title = element_text(size = 24)) +
  ggtitle("Distal High vs Low")
ggsave("Distal_High_vs_Low_GSEA_relabelled_dotplot.tiff", dhlGsea.dotplot, dpi = 800, width = 18, height = 12)


#-----PCA for proximal results
vsd <- vst(dds)
rownames(design) <- design$Sample

PCA <- plotPCA(vsd, intgroup = "URoAClass") +
  geom_text_repel(aes(label = rownames(design)), size = 3, max.overlaps = Inf) +
  ggtitle("Proximal Colon No Medium Producers") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA
ggsave("Proximal_PCA.tiff", PCA, width = 12, height = 8, dpi = 800)


#-----GSEA for proximal
gsea <- readRDS("../October_20th_2023_Proximal_No_Medium/Proximal_NoMedium_High_vs_Low_GSEA_results.rds")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Activated", suppressed = "Suppressed")
)
Gsea.dotplot <- dotplot(gsea, x = "GeneRatio", color = "p.adjust", 
                           showCategory = 20, 
                           label_format = 50, 
                           split = ".sign") +
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") +
  theme(strip.text = element_text(size = 20),
        title = element_text(size = 24)) +
  ggtitle("Proximal High vs Low")
ggsave("Proximal_High_vs_Low_GSEA_relabelled_dotplot.tiff", Gsea.dotplot, dpi = 800, width = 18, height = 12)


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
rownames(results) <- results$Gene
genes <- getEntrez(results)



proximal_heatplot <-heatplot(gsea, showCategory = c("GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
                                "GOMF_HORMONE_RECEPTOR_BINDING",
                                "PID_IL12_STAT4_PATHWAY",
                                "GOMF_CALCIUM_ION_BINDING",
                                "GOBP_NEGATIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
                                "PID_AP1_PATHWAY"), 
         foldChange = genes)
ggsave("proximal_heatplot.pdf", proximal_heatplot, width = 18, height = 12, dpi = 800)

hormone <- heatplot(gsea, showCategory = c("GOMF_HORMONE_RECEPTOR_BINDING",
                                           "GOMF_CALCIUM_ION_BINDING"), foldChange = genes,
                    symbol = "rect")
hormone <- hormone + theme(axis.text.x = element_text(size = 10))
ggsave("proximal_hormoneReceptorBinding.pdf", hormone, width = 18, height = 12, dpi = 800)

cancer <- heatplot(gsea, showCategory = c("WP_GASTRIC_CANCER_NETWORK_2",
                                             "WP_INTEGRATED_CANCER_PATHWAY"), foldChange = genes, symbol = "rect")
cancer <- cancer + theme(axis.text.x = element_text(size = 10))
ggsave("proximal_cancerGeneSets.pdf", cancer, width = 18, height = 12, dpi = 800)


results.filtpadj <- results[order(results$padj, decreasing = FALSE),]
results.filtpadj$Symbol <- gsub("^(.*?) - .*", "\\1", rownames(results.filtpadj))
results.padj.filtered <- results.filtpadj[abs(results.filtpadj$log2FoldChange) > 0.5,]

topHits <- results.padj.filtered$Symbol[1:1000]
topFCHits <- results.fc$Symbol[1:100]

write.csv(topHits, file = "top100padjGenes.ProximalHighvsLow.csv")









