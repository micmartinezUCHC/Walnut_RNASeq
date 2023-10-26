#Script to re-run the walnut RNASeq but this time, use a grouping variable to test differences between location and UroA
#Filtering out for IGVs, ribosomal genes, low counts

#-----Load libraries
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
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(patchwork))


#GSEA analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(msigdbr))

#-----Set working directory
setwd("/users/mikemartinez/Desktop/RNA_Sequencing/October_25th_GroupingVariable_tests/")

#-----Load in the counts
distal <- read.csv("../Counts/Walnut_Distal_Counts.csv", header = TRUE, sep = ",")
proximal <- read.csv("../Counts/Walnut_Proximal_Counts.csv", header = TRUE, sep = ",")
counts <- merge(distal, proximal, by = "X", all = TRUE)

#-----Get ENSEMBL IDs
counts$Ensembl <- mapIds(org.Hs.eg.db, key = counts$X, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first")
counts <- counts[!duplicated(counts$Ensembl),]
counts <- counts[!grepl("^IG", counts$X),]
counts <- counts[!grepl("^RPL", counts$X),]
counts <- counts[!grepl("^RPS", counts$X),]
rownames(counts) <- paste(counts$Ensembl, counts$X, sep = " - ")
counts$X <- NULL
counts$Ensembl <- NULL

#-----Read in the design file and filter out medium samples from the counts and design file
design <- read.csv("../metadata/New_Walnut_Design.csv", header = TRUE, sep = ",")
mediumSamples <- design[design$URoAClass == "Med",]$Sample
counts <- counts[,!(colnames(counts) %in% mediumSamples)]
design <- design[!design$Sample %in% mediumSamples,]

#-----Order the design file by Location, get the sample order, and arrange counts to match this
design <- design %>%
  arrange(Location)
rownames(design) <- design$Sample
sampleOrder <- rownames(design)
counts <- counts[,sampleOrder]

#-----Create a new grouping variable
design$Group <- paste(design$Location, design$URoAClass, sep = "_")
write.csv(design, file = "Combined_ProximalDistal_WalnutSamples_NoMedium_Producers_WithGroupingVar.csv")

#-----Sanity Check: Verify all samples exist between counts and design, and they are in the same order
all(colnames(counts) %in% rownames(design))
all(colnames(counts) == rownames(design))

#-----Filter out low counts
#Filter
minimumCountpergene <- 10
MinSampleWithminimumgeneCounts <- 6

#Filter out low read counts 
counts <- counts[rowSums(data.frame(counts>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]
write.csv(counts, file = "Combined_ProximalDistal_WalnutSamples_NoMedium_Producers.csv")

#-----Set up DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = design,
                                  design = ~ Group)
#-----Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)
save(dds, file = "Grouping_Variable_DDS_results.RDS")


#-----Plot PCA
vsd <- vst(dds)
PCA <- plotPCA(vsd, intgroup = "Group") +
  geom_text_repel(aes(label = rownames(design)), size = 3, max.overlaps = Inf) +
  ggtitle("Walnut by Grouping Variable, No Medium Producers") +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  theme_bw() +
  theme(legend.position = "bottom")
PCA

#-----Let's look at the results for proximal low vs distal low (keep in mine, the reference is the last argument in `contrast`)
lowGroups <- as.data.frame(results(dds, contrast = c("Group", "Distal_Low", "Proximal_Low")))
lowGroups <- lowGroups[order(lowGroups$padj, decreasing = FALSE),]
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(lowGroups))
lowGroups$Symbols <- labels
lowGroupsFilt <- lowGroups[lowGroups$padj < 0.05 & abs(highGroups$log2FoldChange) > 0.5,]
lowGroupsFilt <- lowGroupsFilt[order(lowGroupsFilt$log2FoldChange, decreasing = TRUE),]
lowGenes <- lowGroupsFilt$Symbols

#-----Let's look at the results for proximal high vs proximal low
highGroups <- as.data.frame(results(dds, contrast = c("Group", "Distal_High", "Proximal_High")))
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(highGroups))
highGroups$Symbols <- labels
highGroupsFilt <- highGroups[highGroups$padj < 0.05 & abs(highGroups$log2FoldChange) > 0.5,]
highGroupdFilt <- highGroupsFilt[order(highGroupsFilt$log2FoldChange, decreasing = TRUE),]
highGenes <- highGroupsFilt$Symbols

#-----Find the common genes between the low groups and high groups
common <- intersect(lowGenes, highGenes)

#-----Isolate the counts
counts <- as.data.frame(counts(dds, normalized = TRUE))
labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(counts))
counts$Symbols <- labels
countsCommon <- counts[counts$Symbols %in% common,]
rownames(countsCommon) <- countsCommon$Symbols
countsCommon$Symbols <- NULL

design.ordered <- design %>%
  arrange(Group)
newOrder <- design.ordered$Sample

countsCommon <- countsCommon[,newOrder]

Group <- design$Group
names(Group) <- design$Sample
Group <- as.data.frame(Group)

Group$Group <- factor(Group$Group, levels = c("Distal_Low", "Proximal_Low", "Distal_High", "Proximal_High"))


pheatmap(as.matrix(countsCommon),
         scale = "row",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         annotation_col = Group,
         fontsize_row = 3,
         main = "Distal-Proximal Low, vs Distal-Proximal High")

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

#-----Get the entrez gene IDs for the filered low group and high group
lowCP <- getEntrez(lowGroupsFilt)
highCP <- getEntrez(highGroupsFilt)



#-----Save the entrez IDs as a named list
GeneLists <- list(Low = names(lowCP),
                  High = names(highCP))
FC <- rbind(lowCP, highCP)

#-----Run compareClusters
clusters <- compareCluster(geneCluster = GeneLists, fun = enrichGO , OrgDb = org.Hs.eg.db)
cr <- setReadable(clusters, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save(cr, file = "CompareCluster_HighLow_BetweenLocations_DistalVsProximal.rds")
crdf <- as.data.frame(cr)
write.csv(crdf, file = "OverrepresentationAnalysis_results.csv")

ORAclusters <- dotplot(cr, showCategory = 30, label_format = 50)
ORAclusters <- ORAclusters + ggtitle("Over Representation Analysis: Distal vs Proximal (Low) and Distal vs Proximal (High)")
ggsave("ORA_clusters_between_High_and_Low.pdf", ORAclusters, width = 14, height = 8)

ORAgoCNET <- cnetplot(cr, foldChange = FC, node_label = "all", cex_label_category = 1, cex_label_gene = 1.5,
                      showCategory = 10, circular = FALSE)

organicAcidBinding <- c(" CYP2W1", " GLRA2", " GLUL", " GLUD1", " UGT1A8", " PCCA", " PHYH", " OTC", " FABP3", " TAT",
                        " NR1H4", " GRIN2D", " CYP26B1")
organicAcidResults <- highGroups[highGroups$Symbols %in% organicAcidBinding,]

EnhancedVolcano(organicAcidResults,
                lab = organicAcidResults$Symbols,
                title = "Organic Acid Binding Gene Intersection",
                subtitle = "",
                pCutoff = 0.1,
                FCcutoff = 0.1,
                legendPosition = "bottom",
                x = 'log2FoldChange',
                y = 'pvalue')

#-----Let's try with KEGG now 
clustersKegg <- compareCluster(geneCluster = GeneLists, fun = enrichKEGG)
ck <- setReadable(clustersKegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save(ck, file = "ORA_KEGG_clusters_between_High_and_Low.rds")
ckdf <- as.data.frame(ck)
write.csv(ckdf, file = "ORA_KEGG_between_High_and_Low.csv")

ORAKEGG <- dotplot(ck, showCategory = 17) +
  ggtitle("ORA: KEGG")
ggsave("ORA_KEGG.pdf", ORAKEGG, width = 12, height = 8)

#-----Let's try the same thing but with GSEA now 
lowgsea <- getEntrez(lowGroups)
highgsea <- getEntrez(highGroups)

#-----Save the entrez IDs as a named list
gsea_GeneLists <- list(Low = lowgsea,
                  High = highgsea)

clustersGSEA <- compareCluster(geneCluster = gsea_GeneLists, fun = gseGO , OrgDb = org.Hs.eg.db)
crgsea <- setReadable(clustersGSEA, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save(crgsea, file = "gsea_compareClusters_BetweenLocation_DistalVsProximal.rds")
crgseadf <- as.data.frame(crgsea)
write.csv(crgseadf, file = "GSEA_clusters_between_High_and_Low.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Activated", suppressed = "Suppressed")
)


gseaCluster <- dotplot(crgsea, showCategory = 30, label_format = 90)
gseaCluster <- gseaCluster + ggtitle("GSEA Distal va Proximal (Low) and Distal vs Proximal (High)") + 
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") 

ggsave("GSEA_clusters_betweenLocations_HighLow.pdf", gseaCluster, width = 14, height = 8)

#-----Let's do GSEA with KEGG now 
clustersGSEAkegg <- compareCluster(geneCluster = gsea_GeneLists, fun = gseKEGG) 
crgseakegg <- setReadable(clustersGSEAkegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save(crgseakegg, file = "gseaKEGG_compareClusters_BetweenLocation_DistalVsProximal.rds")
crgseakeggdf <- as.data.frame(crgseakegg)
write.csv(crgseadf, file = "GSEA_clusters_between_High_and_Low.csv")

custom_labels <- labeller(
  .default = label_value,
  .sign = c(activated = "Activated", suppressed = "Suppressed")
)


gseaCluster <- dotplot(crgsea, showCategory = 30, label_format = 90)
gseaCluster <- gseaCluster + ggtitle("GSEA Distal va Proximal (Low) and Distal vs Proximal (High)") + 
  facet_wrap(~ .sign, labeller = custom_labels, scales = "free_y") 

ggsave("GSEA_clusters_betweenLocations_HighLow.pdf", gseaCluster, width = 14, height = 8)


commonGenes <- rownames(countsCommon)
lowResultsCommon <- lowGroups[lowGroups$Symbols %in% commonGenes,]
highResultsCommon <- highGroups[highGroups$Symbols %in% commonGenes,]


EnhancedVolcano(lowGroups,
                lab = lowGroups$Symbols,
                title = "Distal Low vs Proximal Low",
                subtitle = "",
                legendPosition = "bottom",
                x = 'log2FoldChange',
                y = 'pvalue') 

EnhancedVolcano(highGroups,
                lab = highGroups$Symbols,
                title = "Distal High vs Proximal High",
                subtitle = "",
                legendPosition = "bottom",
                x = 'log2FoldChange',
                y = 'pvalue') 



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
library(GSVA)

m_t2gC2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene, gs_subcat) %>%
  filter(gs_subcat %in% c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"))

m_t2gC5 <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene,  gs_subcat) %>%
  filter(gs_subcat %in% c("GO", "GO:BP", "GO:CC", "GO:MF"))

#-----Combine into one large dataframe of gene sets
pathways <- rbind(m_t2gC2, m_t2gC5)




#Set seed
set.seed(03061999)

#Let's get the results of our gene expression data
results <- as.data.frame(results(dds, contrast = c("Group", "Distal_Low", "Proximal_Low")))
results$gene <- rownames(results)

#We need the counts as an input with rownames as entrez IDs
counts <- as.data.frame(counts(dds, normalized = TRUE))
counts$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(counts))
counts$Entrez <- mapIds(org.Hs.eg.db, key = counts$Ensembl,
                         column = "ENTREZID", keytype = "ENSEMBL",
                         multiVals = "first")

#We need to remove duplicate entrez IDs so we can set them as rownames
entrezCol <- counts$Entrez
unique <- duplicated(entrezCol)
counts <- counts[!unique,]
rownames(counts) <- counts$Entrez
counts$Entrez <- NULL
counts$Ensembl <- NULL

#Now we need to get the geneset we want and convert each category into a list
m_t2gC5 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  filter(gs_subcat %in% c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"))
GO_list <- split(m_t2gC5$entrez_gene,
                 m_t2gC5$gs_name)

#Now we run GSVA. Note: you can use Gaussian kcdf if you have tpms or a vsd transformed object, but since we have integer counts, we use Poisson                                            
gsva <- gsva(as.matrix(counts), 
             GO_list,
             kcdf = "Poisson",
             min.sz = 1,
             max.sz = 500)
save(gsva, file = "C2_GSVA_Object.rds")



gsva_results <- as.data.frame(gsva)
write.csv(gsva_results, file = "C2_gsva_results.csv")
gsva_results$Desc <- rownames(gsva_results)


search_term <- "Microbial"
matching_rows <- grep(search_term, gsva_results$Desc, ignore.case = TRUE)

gsva_resultsEstrogen <- gsva_results[matching_rows,]
gsva_resultsEstrogen$Desc <- NULL

gsva_resultsFattyAcids <- gsva_results[matching_rows,]
gsva_resultsFattyAcids$Desc <- NULL

gsva_resultsMicrobial <- gsva_results[matching_rows,]
gsva_resultsMicrobial$Desc <- NULL

Group <- design$URoAClass
names(Group) <- design$Sample
Group <- as.data.frame(Group)
Group$Location <- design$Location

Group <- Group %>%
  arrange(desc(Group))

sample_order <- rownames(Group)
gsva_resultsEstrogen <- gsva_resultsEstrogen[,sample_order]

Group$Group <- factor(Group$Group, levels = c("Distal_Low", "Proximal_Low", "Distal_High", "Proximal_High"))


pheatmap(gsva_resultsEstrogen,
         annotation_col = Group,
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = TRUE,
         main = "GSVA Estrogen Terms")

pheatmap(gsva_resultsFattyAcids[1:11,],
         annotation_col = Group,
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = TRUE,
         main = "GSVA Fatty Acid Terms")


pheatmap(gsva_resultsMicrobial[1,],
         annotation_col = Group,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row",
         show_rownames = TRUE,
         main = "GSVA Microbial Terms")

gsva_results$Desc <- NULL

collapsedGSVA <- as.matrix(gsva_results)
collapsedGSVA <- c(collapsedGSVA)
hist(as.numeric(collapsedGSVA))





