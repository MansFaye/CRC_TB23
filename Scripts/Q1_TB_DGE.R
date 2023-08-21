#############################################
## Differential Gene Expression, ORA, GSEA ##
#############################################

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

# Load counts and metadata
data_annot <- read_csv(file = "../TB_epithelium/data_annot_tb_e.csv")
data_raw <- read_csv(file = "../TB_epithelium/data_raw_tb_e.csv", col_select = c(2:29))

data_annot$subgroup.a <- as.factor(data_annot$subgroup.a)
data_annot <- data_annot %>% column_to_rownames(var="SegmentDisplayName")
data_raw <- data_raw %>% column_to_rownames(var="TargetName")

data_annot[which(data_annot$T == TRUE), 12] <- "A"
data_annot[which(data_annot$T == FALSE), 12] <- "B"
colnames(data_annot)[12] <- "Slide" 

# Exclude TMA X4-Y2 (mislabeled)
data_raw <- data_raw %>% dplyr::select(-"ROI 03_22 Bern A 15/03/23 | TMA-2203-A-X4-Y2-415-TB | T")
data_annot <- data_annot[-4,-33] # Drop CMS column too

# Make all count values integers
data_raw[] <- sapply(data_raw, as.integer)

# Create DESeq2 object
all(colnames(data_raw) %in% rownames(data_annot))
all(colnames(data_raw) == rownames(data_annot)) # check that row and colnames align

design1 <- ~ subgroup.a
design2 <- ~ Slide + subgroup.a # Batch effect

dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = data_annot, design = design1)

# Run the DGE
dds <- DESeq(dds)

# Extract results for each pairwise comparison
contrast_SMSS <- c("subgroup.a", "TSS", "TSM") # TSM is base

res_SMSS_unshrunken <- results(dds, contrast=contrast_SMSS, alpha = 0.05, independentFiltering=T, cooksCutoff = T) # filtering always on
res_SMSS <- lfcShrink(dds, contrast=contrast_SMSS, res=res_SMSS_unshrunken, type = "normal") # Shrink Log2 fold-change
res_SMSS_tbl <- res_SMSS %>% data.frame()

contrast_SMMP <- c("subgroup.a", "TMP", "TSM") # TSM is base

res_SMMP <- results(dds, contrast=contrast_SMMP, alpha = 0.05, independentFiltering=T, cooksCutoff = T)
res_SMMP <- lfcShrink(dds, contrast=contrast_SMMP, res=res_SMMP, type = "normal") 
res_SMMP_tbl <- res_SMMP %>% data.frame()

contrast_SSMP <- c("subgroup.a", "TMP", "TSS") # TSS is base

res_SSMP <- results(dds, contrast=contrast_SSMP, alpha = 0.05, independentFiltering=T, cooksCutoff = T)
res_SSMP <- lfcShrink(dds, contrast=contrast_SSMP, res=res_SSMP, type = "normal") 
res_SSMP_tbl <- res_SSMP %>% data.frame()

summary(res_SMSS)
summary(res_SMMP)
summary(res_SSMP)

# Extract sig DE genes in SMSS 
sig_SMSS <- res_SMSS_tbl %>% filter(padj < 0.05)

# MA Plot, visualize shrinkage
plotMA(res_SMSS_unshrunken, ylim=c(-3,3), main ="TSS vs TSM, unshrunken", colSig="red2")

plotMA(res_SMSS, ylim=c(-3,3), main ="TSS vs TSM, shrunken", colSig="red2")

plotMA(res_SMMP, ylim=c(-3,3), main ="TMP vs TSM, shrunken", colSig="red2")

plotMA(res_SSMP, ylim=c(-3,3), main ="TMP vs TSS, shrunken", colSig="red2")

# Volcano plots
pdf("../TB_epithelium/volcano_SMSS.pdf", width = 9)
EnhancedVolcano(res_SMSS_tbl,
                lab = rownames(res_SMSS_tbl),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DGE in TB: TSS vs TSM (Base)",
                subtitle = element_blank(),
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0.58,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 15,
                pointSize = 1,
                labSize = 4)
dev.off()

pdf("../TB_epithelium/volcano_SMMP.pdf", width = 9)
EnhancedVolcano(res_SMMP_tbl,
                lab = rownames(res_SMMP_tbl),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DGE in TB: TMP vs TSM (Base)",
                subtitle = element_blank(),
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0.58,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 15,
                pointSize = 1,
                labSize = 4)
dev.off()

pdf("../TB_epithelium/volcano_SSMP.pdf", width = 9)
EnhancedVolcano(res_SSMP_tbl,
                lab = rownames(res_SSMP_tbl),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DGE in TB: TMP vs TSS (Base)",
                subtitle = element_blank(),
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0.58,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 15,
                pointSize = 1,
                labSize = 4)
dev.off()

# Sig Genes Boxplots



###### DGE, with batch effect ######
dds2 <- DESeqDataSetFromMatrix(countData = data_raw, colData = data_annot, design = design2)
# Run the DGE
dds2 <- DESeq(dds2)

# Extract results for each pairwise comparison

res_SMSS2 <- results(dds2, contrast=contrast_SMSS, alpha = 0.05, independentFiltering=T, cooksCutoff = T) # filtering always on
res_SMSS2 <- lfcShrink(dds2, contrast=contrast_SMSS, res=res_SMSS2, type = "normal") # Shrink Log2 fold-change

res_SMMP2 <- results(dds2, contrast=contrast_SMMP, alpha = 0.05, independentFiltering=T, cooksCutoff = T)
res_SMMP2 <- lfcShrink(dds2, contrast=contrast_SMMP, res=res_SMMP2, type = "normal") 

res_SSMP2 <- results(dds2, contrast=contrast_SSMP, alpha = 0.05, independentFiltering=T, cooksCutoff = T)
res_SSMP2 <- lfcShrink(dds2, contrast=contrast_SSMP, res=res_SSMP2, type = "normal") 

summary(res_SMSS2)
summary(res_SMMP2)
summary(res_SSMP2)


###### GSEA ######
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Shift gene names to SYMBOL column
res_SMSS_gsea <- res_SMSS_tbl %>% rownames_to_column(var="SYMBOL")

res_SMMP_gsea <- res_SMMP_tbl %>% rownames_to_column(var="SYMBOL")

res_SSMP_gsea <-res_SSMP_tbl %>% rownames_to_column(var="SYMBOL")

# Map SYMBOL name to the Entrez IDs
symbol2enz <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res_SMSS_gsea$SYMBOL, 
                                    columns=("ENTREZID"), # desired format
                                    keytype="SYMBOL")  

res_SMSS_gsea <- dplyr::inner_join(res_SMSS_gsea, symbol2enz, by = "SYMBOL")

res_SMMP_gsea <- dplyr::inner_join(res_SMMP_gsea, symbol2enz, by = "SYMBOL")

res_SSMP_gsea <- dplyr::inner_join(res_SSMP_gsea, symbol2enz, by = "SYMBOL")

# Remove the unmapped genes, and only keep log2FC and ID columns
res_SMSS_gsea <- res_SMSS_gsea[!is.na(res_SMSS_gsea$ENTREZID), ] %>% dplyr::select(ENTREZID, log2FoldChange) 

res_SMMP_gsea <- res_SMMP_gsea[!is.na(res_SMMP_gsea$ENTREZID), ] %>% dplyr::select(ENTREZID, log2FoldChange)

res_SSMP_gsea <- res_SSMP_gsea[!is.na(res_SSMP_gsea$ENTREZID), ] %>% dplyr::select(ENTREZID, log2FoldChange)

# Check duplicates
sum(duplicated(res_SMSS_gsea$ENTREZID))


# Rank according to Log2FC
res_SMSS_gsea <- res_SMSS_gsea[order(-res_SMSS_gsea$log2FoldChange), ] 
log2fc_SMSS <- res_SMSS_gsea$log2FoldChange
names(log2fc_SMSS) <- res_SMSS_gsea[, 1]

res_SMMP_gsea <- res_SMMP_gsea[order(-res_SMMP_gsea$log2FoldChange), ] 
log2fc_SMMP <- res_SMMP_gsea$log2FoldChange
names(log2fc_SMMP) <- res_SMMP_gsea[, 1]

res_SSMP_gsea <- res_SSMP_gsea[order(-res_SSMP_gsea$log2FoldChange), ] 
log2fc_SSMP <- res_SSMP_gsea$log2FoldChange
names(log2fc_SSMP) <- res_SSMP_gsea[, 1]

# Fetch Hallmark gene sets from msigdb
hallmark_gs <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H") %>% dplyr::select(gs_name, entrez_gene)

# Run GSEA
gsea_SMSS <- GSEA(log2fc_SMSS, TERM2GENE = hallmark_gs)

gsea_SMMP <- GSEA(log2fc_SMMP, TERM2GENE = hallmark_gs)

gsea_SSMP <- GSEA(log2fc_SSMP, TERM2GENE = hallmark_gs)

# Visualize GSEA results
dot_SMSS <- dotplot(gsea_SMSS, showCategory=15, split=".sign", font.size=8, size = NULL, label_format = 90) + 
  ggtitle("Hallmark, TSS vs TSM (Base)") +
  facet_grid(~.sign) + 
  theme(plot.title.position = "plot")
# Note: GeneRatio = # of core (contributed to the enrichment score) DE genes in pathway / total count of genes in pathway

dot_SMMP <- dotplot(gsea_SMMP, showCategory=15, split=".sign", font.size=8, size = NULL, label_format = 90) + 
  ggtitle("Hallmark, TMP vs TSM (Base)") +
  facet_grid(~.sign) + 
  theme(plot.title.position = "plot")

dot_SSMP <- dotplot(gsea_SSMP, showCategory=15, split=".sign", font.size=8, size = NULL, label_format = 90) + 
  ggtitle("Hallmark, TMP vs TSS (Base)") +
  facet_grid(~.sign) + 
  theme(plot.title.position = "plot")

pdf("../TB_epithelium/gsea_SMSS_dot.pdf", width = 8)
dot_SMSS
dev.off()

pdf("../TB_epithelium/gsea_SMMP_dot.pdf", width = 8)
dot_SMMP
dev.off()

pdf("../TB_epithelium/gsea_SSMP_dot.pdf", width = 8)
dot_SSMP
dev.off()

# Extract SYMBOL names of genes in core enrichment


####### Overrepresentation Analysis (ORA) #######
#Get ENTREZ IDs of sig. DE genes
sig_SMSS <- sig_SMSS %>% rownames_to_column(var="SYMBOL")
symbol2enz <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=sig_SMSS$SYMBOL, 
                                    columns=("ENTREZID"), # desired format
                                    keytype="SYMBOL")
sig_SMSS <- dplyr::inner_join(sig_SMSS, symbol2enz, by = "SYMBOL")
gene <- sig_SMSS$ENTREZID

# Run ORA and make upset plot
ora_SMSS <- enricher(gene, TERM2GENE=hallmark_gs)

pdf("../TB_epithelium/ora_SMSS.pdf", width = 10)
upsetplot(ora_SMSS) + labs(title = "TSS vs TSM") + theme(plot.title.position = "plot")
dev.off()




