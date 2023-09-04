#####################################################
# Normalization and Quality Control #################
#####################################################

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)


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
all(colnames(data_raw) == rownames(data_annot))

dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = data_annot, design = ~ subgroup.a)

# Normalize the counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file="../TB_epithelium/normalized_counts.tsv", sep="\t", quote=F, col.names=NA)

####### PCA #######
# rlog transform the counts for PCA
rld <- rlog(dds, blind=TRUE)

# Plot the PCA
pdf("../TB_epithelium/Figures/pca_region.pdf")
plotPCA(rld, intgroup="subgroup.a", ntop = 1000) + labs(color = "region")
dev.off()

pdf("../TB_epithelium/Figures/pca_slide.pdf")
plotPCA(rld, intgroup="Slide", ntop = 1000) + labs(color = "Slide")
dev.off()

####### Heatmap ########
# Calculate pairwise correlation between samples
rld_mat <- assay(rld) #extract rlog matrix
rld_cor <- cor(rld_mat)


# Plot the heatmap
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
group_colors <- list(region = c('#F8766D', '#00ba38','#83b0fc'))
names(group_colors$region) <- c('TMP', 'TSM', 'TSS')

pdf("../TB_epithelium/Figures/counts_heatmap.pdf", width = 8)
pheatmap(rld_cor, 
  annotation_col = data.frame(region=data_annot$subgroup.a, row.names = rownames(data_annot)),
  annotation_colors = group_colors,
  fontsize = 9,
  fontsize_row = 6,
  show_colnames = F,
  show_rownames = F
  )
dev.off()

###### Sequencing stats ######
# Plot raw reads + mapped reads

reads_df <- data_annot[, c(1, 5, 12, 19, 21)]
reads_df[, 1] <- as.character(reads_df[, 1])
colnames(reads_df)[5] <- "percent_aligned"

pdf("../TB_epithelium/Figures/read_counts.pdf")
ggplot(reads_df, aes(x=...1, y=RawReads, fill=percent_aligned)) +
  geom_col() + labs( x = "Sample", y = "# of reads sequences", fill = "% aligned") +
  facet_wrap( ~ subgroup.a, strip.position = "bottom", scales = "free_x")
dev.off()

####### GSVA #######
library(org.Hs.eg.db)
library(GSVA)

# Retrieve the gene names in rowname() and make a column named "SYMBOL"
data_gsva <- rld_mat %>% as.data.frame() %>% tibble::rownames_to_column(var = "SYMBOL") 

# Map SYMBOL name to the Entrez IDs
symbol2enz <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=data_gsva$SYMBOL, 
                                    columns=("ENTREZID"), # desired format
                                    keytype="SYMBOL")  

data_gsva <- dplyr::inner_join(data_gsva, symbol2enz, by = "SYMBOL")
data_gsva <- data_gsva[!is.na(data_gsva$ENTREZID),] # remove NA

# Remove duplicates (Keeping the most highly expressed one)
sum(duplicated(data_gsva$ENTREZID))
gene_means <- rowMeans(data_gsva %>% dplyr::select(-SYMBOL, -ENTREZID))
data_gsva <- data_gsva %>%
  dplyr::mutate(gene_means) %>% # add as a column
  dplyr::select(SYMBOL, ENTREZID, gene_means, dplyr::everything())

filtered_data_gsva <- data_gsva %>% 
  dplyr::arrange(dplyr::desc(gene_means)) %>% # sort by mean expression (descending)
  dplyr::distinct(ENTREZID, .keep_all = TRUE) # only keeps the 1st duplicate it encounters

sum(duplicated(filtered_data_gsva$ENTREZID))

# Build the input matrix
filtered_data_gsva_mat <- filtered_data_gsva %>%
  dplyr::select(-SYMBOL, -gene_means) %>%
  tibble::column_to_rownames("ENTREZID") %>%
  as.matrix()

# Fetch Hallmark gene sets from msigdb
hallmark_gs <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H" # Only hallmark gene sets
)

# Make the gene set list
hallmarks_list <- split(
  hallmark_gs$entrez_gene, # The genes we want split into pathways
  hallmark_gs$gs_name # The pathways made as the higher levels of the list
)

# Run GSVA
gsva_results <- gsva(
  filtered_data_gsva_mat,
  hallmarks_list,
  method = "gsva",
  kcdf = "Gaussian", # because rlog does not yield integer counts
  # Minimum gene set size
  min.sz = 15,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

# Save results
gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  write.table(file="../TB_epithelium/gsva.tsv", sep="\t", quote=F, col.names=NA)

# Heatmap of GSVA results
pathway_heatmap <- pheatmap(gsva_results,
  annotation_col = data.frame(region=data_annot$subgroup.a, row.names = rownames(data_annot)), 
  show_colnames = FALSE, # Don't show sample labels
  fontsize_row = 6,
  fontsize = 6,
  legend_labels = c("-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "GSVA score"),
  legend_breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, max(gsva_results))
)

pdf("../TB_epithelium/Figures/gsva_heatmap.pdf", width = 10)
pathway_heatmap
dev.off()

# For paper/presentation version of heatmap: limit gene sets to those of interest









