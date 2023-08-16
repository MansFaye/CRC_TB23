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
data_raw <- data_raw %>% select(-"ROI 03_22 Bern A 15/03/23 | TMA-2203-A-X4-Y2-415-TB | T")
data_annot <- data_annot[-4,-33] # Drop CMS column too

# Make all count values integers
data_raw <- sapply(data_raw, as.integer)

# Create DESeq2 object
all(colnames(data_raw) %in% rownames(data_annot))
all(colnames(data_raw) == rownames(data_annot))

dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = data_annot, design = ~ subgroup.a)

# Normalize the counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, file="../TB_epithelium/normalized_counts.tsv", sep="\t", quote=F, col.names=NA)

####### PCA #######
# rlog transform the normalized counts for PCA
rld <- rlog(dds, blind=TRUE)

# Plot the PCA
pdf("../TB_epithelium/pca_region.pdf")
plotPCA(rld, intgroup="subgroup.a", ntop = 1000)
dev.off()

pdf("../TB_epithelium/pca_slide.pdf")
plotPCA(rld, intgroup="SegmentLabel", ntop = 1000)
dev.off()

####### Heatmap ########
# Calculate pairwise correlation between samples
rld_mat <- assay(rld) #extract rlog matrix
rld_cor <- cor(rld_mat)

# Plot the heatmap
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation_row = data.frame(region=data_annot$subgroup.a, row.names = rownames(data_annot)))

###### Sequencing stats ######
# Plot raw reads + mapped reads

reads_df <- data_annot[, c(1, 5, 12, 19, 21)]
reads_df[, 1] <- as.character(reads_df[, 1])
colnames(reads_df)[5] <- "percent_aligned"

pdf("../TB_epithelium/read_counts.pdf")
ggplot(reads_df, aes(x=...1, y=RawReads, fill=percent_aligned)) +
  geom_col() + labs(title = "Sequenced reads", x = "Sample", y = "# of reads", fill = "% aligned") +
  facet_wrap( ~ subgroup.a, strip.position = "bottom", scales = "free_x")
dev.off()

####### GSVA #######
#





