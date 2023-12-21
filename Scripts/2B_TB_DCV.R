######################################################
###### Spatial Deconvolution of Tumor Bud ROIs #######
######################################################

library(tidyverse)
library(DESeq2)
library(SpatialDecon)
library(GeomxTools)
library(Seurat)


# Load counts, profile matrix and metadata
data_annot <- read_csv(file = "../TB_epithelium/data_annot_tb_e.csv")
data_raw <- read_csv(file = "../TB_epithelium/data_raw_tb_e.csv", col_select = c(2:29))
prof_mtx <- read_csv(file = "../Master_files/custom_crc_top_profileMatrix.csv")

data_annot$subgroup.a <- as.factor(data_annot$subgroup.a)
data_annot <- data_annot %>% column_to_rownames(var="SegmentDisplayName")
data_raw <- data_raw %>% column_to_rownames(var="TargetName")
prof_mtx <- prof_mtx %>% column_to_rownames(var = "...1")

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

# Calculate the expected background
bg = derive_GeoMx_background(norm = normalized_counts,
                             probepool = rep(1, nrow(normalized_counts)),
                             negnames = "NegProbe-WTX")

# Run spatial decon
res <- spatialdecon(
  normalized_counts,
  bg,
  X = as.matrix(prof_mtx),
  raw = as.matrix(data_raw),
  align_genes = TRUE
)

# Visualization
## reshape proportion data to long format
prop_df <- t(res$beta)
prop_df <- data.frame(prop_df) %>% rownames_to_column(var = 'Sample') %>% 
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion")

group_df <- data_annot[, c(1, 5)] %>% rownames_to_column(var = 'Sample')
group_df <- group_df[, c(1, 3)]

prop_df <- prop_df %>% left_join(group_df, by = "Sample")

## order the samples by tissue region 
prop_df <- prop_df %>% arrange(subgroup.a, Sample) %>% mutate(Sample = factor(Sample, levels = unique(Sample)))

## plot barplots
## by region
pdf("../TB_epithelium/dcv_barplot_region.pdf", width = 8)
ggplot(prop_df, aes(x = subgroup.a, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = 'fill') +
  labs(title = "Cell Type Proportions by Region", x = 'Region', y = "Proportion") +
  theme_minimal()
dev.off()

# by sample
ggplot(prop_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell Type abundance scores by Sample", x = "", y = "Abundance score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



















