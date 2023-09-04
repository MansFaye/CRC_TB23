#############################################
####### TME Deconvolution around TBs ########
#############################################

library(tidyverse)
library(SpatialDecon)
library(GeomxTools)
library(DESeq2)
library(pheatmap)
library(viridis)
library(ggpubr)
library(car)
library(ggsignif)

# Load counts and metadata

data_annot <- read_csv(file = "../TB_TME/data_annot_tb_tme.csv")
data_raw <- read_csv(file = "../TB_TME/data_raw_tb_tme.csv", col_select = 2:39)
prof_mtx <- read_csv(file = "../Master_files/custom_crc_mid_profileMatrix.csv")


data_annot$subgroup.a <- as.factor(data_annot$subgroup.a)
data_annot <- data_annot %>% column_to_rownames(var="SegmentDisplayName")
data_raw <- data_raw %>% column_to_rownames(var="TargetName")
prof_mtx <- prof_mtx %>% column_to_rownames(var = "...1")

# Exclude TMA X4-Y2 (mislabeled)
data_raw <- data_raw %>% dplyr::select(-"ROI 03_22 Bern A 15/03/23 | TMA-2203-A-X4-Y2-415-TB | TME")
data_annot <- data_annot[-5,-32] # Drop CMS column too

# Make all count values integers
data_raw[] <- sapply(data_raw, as.integer)

# Create DESeq2 object
all(colnames(data_raw) %in% rownames(data_annot))
all(colnames(data_raw) == rownames(data_annot))

dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = data_annot, design = ~ subgroup.a)

# Normalize the counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# rlog transform the counts for PCA
rld <- rlog(dds, blind=TRUE)

# Plot the PCA
pdf("../TB_TME/pca_region.pdf")
plotPCA(rld, intgroup="subgroup.a", ntop = 1000) + labs(color = "region")
dev.off()

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

# Visualize results
## Heatmap
pheatmap(t(res$beta), 
         annotation_row = data.frame(region=data_annot$subgroup.a, row.names = rownames(data_annot)),
         col = inferno(10),
         show_rownames = F,
         )
## reshape proportion data to long format
prop_df <- t(res$prop_of_nontumor)
prop_df <- data.frame(prop_df) %>% rownames_to_column(var = 'Sample') %>% 
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion")

group_df <- data_annot[, c(1, 5)] %>% rownames_to_column(var = 'Sample')
group_df <- group_df[, c(1, 3)]

prop_df <- prop_df %>% left_join(group_df, by = "Sample")

## order the samples by tissue region 
prop_df <- prop_df %>% arrange(subgroup.a, Sample) %>% mutate(Sample = factor(Sample, levels = unique(Sample)))

## plot barplots
pdf("../TB_TME/dcv_barplot_region.pdf", width = 8)
ggplot(prop_df, aes(x = subgroup.a, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = 'fill') +
  labs(title = "Cell Type Proportions by Region", x = 'Region', y = "Proportion") +
  theme_minimal()
dev.off()

ggplot(prop_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~subgroup.a, nrow = 1, strip.position = "bottom") +
  labs(title = "Cell Type Proportions by Sample", x = "", y = "Proportion") +
  theme(axis.text.x = element_blank())

# Tests and boxplots of proportion for all cell types
comparisons <- list(c('TSM', 'TMP'), c('TSM', 'TSS'), c('TSS', 'TMP'))
v_celltypes <- unique(prop_df$CellType)

PlotCellType <- function(celltype){
  prop_df_cell <- prop_df[prop_df$CellType == celltype, 2:4]
  
  cat(paste0('Cell type: ', celltype), "\n")
  
  # Perform and print statistical tests
  kruskal_result <- compare_means(Proportion ~ subgroup.a, data = prop_df_cell, method = 'kruskal.test')
  wilcox_result <- compare_means(Proportion ~ subgroup.a, data = prop_df_cell, method = 'wilcox.test')
  
  print(kruskal_result)
  print(wilcox_result)
  
  try(
    {
    shapiro_result_tsm <- shapiro.test(prop_df_cell$Proportion[prop_df_cell$subgroup.a == 'TSM'])
    shapiro_result_tmp <- shapiro.test(prop_df_cell$Proportion[prop_df_cell$subgroup.a == 'TMP'])
    shapiro_result_tss <- shapiro.test(prop_df_cell$Proportion[prop_df_cell$subgroup.a == 'TSS'])
  
    print(shapiro_result_tsm)
    print(shapiro_result_tmp)
    print(shapiro_result_tss)
    }
  )
  # Create and save plots
  p_nonp <- ggplot(prop_df_cell, aes(x=subgroup.a, y=Proportion, fill = subgroup.a)) +
    geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
    labs(title = paste0('Proportion of ', celltype), x = 'Region', fill = 'Region')
  
  pdf(paste0("../TB_TME/", celltype,"_nonp.pdf"), width = 8)  
  try(print(p_nonp + stat_compare_means(method = "kruskal.test", label.x.npc = 'left', label.y = 0) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons#, label = 'p.signif'
    ))
  )
  dev.off()
  
  p <- ggplot(prop_df_cell, aes(x=subgroup.a, y=Proportion, fill = subgroup.a)) +
    geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
    labs(title = paste0('Proportion of ', celltype), x = 'Region', fill = 'Region')
  
  pdf(paste0("../TB_TME/", celltype,"_param.pdf"), width = 8)
  try(print(p + stat_compare_means(method = "anova", label.x.npc = 'left', label.y = 0) +
    stat_compare_means(method = 't.test', comparisons = comparisons, label = 'p.adj'
    ))
  )
  dev.off()
}

for(i in 1:length(v_celltypes)){
  PlotCellType(v_celltypes[i])
}




