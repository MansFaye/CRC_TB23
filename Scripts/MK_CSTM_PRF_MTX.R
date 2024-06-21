#######################################################
###### Building a custom profile matrix for CRC #######
#######################################################
# As the single cell dataset from Pelka & al. is extremely large, we recommend running this on 
# an HPC cluster or a cloud instance. R struggles with managing memory when loading the entire count matrix 
  
# Load libraries
library(tidyverse)
library(SpatialDecon)
library(Seurat)


# Load the single-cell data and metadata
annot <- read_tsv(file = "./Master_files/crc10x_tSNE_cl_global.tsv")
annot <- annot[2:nrow(annot),] # remove 2nd headers

sc_data <- Read10X(data.dir = "./Master_files/sc_data/", gene.column = 2) # use SYMBOL names

ctrl <- CreateSeuratObject(counts = sc_data, min.features = 200, min.cells = 10) # Create the Seurat object
rm(sc_data)

# examine cell types in the dataset
table(annot$ClusterTop)
table(annot$ClusterMidway)

# Sift through dataset
annot <- annot %>% dplyr::filter(NAME %in% ctrl@assays$RNA@counts@Dimnames[[2]])

# Randomly sample 250 cells from each Midway cell type
sampled_names <- character(0)

unique_cell_types <- unique(annot$ClusterMidway)

set.seed(123)
for (cell_type in unique_cell_types) { # for 250 samples per cell type
  # Filter the dataframe to get rows with the current cell type
  cell_type_df <- annot[annot$ClusterMidway == cell_type, ]
  
  # Randomly sample n sample names
  sampled_samples <- sample(cell_type_df$NAME, size = 250, replace = FALSE)
  
  # Store the sampled sample names in the list
  sampled_names <- c(sampled_names, sampled_samples)
}

ctrl_250 <- subset(x = ctrl, cells = sampled_names)

rm(ctrl)

# remove gene programs
genes <- as.vector(ctrl_250@assays$RNA@counts@Dimnames[1])[[1]]
genes <- genes[!grepl("^p", genes)]

ctrl_250 <- ctrl_250[rownames(x = ctrl_250) %in% genes, ]

# Evaluate variation and select most variable genes
ctrl_250 <- FindVariableFeatures(ctrl_250, selection.method = "vst", nfeatures = 2000) # 2000 most variable genes

top10_250 <- head(VariableFeatures(ctrl_250), 10)

plot250 <- VariableFeaturePlot(ctrl_250)

LabelPoints(plot = plot250, points = top10_250)

# Extract the dataframe
counts_250 <- GetAssayData(object = ctrl_250, slot = "counts") # "counts" for raw data
counts_250 <- counts_250[which(rownames(counts_250) %in% VariableFeatures(object = ctrl_250)), ]


# Make profile matrix
annot <- annot[annot$NAME %in% sampled_names,]

custom_mtx_top <- create_profile_matrix(mtx = counts_250,            # cell x gene count matrix, raw
                         cellAnnots = annot,  # cell annotations with cell type and cell name as columns 
                         cellTypeCol = "ClusterTop",  # column containing cell type
                         cellNameCol = "NAME",           # column containing cell ID/name
                         matrixName = "custom_crc_top", # name of final profile matrix
                         outDir = "../Master_files/",                    # path to desired output directory, set to NULL if matrix should not be written
                         normalize = TRUE,                # Should data be normalized? 
                         discardCellTypes = TRUE)    

custom_mtx_midway <- create_profile_matrix(mtx = counts_250,           
                                        cellAnnots = annot,   
                                        cellTypeCol = "ClusterMidway",  
                                        cellNameCol = "NAME",
                                        matrixName = "custom_crc_mid", 
                                        outDir = "../Master_files/",                    
                                        normalize = TRUE,                 
                                        discardCellTypes = TRUE)    

