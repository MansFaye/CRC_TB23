# Welcome!
Code accompanying our project on "Identifcation of phenotypic switching of colorectal tumour buds in response to tumour-microenvironment using spatially resolved transcriptomics", abstract available [here](https://link.springer.com/article/10.1007/s00292-023-01249-7).

#### Overview
This repo contains code for all analysis steps, from gene counts to pathway enrichment analyses and spatial deconvolution. 
The data processing steps were handled by Nanostring's `GeoMx NGS pipeline v2.3.3.10`.

<img src="./Master_files/Workflow_Fig.png" alt="drawing" width="650"/>

#### Profile matrices
We built custom cell profile matrices from colorectal cancer single-cell data for the deconvolution tasks. The data is freely accessible from [Pelka & al](https://www.cell.com/cell/fulltext/S0092-8674(20)30870-9)'s article "Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front".
Two matrices are available: the "top" matrix groups together similar cell types, and the "mid" categories are a bit more detailed. Short descriptions for the cell type names are available [here](./Master_files/name_mapping.xlsx).

#### DGE, GSEA and deconvolution


