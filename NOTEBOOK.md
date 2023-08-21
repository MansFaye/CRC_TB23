## Section 1: Tumor Bud GE in TSM, TMP, TSS
### 15/8 - 17/8: Exploration/QC (PCA, GSVA)

We started by manually normalizing the raw reads in order to asses sample-level variation with PCA. We then performed a Gene-set variation analysis, allowing us to see if our sequencing results align with what we would expect.

The normalization method used was DESeq2's median of ratios method, which accounts for sequencing depth and RNA composition. This allows comparisons between samples, but not within them. The counts were log-transformed (using rlog) before PCA to improve the separation/clustering of samples.
The sample `TMA-2203-A-X4-Y2-415-TB` was excluded from the analysis due to a mislabeling problem

We generated PCA plots highlighting both the sample region (TMP, TSM, TSS) and the tissue slide of origin (A or B). A heatmap was also generated to visualize the correlation between samples.

We reviewed the main statistics from the sequencing and data preprocessing, plotting the # reads sequenced and proportion of mapped reads.

The gene-set variation analysis was conducted using the Molecular Signatures DB's hallmark gene sets. We made another heatmap to visualize the GSVA results. 

### 18/8 - 21/8: DGE, GSEA, ORA

We ran the DGE both with and without the batch effect in the design. When extracting results, both the independent filtering and cook's distance cutoff parameters were set to true.    

Without batch effect:   

At alpha cutoff = 0.05, 25 genes were significantly DE between the TSM and TSS samples. For the two pairs involving the TMP samples, 0 genes were sig. DE. 
Looking at both the PCA and GSVA results, the results for TMP are not too surprising (see overlap).

With batch effect:  

At alpha cutoff = 0.05, 8 genes were significantly DE between the TSM and TSS samples. The most significant ones are the same as for the first run. For the two pairs involving the TMP samples, 0 genes were sig. DE.

*Note: When comparing samples from TSM and TSS, the majority of genes are filtered out due to low mean count.This is not the case for the other two comparisons. When looking at the summary, we can see that the mean count cutoff is systematically higher for that pair of conditions than for the others.* 

The volcano plots were made using `EnhancedVolcano`. The axis show log2FC and *p-value* on a log base 10 scale, not the adjusted p-value. The color annotation, however, show which genes are significantly DE based on the adjusted p-value.

For the GSEA, 25 of the 6920 genes analyzed did not map to an ENTREZ ID and were filtered out. None of these were sig. DE or passed the count cutoff filtering for TSM and TSS.  
Log2FC was used to rank the genes, and the GSEA was conducted using the `GSEA` function from the `clusterProfiler` package. Similarily to the GSVA, the gene sets were fetched from the Molecular Signature DB using the package `AnnotationDbi`. Dot plots were generated for each comparison, showing the adjusted p-value, gene count and gene ratio for each significantly enriched pathway.


The ORA was only ran on the TSS/TSM pair, as it is the only one that had differentially expressed genes. We used the `enricher` function from the `clusterProfiler` package. Results can be viewed in an upset plot.

