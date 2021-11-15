# BPDCN
Scripts used to generate results for BPDCN paper (Fall 2021)

## General

### 200116_FunctionsGeneral.R
A set of functions for computing signature scores in data, scaling, and plotting developed by Peter van Galen. This code is called within various scripts in this repository to simplify colorizing cells within a UMAP and TSNE plots by a selected gene signature.

### Seurat_to_h5ad.R
This short script converts a Seurat object to the h5ad format for use with the sceasy package. This was not used in the final results in the paper, but nonetheless may be useful for others.

### seurat-3.0.R
This script is a modification of the the standard clustering Seurat tutorial (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). It also identifies marker genes and colorizes UMAP plots based on the signature scores from the 200116_FunctionsGeneral.R script. The section at the bottom of the script reads in a Seurat.rds object as writes a gene by cell matrix for use with cellHarmony.

### Seurat_Integration3.R
This script takes a list of Seurat.rds for integration using Seurat's Integration function. This script also generates UMAP plots colored by Seurat cluster and signatures scores (see above) to assist with determining the cell type associated with each cluster. Note that cell cycle genes were removed from consideration when defining anchors in line 42. The bottom of the script contains legacy code used to color plots by specific genes.

### T_cell_subcluster.R


## Classification

### random_forest_newref.R

### scPred.R

### SeuratTransferData.R



## Gene Set Enrichment and Pathway Analysis

### gsea.R

### pathway_genes_all.R

### pathway_genes_TNKR.R



## Plot Generation

### CH_heatmaps.R

### heatmap_DE.R

### project_plots.R



## Correlation and Statistical Analyses

### corr_to_exh_Tonly.R

### corr_to_exh.R

### corr_xy.R

### relevant_stats.R


