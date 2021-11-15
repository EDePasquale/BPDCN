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


## Classification

### random_forest_newref.R
This script was used for performing random forest classification on the BPDCN samples using the integrated control samples as the reference. This process involves building a classifier, doing 5-fold cross validaton, classifying BPDCN cells, writing results, and plotting the results in a stacked bar plot.

### scPred.R
This script runs scPred on the BPDCN samples using the integrated control samples as the reference. Various models were tested, but the function defaults were used in the paper. This process involves running scPred and plotting the results in a stacked bar plot.

### SeuratTransferData.R
This script was used for 1) integrating 2 older healthy control samples with 3 younger control samples that were previously integrated, with cells classified, using Seurat Integration and 2) classifying BPDCN patient cells. This process involves running Seurat's FindTransferAnchors() and TransferData() functions followed by plotting the results in a stacked bar plot for BPDCN samples.


## Gene Set Enrichment and Pathway Analysis

### gsea.R
This script calls a custom GSEA function (https://bioinformaticsbreakdown.com/how-to-gsea/) that combines 2 independent GSEA implementations (fgsea() and gage() functions) and returns pathways significantly enriched with both methods. This script also generates a heatmap of enrichment scores, as significant pathways by cell type/sample.

### pathway_genes_all.R
This script reads in a list of sample names and a geneset file downloaded from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=H) to generate a heatmap of genes from this signature as rows and cell types/samples as columns. 

### pathway_genes_TNKR.R
This script is the same as above with two modifications: 1) The cell populations are limited to T and NK cells for the purposes of generating figures for the paper, and 2) The genes are limited to those in the gene set that are significantly differentially expressed between controls and any BPDCN populations with a minimum aboslute expression difference threshold. The unfiltered version was also included in the manuscript as a Supplemental Figure.


## Plot Generation

### CH_heatmaps.R

### heatmap_DE.R

### project_plots.R

### T_cell_subcluster.R
This script was used for generating plots related to further subclustering of the T-cells in the integrated control sample data. It takes in the original integrated Seurat object as well as higher resolution clustering from Seurat, then applies the more granular T-cell clusters to the integrated object, generates plots, and finally saves the new object. 


## Correlation and Statistical Analyses

### corr_to_exh_Tonly.R

### corr_to_exh.R

### corr_xy.R

### relevant_stats.R


