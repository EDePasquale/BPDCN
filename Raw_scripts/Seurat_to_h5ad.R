#################################
#                               #
# Convert Seurat Object to h5ad #
#                               #
#################################

install.packages('Seurat')
devtools::install_github('satijalab/seurat-data')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")


# Load libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# Load Seurat object
immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")

# Save as a h5ad
SaveH5Seurat(immune.combined, filename = "immune.combined.h5Seurat")
Convert("immune.combined.h5Seurat", dest = "h5ad")

#NOTE: needs R 4

#_______________________________

BiocManager::install("zellkonverter")

#NOTE: needs R 4

#_______________________________

BiocManager::install("LoomExperiment")
devtools::install_github("cellgeni/sceasy")
install.packages('reticulate')

library(Seurat)
library(sceasy)
library(reticulate)
library(LoomExperiment)

# Load Seurat object
immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")

setwd("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/")
sceasy::convertFormat(immune.combined, from="seurat", to="anndata",
                      outFile='immune.combined.h5ad')

#NOTE: needs non-system install of Python (miniconda?)
