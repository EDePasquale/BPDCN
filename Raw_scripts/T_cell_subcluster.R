#######################
#                     #
# Sub Cluster T-cells #
#   (Create Plots)    #
#                     #
#######################

library(Seurat)
library(dplyr)

# Pull in new, more granular clustering and seperate out the 2 clusters I want to put into the old clustering scheme
new=read.table("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_1.0_ccrem/Cluster.txt", sep="\t", header=T)
new_redu=new[which(new$x %in% c(0,8)),]

# Pull in the old clustering scheme and seperate out the original CD8 cluster
old=read.table("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/FinalGroups.txt", sep="\t")
old_redu=old[which(old$V2 %in% c(1)),]

# Harmonize cell names
new_redu[,1]=gsub("-", ".", new_redu[,1])

# Check that these two seperated out cell subsets are (nearly) identical
length(intersect(old_redu[,1], new_redu[,1]))/nrow(new_redu) #>99% overlap (1890/1897)

# Get a list of the cell identities that need to change. In this case, it is the original cluster 0 cells that are in the
# new granular cluster 8
new_8=new[which(new$x %in% c(8)),1]

# Pull in Master Seurat object to modify
immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")

# Plot with and without cluster labels
pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/UMAP_Seurat_Labs.pdf"), width = 6, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/UMAP_Seurat_No_Labs.pdf"), width = 6, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE)
dev.off()

# Add in a new level (16) for the identities factor. This needs to be at the end of the list so other numbers aren't changed
levels(immune.combined@active.ident)=c(levels(immune.combined@active.ident), 16)
immune.combined@active.ident[which(names(immune.combined@active.ident) %in% new_8)]<-as.factor(16)

# Plot with and without cluster labels
pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/UMAP_Seurat_TCELL_Labs.pdf"), width = 6, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/UMAP_Seurat_TCELL_No_Labs.pdf"), width = 6, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE)
dev.off()

# Save this Seurat object as Seurat_int_0.5_TCELLSPLIT.rds to designate it as a modified master
saveRDS(immune.combined, file = paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5_TCELLSPLIT.rds"))

