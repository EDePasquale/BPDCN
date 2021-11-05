library(ggplot2)
library(umap)
library(Seurat)

reference=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds")
datasets=c("~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds",
           "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds")
samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
samples_C=#a list of samples with different identifiers (for differing names in files)
cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", 
            "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090") #0,1,2,3,4,5 etc 
cluster_colors2=c("#C0007C", "#2a4bd7", "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090",
                 "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff") #0,1,10,11 etc

for(i in 1:length(datasets)){
  
  query=readRDS(datasets[i])
  anchors = FindTransferAnchors(reference = reference, query = query,  reference.reduction = 'pca' )
  
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
  query <- MapQuery(anchorset = anchors, 
                    reference = reference, 
                    query = query,
                    reference.reduction = "pca", 
                    reduction.model = "umap")
  
  file_path=paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_A[i], "/input/CellClassification/", samples_B[i],"-CellClassification.txt")
  labs=read.table(file_path, header=T, sep="\t")
  labs2=labs$Ref.Partition #actually the numeric label
  labs2=labs2-1 #because Seurat
  names(labs2)=row.names(labs)
  names(labs2)=gsub("[.]", "-", names(labs2))
  labs3=labs2[names(query@active.ident)]
  query@meta.data[["CH"]] <- as.factor(labs3)

  
  p2 <- DimPlot(query, cols = cluster_colors[(as.numeric(names(table(labs2)))+1)], reduction = "ref.umap", group.by = "CH", label = F,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle(samples_C[i])
  pdf(file = paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_B[i], "_project_plots.pdf"), width = 6, height = 6)
  par(mar=c(2, 2, 2, 2))
    p2
  dev.off()
  
}


