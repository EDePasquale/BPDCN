
  # to install the SeuratData package see https://github.com/satijalab/seurat-data
  # library(SeuratData)
  # data("pbmc3k")
  
library(Seurat)

  # for demonstration, split the object into reference and query
  # pbmc.reference <- pbmc3k[, 1:1350]
  # pbmc.query <- pbmc3k[, 1351:2700]
  # 
  # # perform standard preprocessing on each object
  # pbmc.reference <- NormalizeData(pbmc.reference)
  # pbmc.reference <- FindVariableFeatures(pbmc.reference)
  # pbmc.reference <- ScaleData(pbmc.reference)
  # 
  # pbmc.query <- NormalizeData(pbmc.query)
  # pbmc.query <- FindVariableFeatures(pbmc.query)
  # pbmc.query <- ScaleData(pbmc.query)

reference=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")
datasets=c("SAMPLE_NAME", 
           "SAMPLE_NAME")
for(i in datasets){
  query=readRDS(paste0("~/Documents/Projects/Data_files_temp/Donor/", i, "/Seurat.rds"))
  
  # find anchors
  anchors <- FindTransferAnchors(reference = reference, query = query)
  
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = reference$seurat_clusters)
  query <- AddMetaData(object = query, metadata = predictions)
  
  assign(paste0("res_", i), table(query@meta.data[["predicted.id"]]))
}



library(dplyr)
PlotFreq=data.frame(x=0:15, y=0:15)

temp1=as.data.frame(res_SAMPLE_NAME)
temp1[,1]=as.numeric(as.character(temp1[,1]))
colnames(temp1)[2]="SAMPLE_NAME"
PlotFreq=left_join(PlotFreq, temp1, by=c("x"="Var1"))

temp1=as.data.frame(res_SAMPLE_NAME)
temp1[,1]=as.numeric(as.character(temp1[,1]))
colnames(temp1)[2]="SAMPLE_NAME"
PlotFreq=left_join(PlotFreq, temp1, by=c("x"="Var1"))


PlotFreq=PlotFreq[,-(1:2)]
PlotFreq[is.na(PlotFreq)] = 0


PlotFreqNorm.mat <- sweep(PlotFreq, 2, colSums(PlotFreq), "/")*100

setwd("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/")
dir.ch <- "210722_SeuratTransData_10BPDCN_ccrem"
dir.create(dir.ch)
setwd(paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/", dir.ch))

write.table(PlotFreqNorm.mat, "PlotFreqNorm.txt", sep="\t", quote=F)

##### reorder to match peter's ####
ordering_vec=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,4,6)
#ordering_vec=c(11,7,2,6,4,13,12,15,10,9,1,14,8,0,3,5)
renamed=PlotFreqNorm.mat
row.names(renamed)=1:16
renamed=renamed[ordering_vec,]

clus_names=as.data.frame(cbind(0:15, c("CD8T", "B", "MidEryth", "CTL", "Mono", "NK","LateEryth", "EarlyEryth", "CD4T", "PreB",
                                       "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC")))
CellTypeCol.ch=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a")
clus_names_col=cbind(clus_names, CellTypeCol.ch)

renamed2=clus_names_col$CellTypeCol.ch
renamed2=renamed2[ordering_vec]

pdf(paste0("3_CellTypeFrequencies.pdf"), width = 8, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(as.matrix(renamed[nrow(renamed):1,]), col = rev(renamed2), xaxt = "n", ylab = "Population frequency (%)", border = NA)

axis(side = 1, at = seq(1,ncol(renamed))*1.2-0.5, labels = colnames(renamed), las = 2)
legend(x = ncol(renamed)*1.2+0.5, y = 100, legend = clus_names$V2[ordering_vec], fill = renamed2, bty = "n", border = NA)
dev.off()


