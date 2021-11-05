library(Seurat)
library(scPred)
library(dplyr)
library(gplots)

#models_2=c("svmRadial", "svmLinear", "bayesglm", "glmboost", "mda", "rmda", "gbm", "fda", "rf")

reference=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/Seurat_int_0.5.rds")
reference <- getFeatureSpace(reference, "seurat_clusters")
reference <- trainModel(reference)
#reference <- trainModel(reference, model = "gamSpline", reclassify = c("12", "15"))
get_probabilities(reference) %>% head()
get_scpred(reference)
plot_probabilities(reference)

datasets=c("SAMPLE_NAME", 
           "SAMPLE_NAME", 
           "BPDCN181128")

for(i in datasets){
  query=readRDS(paste0("~/Documents/Projects/Data_files_temp/Donor/", i, "/Seurat.rds"))
  query <- scPredict(query, reference, threshold=0)
  query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
  assign(paste0("res_", i),  table(query@meta.data[["scpred_prediction"]]))
}

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
dir.ch <- "210722_scPred_10BPDCN_ccrem"
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







#########################################
query=readRDS("~/Documents/Projects/Data_files_temp/Donor/BPDCN190711/Seurat.rds")
query <- scPredict(query, reference)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

query=readRDS("~/Documents/Projects/Data_files_temp/Donor/BPDCN181128/Seurat.rds")
query <- scPredict(query, reference)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

query=readRDS("~/Documents/Projects/Data_files_temp/Donor/BPDCN180329/Seurat.rds")
query <- scPredict(query, reference)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

preds=cbind(query@meta.data[["scpred_0"]],
            query@meta.data[["scpred_1"]],
            query@meta.data[["scpred_2"]],
            query@meta.data[["scpred_3"]],
            query@meta.data[["scpred_4"]],
            query@meta.data[["scpred_5"]],
            query@meta.data[["scpred_6"]],
            query@meta.data[["scpred_7"]],
            query@meta.data[["scpred_8"]],
            query@meta.data[["scpred_9"]],
            query@meta.data[["scpred_10"]],
            query@meta.data[["scpred_11"]],
            query@meta.data[["scpred_12"]],
            query@meta.data[["scpred_13"]],
            query@meta.data[["scpred_14"]],
            query@meta.data[["scpred_15"]])

breaks=seq(0, #start point of color key
           1,  #end point of color key
           by=0.01) #length of sub-division
mycol=colorpanel(n=length(breaks)-1,low="white", high="blue") #heatmap colors
heatmap.2(preds, col=mycol, trace="none", Rowv = FALSE, Colv = FALSE, dendrogram = 'none')

r_avg=AverageExpression(reference, "RNA")
q_avg=AverageExpression(query, "RNA")
q=as.matrix(query@assays[["RNA"]]@counts)
r_avg_matrix=r_avg[["RNA"]]

genes=intersect(row.names(r_avg_matrix), row.names(q))
r_avg_matrix2=r_avg_matrix[genes,]
q2=q[genes,]

corr=cor(r_avg_matrix2, q2)
pheatmap(corr)

P=as.data.frame(query@active.ident)
identsP=cbind(row.names(P), P)
pheatmap(corr, cluster_cols = F, cluster_rows = F, show_colnames=F, annotation_col=identsP)


