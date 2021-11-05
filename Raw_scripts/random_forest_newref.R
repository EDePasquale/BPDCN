################################
#                              #
# Random Forest Classification #
#      Erica DePasquale        #
#       June 22, 2021          #
#                              #
################################

#MODS on lines: 16-19, 22, 29, 32, 77-81, 102-111
# Load libraries
library(Seurat)
library(dplyr)
#library(randomcoloR)
library(gplots)

library(gdata)
library(randomForest)
library(Seurat)
library(tidyverse)
library(gplots)
library(data.table)

# Set working directory
# setwd("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5/")
# dir.ch <- "210628_RandomForest_10BPDCN"
# dir.create(dir.ch)
# setwd(paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5/", dir.ch))
setwd("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/")
dir.ch <- "210726_RandomForest_10BPDCN_ccrem"
dir.create(dir.ch)
setwd(paste0("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/", dir.ch))

# Pull in functions
source("~/Dropbox (Partners HealthCare)/Single-cell_BPDCN/AnalysisPeter/SingleCell_BPDCN_Functions.R")

#==================================
# Load and process reference data
#==================================

# Read in Seurat Integrated data
#immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5/Seurat_int_0.5.rds")
immune.combined=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_TCELLSPLIT.rds")

# Cluster - cell type associations
# clus_names=as.data.frame(cbind(0:16, c("CD8T", "B", "MidEryth", "CTL", "NK", "Mono", "LateEryth", "EarlyEryth", "CD4T", "PreB",
#                                        "cDC", "HSCProg", "ProBPreB", "ncMono", "Plasma", "ProB", "pDC")))
clus_names=as.data.frame(cbind(0:16, c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK","LateEryth", "EarlyEryth", "CD4MemoryT", "PreB",
                                       "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC", "CD8NaiveT")))


# Find marker genes
M.markers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct=0.25, logfc.threshold = 0.25)
M.markers <- left_join(M.markers, clus_names, by=c("cluster" = "V1"))
top50 <- M.markers %>% group_by(V2) %>% top_n(n = 50, wt = avg_log2FC)
selected.genes <- unique( unlist(top50$gene, use.names = F) ) #just unique genes

# Make Seurat's count matrix easier to access (log, transcript per 10K)
#bm.cm <- as.matrix( GetAssayData(object = immune.combined, slot = "data") )
bm.cm=as.matrix(immune.combined@assays[["integrated"]]@data)

# Plot named clusters
#CellTypeCol.ch = distinctColorPalette(nrow(clus_names))
# CellTypeCol.ch=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#BE8A66", "#ffee33", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
#                 "#1d6914", "#000000", "#708090", "#29d0d0", "#575757", "#8126c0", "#81c57a")
CellTypeCol.ch=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090")




clus_names_col=cbind(clus_names, CellTypeCol.ch)

DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.5, cols = CellTypeCol.ch) + #CellTypeCol.ch is colors per cell type
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")



#==================================
# Load and process patient data
#==================================

# Function to load and process 10X data
# Process_10X <- function(RData, patient_id) {
#   load(RData, verbose = T)
#   stopifnot(all(colnames(CM.dgm) == stats.dt$cell))
#   CM.df <- as.matrix( CM.dgm[rownames(bm),] )
#   colnames(CM.df) <- colnames(CM.df) %>% str_replace("-1", str_c("-", patient_id, ".1")) %>%
#     str_replace("-2", str_c("-", patient_id, ".2")) %>% str_replace("-3", str_c("-", patient_id, ".3")) %>%
#     str_replace("-4", str_c("-", patient_id, ".4"))
#   seu <- CreateSeuratObject(counts = CM.df, project = patient_id)
#   seu <- NormalizeData(seu) # (log, transcript per 10K)
#   seu$replicate <- cutf(colnames(seu), d = "-", f = 2)
#   seu$tech <- "TenX"
#   seu$my.tSNE.x <- stats.dt$tSNEx
#   seu$my.tSNE.y <- stats.dt$tSNEy
#   return(seu)
# }

# Create tibble with information on BPDCN samples that were analyzed with 10X
TenX_samples.tib <- tribble(~Patient_ID, ~Rds, ~Note,
                            "SAMPLE_NAME", "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds", "Erica's project",
                            "SAMPLE_NAME", "~/Documents/Projects/Data_files_temp/Donor/SAMPLE_NAME/Seurat.rds", "Erica's project",


# Create list of processed Seurat objects
seu_10X.ls <- vector(mode = "list", length = nrow(TenX_samples.tib))
for (x in 1:nrow(TenX_samples.tib)) {
  #seu_10X.ls[[x]] <- Process_10X(RData = TenX_samples.tib[[x,"RData"]],
  #                               patient_id = TenX_samples.tib[[x,"Patient_ID"]])
  seu_10X.ls[[x]] <- readRDS(as.character(TenX_samples.tib[x,2]))
}
names(seu_10X.ls) <- TenX_samples.tib$Patient_ID
#seu.ls <- c(Pt9Dx = list(Pt9Dx.seu), seu_10X.ls)
#cm.ls <- lapply(seu_10X.ls, function(x) as.matrix(GetAssayData(x), slot = "data"))
cm.ls <- lapply(seu_10X.ls, function(x) as.matrix(x@assays[["RNA"]]@data))

# Merge all gene expression data & check if it all makes sense.
#all.mat <- cbind(bm.cm, do.call(cbind, cm.ls))
#table(cutf(colnames(all.mat), d = "-", f = 2))
#table(round(colSums(exp(all.mat)-1)))
all(colnames(bm.cm) == names(immune.combined@active.ident))

# #select_genes_common=intersect(intersect(intersect(row.names(bm.cm), row.names(cm.ls[[1]])), row.names(cm.ls[[2]])), row.names(cm.ls[[3]]))
#    #limit to genes found in the BPDCN data, because of genes filtered out in previous steps. Need to train on only the common genes
# select_genes_common=intersect(
#                       intersect(
#                         intersect(
#                           intersect(
#                             intersect(row.names(bm.cm),
#                               row.names(cm.ls[[1]])),
#                             row.names(cm.ls[[2]])),
#                           row.names(cm.ls[[3]])),
#                         row.names(cm.ls[[4]])),
#                       row.names(cm.ls[[5]]))

select_genes_common=Reduce(intersect, list(row.names(bm.cm),
                       row.names(cm.ls[[1]]),
                       row.names(cm.ls[[2]]),
                       row.names(cm.ls[[3]]),
                       row.names(cm.ls[[4]]),
                       row.names(cm.ls[[5]]),
                       row.names(cm.ls[[6]]),
                       row.names(cm.ls[[7]]),
                       row.names(cm.ls[[8]]),
                       row.names(cm.ls[[9]]),
                       row.names(cm.ls[[10]])))

#==================================
# Build classifier
#==================================

# Plant classification trees
set.seed(123)

#rf <- randomForest(x = t(bm.cm[selected.genes,]),    # matrix of predictors
rf <- randomForest(x = t(bm.cm[select_genes_common,]),    # matrix of predictors
                   y = immune.combined@active.ident, # response vector
                   sampsize = rep(50, length(levels(immune.combined@active.ident))),
                   ntree = 1000,
                   do.trace = 100)

# Plot confusion matrix (based on out-of-bag data), with colors normalized for the total cell number in each population
Conf.mat <- rf$confusion[, rownames(rf$confusion)]
NormConf.mat <- Conf.mat / rowSums(Conf.mat)

pdf(paste0("1_ConfusionMatrix.pdf"), width = 8, height = 8)
heatmap.2(NormConf.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.mat)), rowsep = c(0, nrow(NormConf.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.mat)) / sum(Conf.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.mat), each=nrow(NormConf.mat)),
                          rep(ncol(NormConf.mat):1, nrow(NormConf.mat)), Conf.mat))
dev.off()

#==================================
# Five-fold cross-validation
#==================================

# Split dataset in five parts
cv <- split(colnames(bm.cm), rep(1:5, 1E6)[1:ncol(bm.cm)])

# Build five forests, each with 4/5 of the data
rf.cv <- lapply(cv, function(n) {
  set.seed(123)
  randomForest(x = t(bm.cm[selected.genes,setdiff(colnames(bm.cm), n)]),
               y = immune.combined@active.ident[! colnames(bm.cm) %in% n],
               sampsize = sapply(table(immune.combined@active.ident[! colnames(bm.cm) %in% n]), min, 50), # check that it's not too low
               ntree = 1000,
               do.trace = 100)
})

# Predict the sets that were not used for training
rf.cv.predict.prob <- lapply(rf.cv, function(rf) {
  predict(rf, t(bm.cm[selected.genes, setdiff(colnames(bm.cm), names(rf$y))]), type = "prob")
})

# Maximum probability
rf.cv.predict <- lapply(rf.cv.predict.prob, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))  # maximum sore
  names(y) <- rownames(x)
  y
})
head(rf.cv.predict[[1]])

# Confusion matrix
Conf.cv.mat <- table(immune.combined@active.ident[unlist(lapply(rf.cv.predict, names))], unlist(rf.cv.predict))
NormConf.cv.mat <- Conf.cv.mat / rowSums(Conf.cv.mat)

pdf(paste0("2_CrossValidation.pdf"), width = 8, height = 8)
heatmap.2(NormConf.cv.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.cv.mat)), rowsep = c(0, nrow(NormConf.cv.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.cv.mat)) / sum(Conf.cv.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.cv.mat), each=nrow(NormConf.cv.mat)),
                          rep(ncol(NormConf.cv.mat):1, nrow(NormConf.cv.mat)), Conf.cv.mat))
dev.off()

#==================================
# Classify cells from patients
#==================================

# Reduce to common genes
bm.cm2 <- bm.cm[select_genes_common,]
cm.ls[[1]] <- cm.ls[[1]][select_genes_common,]
cm.ls[[2]] <- cm.ls[[2]][select_genes_common,]
cm.ls[[3]] <- cm.ls[[3]][select_genes_common,]
cm.ls[[4]] <- cm.ls[[4]][select_genes_common,]
cm.ls[[5]] <- cm.ls[[5]][select_genes_common,]
cm.ls[[6]] <- cm.ls[[6]][select_genes_common,]
cm.ls[[7]] <- cm.ls[[7]][select_genes_common,]
cm.ls[[8]] <- cm.ls[[8]][select_genes_common,]
cm.ls[[9]] <- cm.ls[[9]][select_genes_common,]
cm.ls[[10]] <- cm.ls[[10]][select_genes_common,]

# Predict cell types in each of the BPDCN samples (ties are broken at random)
predictions.mat.ls <- lapply(c(BM = list(bm.cm2), cm.ls), function(x) predict(rf, t(x[select_genes_common,]), type = "prob"))
#predictions.mat.ls <- lapply(c(BM = list(bm.cm), cm.ls), function(x) predict(rf, t(x[intersect(row.names(x), selected.genes),]), type = "prob"))

# Maximum probability
CellTypes.ls <- lapply(predictions.mat.ls, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))
  names(y) <- rownames(x)
  y
})

for(i in 1:length(CellTypes.ls)){
  write.table(CellTypes.ls[[i]], names(CellTypes.ls[i]), sep="\t", col.names = F, quote=F)
}


#TODO not working!
# Plot the clustered cell type frequencies in BM (per donor) and the predicted cell type frequencies in BPDCN
BMfreq.mat <- do.call(cbind, lapply(split(immune.combined@meta.data, f = cutf(immune.combined@meta.data$orig.ident, d = "\\.")), function(x) table(x$CellType)))
PredictFreq.mat <- do.call(cbind, lapply(CellTypes.ls, table))
# Percent doublets. For healthy donor doublets, see bm$CellType
#apply(PredictFreq.mat, 2, function(x) x["Doublets"] / sum(x) * 100)
#PlotFreq.mat <- cbind(BMfreq.mat, PredictFreq.mat[,-1])[-match("Doublets", rownames(PredictFreq.mat)),]
#PlotFreq.mat <- cbind(BMfreq.mat, PredictFreq.mat[,-1])
PlotFreq.mat<-PredictFreq.mat

# Normalize to 100
PlotFreqNorm.mat <- sweep(PlotFreq.mat, 2, colSums(PlotFreq.mat), "/")*100

write.table(PlotFreqNorm.mat, "PlotFreqNorm.txt", sep="\t", quote=F)

##### reorder to match peter's ####
ordering_vec=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,4,6)
#ordering_vec=c(11,7,2,6,4,13,12,15,10,9,1,14,8,0,3,5)
renamed=PlotFreqNorm.mat
row.names(renamed)=1:17
renamed=renamed[ordering_vec,]
renamed2=clus_names_col$CellTypeCol.ch
renamed2=renamed2[ordering_vec]

pdf(paste0("3_CellTypeFrequencies.pdf"), width = 8, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(renamed[nrow(renamed):1,], col = rev(renamed2), xaxt = "n", ylab = "Population frequency (%)", border = NA)

axis(side = 1, at = seq(1,ncol(renamed))*1.2-0.5, labels = colnames(renamed), las = 2)
legend(x = ncol(renamed)*1.2+0.5, y = 100, legend = clus_names$V2[ordering_vec], fill = renamed2, bty = "n", border = NA)
dev.off()



# pdf(paste0("3_CellTypeFrequencies.pdf"), width = 8, height = 6)
# par(mar = c(8,4,8,12), xpd = T)
# 
# barplot(PlotFreqNorm.mat[nrow(PlotFreqNorm.mat):1,], col = rev(clus_names_col$CellTypeCol.ch), xaxt = "n", ylab = "Population frequency (%)", border = NA)
# #barplot(PlotFreqNorm.mat, col = rev(clus_names_col$CellTypeCol.ch), xaxt = "n", ylab = "Population frequency (%)", border = NA)
# 
# 
# axis(side = 1, at = seq(1,ncol(PlotFreqNorm.mat))*1.2-0.5, labels = colnames(PlotFreqNorm.mat), las = 2)
# legend(x = ncol(PlotFreqNorm.mat)*1.2+0.5, y = 100, legend = clus_names$V2, fill = clus_names_col$CellTypeCol.ch, bty = "n", border = NA)
# dev.off()
# 




#TODO finish this

#==================================
# Project cell types & save
#==================================

# Coordinates of normal cells (change UMAP_1 and UMAP_2 to facilitate workflow)
bm.umap <- data.frame(immune.combined@reductions$umap@cell.embeddings)

# Project and save the BPDCN samples
for ( Patient_ID in names(cm.ls) ) {
  #Patient_ID <- names(predictions.mat.ls)[2]
  message(Patient_ID)
  
  # Correlate BPDCN cell prediction scores to BM prediction scores
  cor.mat <- cor(t(predictions.mat.ls[["BM"]]), t(predictions.mat.ls[[Patient_ID]]))
  cor.mat[1:3,1:3]
  cor.max <- apply(cor.mat, 2, which.max)
  cor.id <- rownames(predictions.mat.ls[["BM"]])[cor.max]
  
  # Plot single nearest cells
  bpdcn.project.umap <- data.frame(row.names = colnames(cor.mat),
                                   project.umap.x = bm.umap[cor.id, "UMAP_1"],
                                   project.umap.y = bm.umap[cor.id, "UMAP_2"],
                                   CellType = CellTypes.ls[[Patient_ID]],
                                   CellTypeCol = popcol.df[as.character(CellTypes.ls[[Patient_ID]]),"hex"])
  
  pdf(paste0(Patient_ID, "_cor.predict.pdf"), width = 6, height = 6)
  par(mar=c(4, 4, 4, 4))
  
  plotTSNE(bm.umap, cex = 0.3, col = "#DDDDDD", main = paste(Patient_ID, "on normal BM (grey)"))
  points(bpdcn.project.umap[,1:2], pch = 16, cex = 0.3, col = bpdcn.project.umap$CellTypeCol)
  
  # How does it look for the 2, 5, 10 and 20 most similar cells?
  for (n in c(2, 5, 10, 20)) {
    #n <- 2
    plotTSNE(bm.umap, cex = 0.3, col = "#DDDDDD", main = paste(n, "nearest cells"))
    
    cor.n <- t(apply(cor.mat, 2, function(x) rownames(cor.mat)[ order(x, decreasing = T)[1:n] ] ))
    
    bpdcn.project.umap.n <- data.frame(row.names = rownames(cor.n),
                                       project.umap.x = apply(cor.n, 1, function(x) mean(bm.umap[x,"UMAP_1"])),
                                       project.umap.y = apply(cor.n, 1, function(x) mean(bm.umap[x,"UMAP_2"])),
                                       CellType = CellTypes.ls[[Patient_ID]],
                                       CellTypeCol = popcol.df[as.character(CellTypes.ls[[Patient_ID]]),"hex"])
    
    points(bpdcn.project.umap.n[,1:2], pch = 16, cex = 0.3, col = bpdcn.project.umap$CellTypeCol)
  }
  
  dev.off()
  
  # Save Seurat object will all relevant data
  s <- seu.ls[[Patient_ID]]
  
  # Check
  stopifnot(all(rownames(bpdcn.project.umap) == colnames(s)))
  stopifnot(all(rownames(predictions.mat.ls[[Patient_ID]]) == colnames(s)))
  
  # Add prediction scores for each cell type
  for (x in colnames(predictions.mat.ls[[Patient_ID]])) {
    s <- AddMetaData(s, predictions.mat.ls[[Patient_ID]][,x], col.name = paste0("Predict.", x))
  }
  
  # Add other metadata, exclude doublets, save
  s$CellType <- bpdcn.project.umap$CellType
  s@active.ident <- s$CellType
  s$project.umap.x <- bpdcn.project.umap$project.umap.x
  s$project.umap.y <- bpdcn.project.umap$project.umap.y
  s@commands <- list()
  s <- subset(s, subset = CellType != "Doublets")
  saveRDS(s, file = paste0(Patient_ID, "_Seurat_Predict.rds"))
  
}


#==================================
# Save normal BM object
#==================================

# Add prediction scores to bm & save
for (x in colnames(predictions.mat.ls[["BM"]])) {
  immune.combined <- AddMetaData(immune.combined, predictions.mat.ls[["BM"]][,x], col.name = paste0("Predict.", x))
}

# Some additional wrangling
immune.combined <- subset(immune.combined, subset = CellType != "Doublets")
immune.combined <- DietSeurat(immune.combined)

saveRDS(immune.combined, file = paste0("BM_Seurat_Predict.rds"))






clus_names=as.data.frame(cbind(0:15, c("CD8T", "B", "MidEryth", "CTL", "Mono", "NK","LateEryth", "EarlyEryth", "CD4T", "PreB",
                                       "ProB", "HSCProg", "cDC", "ncMono", "Plasma",  "pDC")))
CellTypeCol.ch=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a")
clus_names_col=cbind(clus_names, CellTypeCol.ch)

setwd("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem")
PlotFreqNorm.mat=read.table("update_nums_5.txt", sep="\t", header=T, row.names = 1)
for(i in 1:ncol(PlotFreqNorm.mat)){
  ssum=sum(PlotFreqNorm.mat[,i])
  PlotFreqNorm.mat[,i]=PlotFreqNorm.mat[,i]/ssum
}
ordering_vec=c(12,8,3,7,5,14,13,16,11,10,2,15,9,1,4,6)
renamed=PlotFreqNorm.mat
row.names(renamed)=1:16
renamed=renamed[ordering_vec,]
renamed2=clus_names_col$CellTypeCol.ch
renamed2=renamed2[ordering_vec]
