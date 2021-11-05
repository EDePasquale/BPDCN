library(plyr)
library(dplyr)
library(Seurat)

samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
mycolors_master1=c("#AD82BA", "#88D1D2", "#B4D458", "#CFE1B6", "#D7B063")

for(i in 1:(length(samples_A))-1){
  
  #edit
  S_PT=readRDS(paste0("~/Documents/Projects/Data_files_temp/Donor/", samples_B[i], "/Seurat.rds"))
  CH_PT=read.table(paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_A[i], "/input/CellClassification/", samples_B[i], "-CellClassification.txt"), header=T, sep="\t")
  E_PT=read.table(paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisPeter/01_SignatureScores/", samples_B[i] ,"_Texh.txt"), header=T, sep="\t")
  T_PT=read.csv(paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisMarina/", samples_B[i], "_barcode_results.csv"), header=T)
  result_name=samples_B[i]
  
  # run
  A=S_PT@reductions[["umap"]]@cell.embeddings
  row.names(A)=gsub("-", ".", row.names(A))
  A=cbind(names=row.names(A),A)
  B=cbind(names=row.names(CH_PT),CH_PT)
  
  M_PT=left_join(as.data.frame(A), B)
  M_PT=M_PT[,-(4:8)]
  
  E_PT[,1]=gsub("-", ".", E_PT[,1])
  colnames(E_PT)[1]="names"
  M_PT=left_join(M_PT, E_PT)
  colnames(M_PT)[5]="Exh"
  
  D=T_PT[which(T_PT$RNAannotation != "Not_annotated"),]
  E=T_PT[!(is.na(T_PT$TRB_CDR3_CloneSize)),]
  colnames(E)[1]="strip"
  
  M_PT=cbind(strip=M_PT[,1], M_PT)
  M_PT[,1]=gsub("\\..*", "", gsub("_.*", "", M_PT[,1]))
  
  M_PT=left_join(M_PT, E)
  M_PT=M_PT[,-c(7, 9:34, 36:37)]
  M_PT_2=M_PT
  
  M_PT[is.na(M_PT)] = 0
  
  M_PT_3=M_PT_2[which(!is.na(M_PT_2$TRB_CDR3_CloneSize)),]
  
  
  # clone size and clone prevalence
  M_PT_4=M_PT_3[M_PT_3$Label %in% c("CD4MemoryT", "CD4NaiveT", "CD8MemoryT", "CD8NaiveT"),]
  mean_exh=aggregate(x = M_PT_4$Exh,                # Specify data column
                     by = list(M_PT_4$TRB_CDR3),              # Specify group indicator
                     FUN = mean)  
  mean_CS=aggregate(x = M_PT_4$TRB_CDR3_CloneSize,                # Specify data column
                    by = list(M_PT_4$TRB_CDR3),              # Specify group indicator
                    FUN = mean) 
  table(mean_CS$x)
  round((table(mean_CS$x)/nrow(mean_CS))*100, 2)
  table(M_PT_4$TRB_CDR3_CloneSize)
  round((table(M_PT_4$TRB_CDR3_CloneSize)/nrow(M_PT_4))*100, 2)
  
  # barplot (mean exhaustion per clonotype)
  vals=cbind(CS=mean_CS$x, EXH=mean_exh$x)
  unique(vals[,1])
  x=sort(as.numeric(unique(vals[,1])))
  PROP=round((x/nrow(mean_CS)*100),2)
  y=paste0(x," (", PROP, "%)")
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_exhperclono.pdf"), width = 10, height = 5)
    boxplot(EXH ~ CS, data=vals,  names=y, xlab="Clone Size (Proportion)", ylab = "Mean Exhaustion Score", col=mycolors_master1[i], cex.axis=0.5, las=2)
  dev.off()
  
  # barplot (exhaustion per cell)
  vals=cbind(CS=M_PT_4$TRB_CDR3_CloneSize, EXH=M_PT_4$Exh)
  unique(vals[,1])
  x=sort(as.numeric(unique(vals[,1])))
  PROP=round((x/nrow(M_PT_4)*100),2)
  y=paste0(x," (", PROP, "%)")
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_exhpercell.pdf"), width = 10, height = 5)
    boxplot(EXH ~ CS, data=vals,  names=y, xlab="Clone Size (Proportion)", ylab = "Exhaustion Score", col=mycolors_master1[i], cex.axis=0.5, las=2)
  dev.off()
  
  
  # barplot (mean exhaustion per clonotype) + points
  vals=cbind(CS=mean_CS$x, EXH=mean_exh$x)
  unique(vals[,1])
  x=sort(as.numeric(unique(vals[,1])))
  PROP=round((x/nrow(mean_CS)*100),2)
  y=paste0(x," (", PROP, "%)")
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_exhperclono_points.pdf"), width = 10, height = 5)
    boxplot(EXH ~ CS, data=vals,  names=y, col="white", xlab="Clone Size (Proportion)", ylab = "Mean Exhaustion Score", cex.axis=0.5, las=2)
    stripchart(EXH ~ CS,  
              data=vals,# Data
               method = "jitter", # Random noise
               pch = 16, # Pch symbols
               cex=0.5,
               col = mycolors_master1[i],           # Color of the symbol
               vertical = TRUE,   # Vertical mode
               add = TRUE)        # Add it over
  dev.off()
  
  # barplot (exhaustion per cell) + points
  vals=cbind(CS=M_PT_4$TRB_CDR3_CloneSize, EXH=M_PT_4$Exh)
  unique(vals[,1])
  x=sort(as.numeric(unique(vals[,1])))
  PROP=round((x/nrow(M_PT_4)*100),2)
  y=paste0(x," (", PROP, "%)")
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_exhpercell_points.pdf"), width = 10, height = 5)
  boxplot(EXH ~ CS, data=vals,  names=y, col="white", xlab="Clone Size (Proportion)", ylab = "Exhaustion Score", cex.axis=0.5, las=2)
  stripchart(EXH ~ CS,  
             data=vals,# Data
             method = "jitter", # Random noise
             pch = 16, # Pch symbols
             cex=0.5,
             col = mycolors_master1[i],           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  dev.off()
  
  
  
  foo <- mapvalues(M_PT_4$TRB_CDR3_CloneSize, from=sort(unique(M_PT_4$TRB_CDR3_CloneSize), decreasing=T), to=seq(1, length(unique(M_PT_4$TRB_CDR3_CloneSize))))
  M_PT_4=cbind(M_PT_4, TRB_CDR3_CloneRank=foo)
  
  # correlation - cell level
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_Tcell.pdf"), width = 7, height = 7)
    plot(M_PT_4$Exh, M_PT_4$TRB_CDR3_CloneRank, pch=16, xlab="Clone Rank", ylab="Exhaustion Score", col=mycolors_master1[i], ylim=c(max(foo), 1))
    x=cor.test(M_PT_4$Exh, M_PT_4$TRB_CDR3_CloneSize, method="spearman", exact=FALSE)
    text(max(M_PT_4$Exh), (max(M_PT_4$TRB_CDR3_CloneRank)/1.1), labels=paste0("R = ", round(x[["estimate"]][["rho"]], 4)), pos=2)
    text(max(M_PT_4$Exh), (max(M_PT_4$TRB_CDR3_CloneRank)/1.05), labels=paste0("P = ", formatC(x[["p.value"]], format = "e", digits = 2)), pos=2)
  dev.off()

  # correlation - cell level
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_Tcell_trans.pdf"), width = 7, height = 7)
    plot(M_PT_4$TRB_CDR3_CloneRank, M_PT_4$Exh, pch=16, xlab="Clone Rank", ylab="Exhaustion Score", col=mycolors_master1[i], xlim=c(max(foo), 1))
    x=cor.test(M_PT_4$Exh, M_PT_4$TRB_CDR3_CloneSize, method="spearman", exact=FALSE)
    text(min(M_PT_4$TRB_CDR3_CloneRank), (min(M_PT_4$Exh)/1.3), labels=paste0("R = ", round(x[["estimate"]][["rho"]], 4)), pos=2)
    text(min(M_PT_4$TRB_CDR3_CloneRank), (min(M_PT_4$Exh)/1.05), labels=paste0("P = ", formatC(x[["p.value"]], format = "e", digits = 2)), pos=2)
  dev.off()
  
  
  mean_exh=aggregate(x = M_PT_3$Exh,                # Specify data column
            by = list(M_PT_3$TRB_CDR3),              # Specify group indicator
            FUN = mean)  
  mean_CS=aggregate(x = M_PT_3$TRB_CDR3_CloneSize,                # Specify data column
                     by = list(M_PT_3$TRB_CDR3),              # Specify group indicator
                     FUN = mean) 
  
  
  # box plots for clone size vs exhaustion
  percentile <- function(p) function(x) quantile(x, seq(0, 1, by = 1/p)) 
  quartile <- percentile(2)
  bin <- function(x, q, ...) cut(x, q(x), include.lowest = TRUE, ...)
  vals=cbind(CS=mean_CS$x, EXH=mean_exh$x)
  vals=vals[which(vals[,1] >=3),]
  vals=cbind(vals, QU=bin(vals[,1], quartile))
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_box.pdf"), width = 7, height = 7)
    boxplot(EXH ~ QU, data=vals, names=c("[3-4]","(4-27]"), xlab="Clone Size (Quantiles)", ylab = "Mean Exhaustion Score", col=mycolors_master1[i])
  dev.off()
  
  # correlation - by clone
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_mean.pdf"), width = 7, height = 7)
    plot(mean_exh$x, mean_CS$x, pch=16, xlab="Mean Exhaustion Score", ylab="Clone Size", col=mycolors_master1[i])
    x=cor.test(mean_exh$x, mean_CS$x)
    text(max(mean_exh$x), (max(mean_CS$x)/10), labels=paste0("R = ", round(x[["estimate"]][["cor"]], 4)), pos=2)
    text(max(mean_exh$x), (max(mean_CS$x)/15), labels=paste0("P = ", formatC(x[["p.value"]], format = "e", digits = 2)), pos=2)
  dev.off()
  
  
  mean_exh=aggregate(x = M_PT_3$Exh,                # Specify data column
                     by = list(M_PT_3$Label),              # Specify group indicator
                     FUN = mean)  
  mean_CS=aggregate(x = M_PT_3$TRB_CDR3_CloneSize,                # Specify data column
                    by = list(M_PT_3$Label),              # Specify group indicator
                    FUN = mean) 
  
  # correlation - by cluster
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_mean_lab.pdf"), width = 7, height = 7)
  plot(mean_exh$x, mean_CS$x, pch=16, xlab="Mean Exhaustion Score", ylab="Clone Size", col=mycolors_master1[i])
  x=cor.test(mean_exh$x, mean_CS$x)
  text(mean_exh$x, (mean_CS$x)-1,labels=mean_exh$Group.1, col="darkgray")
  dev.off()
  
  # correlation - by cluster, scaled points by clone size
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corr_xy_mean_lab2.pdf"), width = 7, height = 7)
  plot(mean_exh$x, mean_CS$x, pch=16, xlab="Mean Exhaustion Score", ylab="Clone Size", col=mycolors_master1[i], cex=(0.005*as.numeric(table(M_PT_3$Label))))
  x=cor.test(mean_exh$x, mean_CS$x)
  text(mean_exh$x, (mean_CS$x)-3,labels=mean_exh$Group.1, col="darkgray", cex=(0.002*as.numeric(table(M_PT_3$Label))))
  dev.off()
  ################
  
  S_PT@meta.data[["cellHarmony"]]<-M_PT$Label
  S_PT@meta.data[["exhaustion"]]<-M_PT$Exh
  S_PT@meta.data[["clone_size"]]<-M_PT$TRB_CDR3_CloneSize
  S_PT@meta.data[["clone_norm"]]<-M_PT$TRB_CDR3_CloneSize_Norm
  
  S_PT <- subset(S_PT, subset = (cellHarmony ==  "CD8MemoryT" | cellHarmony == "CD4MemoryT" | cellHarmony == "CD4NaiveT" | cellHarmony == "CD8NaiveT") )
  S_PT <- RunUMAP(S_PT, dims = 1:30)
  
  pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, ".pdf"), width = 10, height = 10)
  par(mar=c(1, 1, 1, 1))
    layout(matrix(1:3, ncol=3, nrow=1, byrow=TRUE))
    DimPlot(S_PT, group.by = "cellHarmony")
    FeaturePlot(S_PT, features = 'exhaustion', cols=c("lightgrey", "red"))
    FeaturePlot(S_PT, features = 'clone_size', cols=c("lightgrey", "blue"))
  dev.off()
  
  # counts=as.data.frame(S_PT@assays[["RNA"]]@counts)
  # which((gsub("-", ".", colnames(counts)) == M_PT$names) == F) #0
  # 
  # cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
  # 
  # corrs=data.frame(genes=row.names(counts)) #blank data frame for the results
  # corrs=as.data.frame(matrix(nrow=nrow(counts), ncol=length(cell_type)))
  # row.names(corrs)=row.names(counts)
  # colnames(corrs)=cell_type
  # for(j in 1:length(cell_type)){
  #   ind=which(M_PT$Label == cell_type[j])
  #   E_temp=M_PT$TRB_CDR3_CloneSize[ind] #only place I need to change for TCR score
  #   C_temp=counts[,ind]
  #   if(!is.null(nrow(C_temp))){
  #     for(i in 1:nrow(C_temp)){
  #       corrs[i,j]=cor(as.numeric(C_temp[i,]), E_temp)
  #     }
  #   }
  # }
  # write.table(corrs, paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_corrs_full.txt"), sep="\t", quote=F)
  # corrs2=corrs[,-10]
  # corrs2[is.na(corrs2)] = 0
  # 
  # tnk=c(1,4,9,16,6)
  # corrs3=corrs2[,tnk]
  # corrs3=corrs3[order(corrs3$CD8MemoryT, decreasing=T),]
  # 
  # library(pheatmap)
  # #
  # pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/", result_name, "_plot.pdf"), width = 5, height = 15)
  # par(mar=c(1, 1, 1, 1))
  # pheatmap(corrs3[1:100,1:5], cluster_rows = F, cluster_cols = F)
  # dev.off()
  # 
  # 
  # x=cbind(as.data.frame(row.names(corrs3)), corrs3$CD8MemoryT)
  # write.table(x, paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/",result_name , "_corrs_CD8.txt"), sep="\t", quote=F, row.names = F)
  # 
}

