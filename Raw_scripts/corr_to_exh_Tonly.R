library(dplyr)
library(Seurat)

i=5

samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)

#edit
S_PT=readRDS(paste0("~/Documents/Projects/Data_files_temp/Donor/", samples_B[i], "/Seurat.rds"))
CH_PT=read.table(paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_A[i], "/input/CellClassification/", samples_B[i], "-CellClassification.txt"), header=T, sep="\t")
E_PT=read.table(paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisPeter/01_SignatureScores/", samples_B[i], "_Texh.txt"), header=T, sep="\t")
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
M_PT=M_PT[,-(7:31)]

M_PT[is.na(M_PT)] = 0

S_PT@meta.data[["cellHarmony"]]<-M_PT$Label
S_PT@meta.data[["exhaustion"]]<-M_PT$Exh
S_PT@meta.data[["clone_size"]]<-M_PT$TRB_CDR3_CloneSize
S_PT@meta.data[["clone_norm"]]<-M_PT$TRB_CDR3_CloneSize_Norm

S_PT <- subset(S_PT, subset = (cellHarmony ==  "CD8MemoryT" | cellHarmony == "CD4MemoryT" | cellHarmony == "CD4NaiveT" | cellHarmony == "CD8NaiveT") )
#DimPlot(S_PT, group.by = "cellHarmony")

S_PT <- RunUMAP(S_PT, dims = 1:30)

pdf(file = paste0("~/Dropbox (Partners HealthCare)/More_BPDCN/AnalysisErica/umaps/", result_name, "_Tonly.pdf"), width = 10, height = 10)
par(mar=c(1, 1, 1, 1))
layout(matrix(1:3, ncol=3, nrow=1, byrow=TRUE))
DimPlot(S_PT, group.by = "cellHarmony")
FeaturePlot(S_PT, features = 'exhaustion', cols=c("lightgrey", "red"))
FeaturePlot(S_PT, features = 'clone_size', cols=c("lightgrey", "blue"))
dev.off()







# counts=as.data.frame(S_PT@assays[["RNA"]]@counts)
# #counts=rbind(M_PT$Label, counts)
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
