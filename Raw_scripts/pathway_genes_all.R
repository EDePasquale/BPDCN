


library(data.table)
library(dplyr)
library(pheatmap)
library(randomcoloR)


IFN_genes=read.table("~/Documents/Projects/Data_files_temp/Donor/GSEA/geneset_IFNG.txt", sep="\t") #IFN alpha
IFN_genes=c(IFN_genes[-(1:2),])

reference=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds")


samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
samples_ord=#a list of samples with different identifiers (for differing names in files)
samples_no=c(1,2,3,4,5)

cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
cell_type_ord=c("N.CD4NaiveT", "K.B", "C.MidEryth", "P.CD8MemoryT", "E.Mono", "Q.NK", "D.LateEryth", "B.EarlyEryth", "M.CD4MemoryT", "J.PreB", "I.ProB", "A.HSCProg", "G.cDC", "F.ncMono", "L.Plasma", "H.pDC", "O.CD8NaiveT")
cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090") #0,1,2,3,4,5 etc


annos=as.data.frame(matrix(nrow=2, ncol=0))
master_avgs=as.data.frame(matrix(nrow=length(IFN_genes), ncol=1))
master_avgs[,1]=IFN_genes
colnames(master_avgs)="genes"

#########
# REF
#########

labs3=reference@active.ident
labs2=as.numeric(as.character(labs3))
names(labs2)=names(labs3)
labs2=labs2+1
temp1=labs2
map = setNames(cell_type, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
labs2[] <- map[unlist(labs2)]
labs=as.data.frame(cbind(temp1, labs2))
colnames(labs)=c("Cluster_num", "Cluster_name")
row.names(labs)=gsub("-", ".", row.names(labs))

for(j in 1:length(cell_type)){
  data=reference@assays[["RNA"]]@data[which(row.names(reference@assays[["RNA"]]@data) %in% IFN_genes),]
  #data=reference@assays[["RNA"]]@counts[which(row.names(reference@assays[["RNA"]]@counts) %in% IFN_genes),]
  data=as.data.frame(data)
  rnames=row.names(data)
  colnames(data)=gsub("-", ".", colnames(data))
  
  #subset to cell type
  cell_ind=row.names(labs)[which(labs$Cluster_num==j)]
  if(length(cell_ind)!=0){ #only move forward if cells exist
    data=data[,which(colnames(data) %in% cell_ind)]
    data=as.matrix(data)
    data2=matrix(as.numeric(data),    # Convert to numeric matrix
                 ncol = ncol(data))
    row.names(data2)=rnames
    colnames(data2)=colnames(data)
    temp=cbind(row.names(data2),as.data.frame(rowMeans(data2)))
    colnames(temp)=c("genes", paste0("R_", j))
    master_avgs=left_join(master_avgs, temp, by=c("genes"="genes"))
    annos=cbind(annos, rbind("0.Reference", cell_type_ord[j]))
    rm(data2)
  }
  rm(data)
}
rm(labs)

#########
# QUERY
#########

for(i in 1:length(samples_A)){
  HM=readRDS(paste0("~/Documents/Projects/Data_files_temp/Donor/", samples_B[i], "/Seurat.rds"))
  labs=read.table(paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_A[i], "/input/CellClassification/", samples_B[i], "-CellClassification.txt"), sep="\t")
  labs=labs[,-(1:4)]
  colnames(labs)=c("Cluster_num", "Cluster_name")
  for(j in 1:length(cell_type)){
    data=HM@assays[["RNA"]]@data[which(row.names(HM@assays[["RNA"]]@data) %in% IFN_genes),]
    #data=HM@assays[["RNA"]]@counts[which(row.names(HM@assays[["RNA"]]@counts) %in% IFN_genes),]
    data=as.data.frame(data)
    rnames=row.names(data)
    colnames(data)=gsub("-", ".", colnames(data))
    
    #subset to cell type
    cell_ind=row.names(labs)[which(labs$Cluster_num==j)]
    if(length(cell_ind)!=0){ #only move forward if cells exist
      data=data[,which(colnames(data) %in% cell_ind)]
      data=as.matrix(data)
      data2=matrix(as.numeric(data),    # Convert to numeric matrix
                   ncol = ncol(data))
      row.names(data2)=rnames
      colnames(data2)=colnames(data)
      temp=cbind(row.names(data2),as.data.frame(rowMeans(data2)))
      colnames(temp)=c("genes", paste0(i, "_", j))
      master_avgs=left_join(master_avgs, temp, by=c("genes"="genes"))
      annos=cbind(annos, rbind(samples_ord[i], cell_type_ord[j]))
      rm(data2)
    }
    rm(data)
  }
  rm(HM)
  rm(labs)
}

# Clean up missing values
master_avgs[is.na(master_avgs)] = 0
master_avgs=master_avgs[,-1]
row.names(master_avgs)=IFN_genes

# mycolors_master1=c("#9657A3", "#CFE1B6", "#B4D458", "#D7B063", "#AD82BA", "#88D1D2", "#D46483")
mycolors_master1=c("#262262", "#AD82BA", "#88D1D2",  "#B4D458", "#CFE1B6", "#D7B063")
names(mycolors_master1) <- c("0.Reference", samples_ord)

mycolors_master2=cluster_colors
names(mycolors_master2) <- cell_type_ord  

# Clean up annotations
annos=as.data.frame((t(annos)))
row.names(annos)=colnames(master_avgs)
colnames(annos)=c("Sample", "Cell")

myannocolors1=mycolors_master1[unique(annos$Sample)]
myannocolors2=mycolors_master2[unique(annos$Cell)]
myannocolors <- list(Sample = myannocolors1, Cell = myannocolors2)


annos2=arrange(annos, Cell, Sample) #arrange non-clustered heatmap by cell type
master_avgs2=master_avgs[,row.names(annos2)]


pheatmap(mat = log(master_avgs2+1),
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         filename = paste0("GSEA/IFNG_log_5_cell.pdf"))
pheatmap(mat = master_avgs2,
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         filename = paste0("GSEA/IFNG_5_cell.pdf"))

