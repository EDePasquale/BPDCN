


library(data.table)
library(dplyr)
library(pheatmap)
library(randomcoloR)
library(matrixStats)
#install.packages("BSDA")
library(BSDA)

geneset="TNFA"

IFN_genes=read.table(paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/geneset_",geneset,".txt"), sep="\t") #IFN alpha
IFN_genes=c(IFN_genes[-(1:2),])

reference=readRDS("~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem/TCELL_SUB/Seurat_int_0.5_5BM.rds")

samples_A=#a list of samples\
samples_B=#a list of samples with different identifiers (for differing names in files)
samples_ord=#a list of samples with different identifiers (for differing names in files)
samples_no=c(1,2,3,4,5)


cell_type=c("CD4NaiveT", "CD8MemoryT",  "NK",  "CD4MemoryT", "CD8NaiveT")
cell_type_ord=c("B.CD4NaiveT",  "D.CD8MemoryT",  "E.NK",  "A.CD4MemoryT", "C.CD8NaiveT")
cluster_colors=c("#C0007C", "#e9debb",  "#BE8A66",  "#ffcdf3", "#708090") #0,1,2,3,4,5 etc
cluster_no=c(1,4,6,9,17)

annos=as.data.frame(matrix(nrow=2, ncol=0))
master_avgs=as.data.frame(matrix(nrow=length(IFN_genes), ncol=1))
master_avgs[,1]=IFN_genes
colnames(master_avgs)="genes"

master_stds=as.data.frame(matrix(nrow=length(IFN_genes), ncol=1))
master_stds[,1]=IFN_genes
colnames(master_stds)="genes"

master_n=as.data.frame(matrix(nrow=length(IFN_genes), ncol=1))
master_n[,1]=IFN_genes
colnames(master_n)="genes"

#########
# REF
#########

labs3=reference@active.ident
labs2=as.numeric(as.character(labs3))
names(labs2)=names(labs3)
labs2=labs2+1
temp1=labs2
map = setNames(cell_type, c(1,2,3,4,5))
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
  cell_ind=row.names(labs)[which(labs$Cluster_num==cluster_no[j])]
  if(length(cell_ind)!=0){ #only move forward if cells exist
    data=data[,which(colnames(data) %in% cell_ind)]
    data=as.matrix(data)
    data2=matrix(as.numeric(data),    # Convert to numeric matrix
                 ncol = ncol(data))
    row.names(data2)=rnames
    colnames(data2)=colnames(data)
    temp=cbind(row.names(data2),as.data.frame(rowMeans(data2)))
    colnames(temp)=c("genes", paste0("R_", cluster_no[j]))
    master_avgs=left_join(master_avgs, temp, by=c("genes"="genes"))
    temp2=cbind(row.names(data2),as.data.frame(rowSds(data2)))
    colnames(temp2)=c("genes", paste0("R_", cluster_no[j]))
    master_stds=left_join(master_stds, temp2, by=c("genes"="genes"))
    temp3=cbind(row.names(data2),as.data.frame(ncol(data2)))
    colnames(temp3)=c("genes", paste0("R_", cluster_no[j]))
    master_n=left_join(master_n, temp3, by=c("genes"="genes"))
    annos=cbind(annos, rbind("0.Reference", cell_type_ord[j]))
    #assign(paste0("R_", cell_type[j]), master_avgs)
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
    cell_ind=row.names(labs)[which(labs$Cluster_num==cluster_no[j])]
    if(length(cell_ind)!=0){ #only move forward if cells exist
      data=data[,which(colnames(data) %in% cell_ind)]
      data=as.matrix(data)
      data2=matrix(as.numeric(data),    # Convert to numeric matrix
                   ncol = ncol(data))
      row.names(data2)=rnames
      colnames(data2)=colnames(data)
      temp=cbind(row.names(data2),as.data.frame(rowMeans(data2)))
      colnames(temp)=c("genes", paste0(i, "_", cluster_no[j]))
      master_avgs=left_join(master_avgs, temp, by=c("genes"="genes"))
      temp2=cbind(row.names(data2),as.data.frame(rowSds(data2)))
      colnames(temp2)=c("genes", paste0(i, "_", cluster_no[j]))
      master_stds=left_join(master_stds, temp2, by=c("genes"="genes"))
      temp3=cbind(row.names(data2),as.data.frame(ncol(data2)))
      colnames(temp3)=c("genes", paste0(i, "_", cluster_no[j]))
      master_n=left_join(master_n, temp3, by=c("genes"="genes"))
      annos=cbind(annos, rbind(samples_ord[i], cell_type_ord[j]))
      #assign(paste0(samples_B[i],"_", cell_type[j]), master_avgs)
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

master_stds[is.na(master_stds)] = 0
master_stds=master_stds[,-1]
row.names(master_stds)=IFN_genes

master_n[is.na(master_n)] = 0
master_n=master_n[,-1]
row.names(master_n)=IFN_genes

##########
# Significance
##########

#cluster_no=c(1,4,6,9,17)
res=as.data.frame(matrix(nrow=length(IFN_genes), ncol=1))
res[,1]=IFN_genes
colnames(res)="genes"

for(j in cluster_no){

  inds=grep(paste0("_", j, "$"), colnames(master_avgs))
  avgs_temp=master_avgs[,inds]
  stds_temp=master_stds[,inds]
  n_temp=master_n[,inds]
  for(k in 2:length(inds)){
    vec=NULL
    for(m in 1:nrow(avgs_temp)){
      x=tsum.test(mean.x = avgs_temp[m,1],
                  s.x = stds_temp[m,1],
                  n.x = n_temp[m,1],
                  mean.y = avgs_temp[m,k],
                  s.y = stds_temp[m,k],
                  n.y = n_temp[m,k])
      vec=c(vec, x$p.value)
    }
    y=as.data.frame(vec)
    colnames(y)=colnames(avgs_temp)[k]
    res=cbind(res, y)
    
  }
  
}
res[is.na(res)] = 1
row.names(res)=res[,1]
res=res[,-1]
res2=data.matrix(res)
res3=matrix(p.adjust(res2, method="bonferroni"),    # Convert to numeric matrix
             ncol = ncol(res2))
mins=rowMins(as.matrix(res3))
row_inds=which(mins<=0.05)

##########







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


#annos2=arrange(annos, Cell, Sample) #arrange non-clustered heatmap by cell type
annos2=arrange(annos, Sample, Cell) #arrange non-clustered heatmap by sample
master_avgs2=master_avgs[,row.names(annos2)]


pheatmap(mat = log(master_avgs2+1),
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_log_5_sample_TNK.pdf"))
pheatmap(mat = master_avgs2,
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_5_sample_TNK.pdf"))

#sig
pheatmap(mat = log(master_avgs2[row_inds,]+1),
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_log_5_sample_TNK_sig.pdf"))
pheatmap(mat = master_avgs2[row_inds,],
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_5_sample_TNK_sig.pdf"))

#sig and min
thresh=0.5 #random
master_avgs3=data.matrix(master_avgs2)
master_avgs3=matrix(as.numeric(master_avgs3),    # Convert to numeric matrix
            ncol = ncol(master_avgs3))
maxs=rowMaxs(master_avgs3)
thresh_inds=which(maxs>=thresh)
row_inds2=intersect(row_inds, thresh_inds)
pheatmap(mat = log(master_avgs2[row_inds2,]+1),
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_log_5_sample_TNK_sig_0.25min.pdf"))
pheatmap(mat = master_avgs2[row_inds2,],
         show_colnames = T,
         fontsize_row = 5,
         annotation_colors = myannocolors,
         cluster_cols = F,
         annotation_col = annos2,
         #scale="row",
         width=10,
         height=12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/GSEA/",geneset,"_5_sample_TNK_sig_0.25min.pdf"))

