library(data.table)
library(dplyr)
library(pheatmap)
library(randomcoloR)

######
sample_ten=TRUE
######

samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
samples_C=#a list of samples with different identifiers (for differing names in files)
cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", 
            "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
cluster_colors=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090") #0,1,2,3,4,5 etc 
cluster_colors2=c("#C0007C", "#2a4bd7", "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090",
                  "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff") #0,1,10,11 etc


annos=as.data.frame(matrix(nrow=2, ncol=0))
heatmaps=list()
for(i in 1:length(samples_A)){
  file_path=paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/", samples_A[i], "/input/cellHarmony/heatmaps/heatmap-all-cells-combined.txt")
  HM=as.data.frame(fread(file_path, header=T, sep="\t"))
  row.names(HM)=HM[,1]
  HM=HM[,-(1:2)]
  if(sample_ten==T){
    HM=HM[,sample(colnames(HM), ncol(HM)/10, replace=F)] #downsample for plot
  }
  if(i!=1){
    HM=HM[,-(grep(".Reference", colnames(HM)))]
  }
  heatmaps=c(list(HM), heatmaps)
  annos=cbind(annos, rbind(colnames(HM), HM[1,]))
  rm(HM)      
}
select_genes_common=Reduce(intersect, lapply(heatmaps, row.names))
select_genes_common=select_genes_common[-1]

for(i in 1:length(samples_A)){
  heatmaps[[i]]=heatmaps[[i]][select_genes_common,]
  if(i==1){
    masterHM=heatmaps[[i]]
  }else{
    masterHM=cbind(masterHM, heatmaps[[i]])
  }
}


#colors
mycolors_master1=c("#262262", "#9657A3", "#DCC2CF", "#CFE1B6", "#B4D458", "#D7B063", "#AD82BA", "#80C891", "#88D1D2", "#D46483", "#7FA7D3")
names(mycolors_master1) <- c("Reference", samples_B)

mycolors_master2=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                   "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090")
names(mycolors_master2) <- cell_type

#annotations
annos=t(annos)
annos[,1]=gsub("cellHarmony-", "",gsub(":.*", "", annos[,1]))
colnames(annos)=c("Sample", "Cell")
annos=as.data.frame(annos)
myannocolors1=mycolors_master1[unique(annos$Sample)]
myannocolors2=mycolors_master2[unique(annos$Cell)]
myannocolors <- list(Sample = myannocolors1, Cell = myannocolors2)

myColor=colorRampPalette(c("cyan", "black", "yellow"))(100)
masterHM2<- matrix(as.numeric(unlist(masterHM)),    # Convert to numeric matrix
                   ncol = ncol(masterHM))
colnames(masterHM2)=colnames(masterHM)
row.names(masterHM2)=row.names(masterHM)

myBreaks <- c(seq(min(masterHM2), 0, length.out=ceiling(100/2) + 1), 0,
              seq(max(masterHM2)/100, max(masterHM2), length.out=floor(100/2)))
myBreaks2 = myBreaks + seq_along(myBreaks) * .Machine$double.eps

#annos2=arrange(annos, Sample, Cell) #arrange non-clustered heatmap by sample
annos2=arrange(annos, Cell, Sample) #arrange non-clustered heatmap by cell type
masterHM22=masterHM2[,row.names(annos2)]
mygaps=NULL
for(i in unique(annos2$Cell)){
  mygaps=c(mygaps, which(annos2$Cell==i)[1])
}

pheatmap(mat = masterHM22,
         color=myColor,
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         fontsize_row = 5,
         width=20,
         gaps_col=mygaps,
         annotation_col = annos2,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = 12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/CH_heatmap_noclus_tiny_celltype_break_halfheight.pdf"))
pheatmap(mat = masterHM2,
         color=myColor,
         cluster_cols = T,
         cluster_rows = F,
         show_colnames = F,
         fontsize_row = 5,
         width=20,
         annotation_col = annos,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = 12,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/CH_heatmap_clus_tiny_halfheight.pdf"))
pheatmap(mat = masterHM2,
         color=myColor,
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = F,
         fontsize_row = 5,
         width=20,
         annotation_col = annos,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = 25,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/CH_heatmap_clusboth_tiny.pdf"))
pheatmap(mat = masterHM22,
         color=myColor,
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         fontsize_row = 5,
         width=20,
         gaps_col=mygaps,
         annotation_col = annos2,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = 25,
         filename = paste0("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4/CH_heatmap_clusrow_tiny_celltype_break.pdf"))

#IGKC MZB1 JCHAIN IRF8 TCF4 CCDC50 UGCG FCER1A ITM2C PLD4 SOX4 TCL1A IGLL1 TUBB C12orf75 DSTN PDIA6 NUCB2 PPP1R14B STMN1 PTGDS PLEK PLAC8
