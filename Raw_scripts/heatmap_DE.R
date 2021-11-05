# Load libraries
library(dplyr)
library(pheatmap)
library(data.table)
library(randomcoloR)
library(gplots)

# Set directory and parameters
setwd("~/Documents/Projects/Data_files_temp/Donor/CH_ccrem_5BM_0.4")
samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
controls=c("BM180221", "BM191119", "BM191227", "SMB1014", "SMB1064")

# Pull out Up-regulated, Down-regulated, and  All regulated genes by cell type across all samples (union)
keep_genes=NULL
for(j in 1:length(cell_type)){
  temp_Up=NULL
  temp_Down=NULL
  for(i in 1:length(samples_A)){
    file_path=paste0(samples_A[i], "/input/cellHarmony/DifferentialExpression_Fold_1.5_adjp_0.05/GE.", cell_type[j], "_", samples_B[i], "_vs_", cell_type[j], "_cellHarmony-Reference.txt")
    try({
      file_1=read.table(file_path, sep="\t", header=T)
      temp_Up=c(temp_Up, file_1[which(file_1$LogFold >=1), "Symbol"])
      temp_Down=c(temp_Down, file_1[which(file_1$LogFold <=1), "Symbol"])
    }, silent=TRUE)
  }
  # assign(paste0(cell_type[j], "_Up"), unique(temp_Up))
  # assign(paste0(cell_type[j], "_Down"), unique(temp_Down))
  assign(paste0(cell_type[j], "_All"), union(temp_Up, temp_Down))
  keep_genes=c(keep_genes, union(temp_Up, temp_Down))
}

# Read in BPDCN expression files and keep only All regulated genes from above (to cut down on size)
# Read in cellHarmony classifications for all BPDCN samples
groups=data.frame(barcode="trash", cluster="trash", donor="trash")
for(i in 1:length(samples_A)){
  temp_exp=as.data.frame(fread(paste0(samples_A[i], "/input/exp.", samples_B[i], ".txt")))
  assign(samples_B[i], temp_exp[which(temp_exp[,1] %in% keep_genes),])
  rm(temp_exp)
  temp_groups=read.table(paste0(samples_A[i], "/input/cellHarmony/QueryGroups.cellHarmony.txt"))
  temp_groups[,3]=rep(samples_B[i], nrow(temp_groups))
  colnames(temp_groups)=c("barcode", "cluster", "donor")
  groups=rbind(groups, temp_groups)
  rm(temp_groups)
}

# Read in control expression file (cellHarmony input file) and keep only All regulated genes from above
cont_exp=as.data.frame(fread(paste0(samples_A[i], "/reference/ExpressionInput/exp.Controls.txt")))
cont_exp=cont_exp[which(cont_exp[,1] %in% keep_genes),]
cont_exp=cont_exp[,-2]
cont_groups=as.data.frame(fread(paste0(samples_A[i], "/reference/ExpressionInput/FinalGroups.txt")))
cont_groups[,2]=cont_groups[,3] #shift around columns
cont_groups[,3]=cont_groups[,1]
cont_groups[,3]=gsub("[^_]*_", "", cont_groups[,3]) #pull out donor number
map = setNames(controls, unique(cont_groups[,3])) # has the potential to be dangerous!
cont_groups[,3] <- map[unlist(cont_groups[,3])]
colnames(cont_groups)=c("barcode", "cluster", "donor")

# Seperate out control exp by donor
sum(colnames(cont_exp)[2:ncol(cont_exp)]!=cont_groups[,1]) #0 <- safety check
for(k in 1:length(controls)){
  inds=which(cont_groups[,3]==controls[k])
  temp_exp=cont_exp[,c(1, inds+1)]
  assign(controls[k], temp_exp)
  rm(temp_exp)
}
groups=rbind(groups, cont_groups) #combine groups files
samples_A=c(controls, samples_A)
samples_B=c(controls, samples_B)

# Make cell type matrices of All regulated genes by cell type
for(j in 1:length(cell_type)){
  print(cell_type[j])
  message("Processing Cell Type: ", cell_type[j])
  temp_genes=get(paste0(cell_type[j], "_All")) #regulated genes for this cell type
  if(!is.null(temp_genes)){
    temp_df=data.frame(genes=temp_genes) #blank data frame for the results
    donor=NULL #blank list for keeping track of donor rows
    for(i in 1:length(samples_B)){
      print(samples_B[i])
      cells=groups$barcode[intersect(which(groups$donor==samples_B[i]), which(groups$cluster==cell_type[j]))] #cells that match the donor AND cell type
      temp_df_samp=get(samples_B[i])[,c("UID", cells)]
      if(!is.null(ncol(temp_df_samp))){ #if there is anything returned (aka, cells classified to that type in the sample)
        temp_df=left_join(temp_df, temp_df_samp, by=c("genes"="UID"))
        donor=c(donor, rep(samples_B[i], length(cells)))
      }
    }
    temp_df=rbind(c("Donor", donor), temp_df)
    row.names(temp_df)=temp_df[,1]
    temp_df=temp_df[,-1]
    temp_df[is.na(temp_df)] = 0 # NA's of all zero values from join, change to zeros
    assign(paste0(cell_type[j], "_mtx"), temp_df)
    rm(temp_df, temp_genes, temp_df_samp, donor)
  }else{
    message("Skipping ", cell_type[j], "! No cells classified as that type")
  }
}

# Generate heatmaps by cell type

mycolors_master=distinctColorPalette(length(samples_A)) #set Donor colors
names(mycolors_master) <- samples_B

###########
# # No donor clustering
# dir.create("Heatmaps_noclus")
# for(j in 1:length(cell_type)){
#   try({
#     mini=get(paste0(cell_type[j], "_mtx"))
#     mini=mini[-1,]
#     names=row.names(mini)
#     mini=as.data.frame(sapply(mini, as.numeric))
#     row.names(mini)=names
#     anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
#     
#     mycolors=mycolors_master[unique(anno_cols$Donor)]
#     mycolors <- list(Donor = mycolors)
#     
#     pdf(paste0("Heatmaps_noclus/", cell_type[j], "_heatmap.pdf"), width = 20, height = 10)
#     par(mar = c(8,4,8,12), xpd = T)
#     pheatmap(mat = mini,
#              color=colorRampPalette(c("black", "yellow"))(100),
#              cluster_cols = F,
#              show_colnames = F,
#              annotation_col = anno_cols,
#              annotation_colors = mycolors)
#     dev.off()
#   })
# }
# 
# # With donor clustering
# dir.create("Heatmaps_clus")
# for(j in 1:length(cell_type)){
#   try({
#     mini=get(paste0(cell_type[j], "_mtx"))
#     mini=mini[-1,]
#     names=row.names(mini)
#     mini=as.data.frame(sapply(mini, as.numeric))
#     row.names(mini)=names
#     anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
#     
#     mycolors=mycolors_master[unique(anno_cols$Donor)]
#     mycolors <- list(Donor = mycolors)
#     
#     pdf(paste0("Heatmaps_clus/", cell_type[j], "_heatmap_cellclus.pdf"), width = 20, height = 10)
#     par(mar = c(8,4,8,12), xpd = T)
#     pheatmap(mat = mini,
#             color=colorRampPalette(c("black", "yellow"))(100),
#             cluster_cols = T,
#             show_colnames = F,
#             annotation_col = anno_cols,
#             annotation_colors = mycolors)
#     dev.off()
#   })
# }
# 
# 
# # Exclude RPS and RPL genes (with donor clustering)
# dir.create("Heatmaps_clus_noribo")
# for(j in 1:length(cell_type)){
#   try({
#     mini=get(paste0(cell_type[j], "_mtx"))
#     mini=mini[-1,]
#     names=row.names(mini)
#     mini=as.data.frame(sapply(mini, as.numeric))
#     row.names(mini)=names
#     
#     rp.ch <- rownames(mini)[grepl("^RPS|^RPL", rownames(mini))]
#     mini <- mini[setdiff(rownames(mini), rp.ch),]
#     
#     anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
#     
#     mycolors=mycolors_master[unique(anno_cols$Donor)]
#     mycolors <- list(Donor = mycolors)
#     
#     pdf(paste0("Heatmaps_clus_noribo/", cell_type[j], "_heatmap_cellclus_noribo.pdf"), width = 20, height = 10)
#     par(mar = c(8,4,8,12), xpd = T)
#     pheatmap(mat = mini,
#              color=colorRampPalette(c("black", "yellow"))(100),
#              cluster_cols = T,
#              show_colnames = F,
#              annotation_col = anno_cols,
#              annotation_colors = mycolors)
#     dev.off()
#   })
# }
# 
# # Reshape (with donor clustering and no ribo genes) at 14 genes per inch plus 1 inch for header
# dir.create("Heatmaps_clus_noribo_reshape")
# for(j in 1:length(cell_type)){
#   try({
#     mini=get(paste0(cell_type[j], "_mtx"))
#     mini=mini[-1,]
#     names=row.names(mini)
#     mini=as.data.frame(sapply(mini, as.numeric))
#     row.names(mini)=names
#     
#     rp.ch <- rownames(mini)[grepl("^RPS|^RPL", rownames(mini))]
#     mini <- mini[setdiff(rownames(mini), rp.ch),]
#     
#     anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
#     
#     mycolors=mycolors_master[unique(anno_cols$Donor)]
#     mycolors <- list(Donor = mycolors)
#     
#     pheatmap(mat = mini,
#              color=colorRampPalette(c("black", "yellow"))(100),
#              cluster_cols = T,
#              show_colnames = F,
#              annotation_col = anno_cols,
#              annotation_colors = mycolors,
#              fontsize_row = 5,
#              height = (nrow(mini)/14)+1,
#              filename = paste0("Heatmaps_clus_noribo_reshape/", cell_type[j], "_heatmap_cellclus_noribo_shape.pdf"))
#   })
# }
########### OLD PLOTS


# Reshape (with donor clustering and no ribo genes) at 14 genes per inch plus 1 inch for header PLUS CONTROLS
dir.create("HeatmapsFC_clus_noribo_reshape_pluscontrols")
for(j in 1:length(cell_type)){
  try({
    mini=get(paste0(cell_type[j], "_mtx"))
    mini=mini[-1,]
    names=row.names(mini)
    mini=as.data.frame(sapply(mini, as.numeric))
    row.names(mini)=names
    
    rp.ch <- rownames(mini)[grepl("^RPS|^RPL", rownames(mini))]
    mini <- mini[setdiff(rownames(mini), rp.ch),]
    
    anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
    
    mycolors=mycolors_master[unique(anno_cols$Donor)]
    mycolors <- list(Donor = mycolors)
    
    pheatmap(mat = mini,
             color=colorRampPalette(c("black", "yellow"))(100),
             cluster_cols = T,
             show_colnames = F,
             annotation_col = anno_cols,
             annotation_colors = mycolors,
             fontsize_row = 5,
             height = (nrow(mini)/14)+1,
             filename = paste0("HeatmapsFC_clus_noribo_reshape_pluscontrols/", cell_type[j], "_heatmapFC_cellclus_noribo_shape.pdf"))
  })
}

# Reshape (with NO donor clustering and no ribo genes) at 14 genes per inch plus 1 inch for header PLUS CONTROLS
dir.create("HeatmapsFC_noclus_noribo_reshape_pluscontrols")
for(j in 1:length(cell_type)){
  try({
    mini=get(paste0(cell_type[j], "_mtx"))
    mini=mini[-1,]
    names=row.names(mini)
    mini=as.data.frame(sapply(mini, as.numeric))
    row.names(mini)=names
    
    rp.ch <- rownames(mini)[grepl("^RPS|^RPL", rownames(mini))]
    mini <- mini[setdiff(rownames(mini), rp.ch),]
    
    anno_cols=as.data.frame(t(get(paste0(cell_type[j], "_mtx"))[1,]))
    
    mycolors=mycolors_master[unique(anno_cols$Donor)]
    mycolors <- list(Donor = mycolors)
    
    pheatmap(mat = mini,
             color=colorRampPalette(c("black", "yellow"))(100),
             cluster_cols = F,
             show_colnames = F,
             annotation_col = anno_cols,
             annotation_colors = mycolors,
             fontsize_row = 5,
             height = (nrow(mini)/14)+1,
             filename = paste0("HeatmapsFC_noclus_noribo_reshape_pluscontrols/", cell_type[j], "_heatmapFC_noellclus_noribo_shape.pdf"))
  })
}





