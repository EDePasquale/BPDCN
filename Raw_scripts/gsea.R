########
#      #
# GSEA #
#      #
########

# Load packages
library(dplyr)
library(gage)
library(fgsea)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(randomcoloR)

# Read in GO_file and set pval parameters
setwd("~/Documents/Projects/Data_files_temp/Donor")
#GO_file="GSEA/h.all.v7.4.symbols.gmt"
GO_file="GSEA/c2.all.v7.4.symbols.gmt"
#GO_file="GSEA/msigdb.v7.4.symbols.gmt"
pval=0.05

# Read in log fold change values for the gene_list parameter (from cellHarmony output)

samples_A=#a list of samples
samples_B=#a list of samples with different identifiers (for differing names in files)
samples_ord=#a list of samples with different identifiers (for differing names in files)
cell_type=c("CD4NaiveT", "B", "MidEryth", "CD8MemoryT", "Mono", "NK", "LateEryth", "EarlyEryth", "CD4MemoryT", "PreB", "ProB", "HSCProg", "cDC", "ncMono", "Plasma", "pDC", "CD8NaiveT")
cell_type_ord=c("N.CD4NaiveT", "K.B", "C.MidEryth", "P.CD8MemoryT", "E.Mono", "Q.NK", "D.LateEryth", "B.EarlyEryth", "M.CD4MemoryT", "J.PreB", "I.ProB", "A.HSCProg", "G.cDC", "F.ncMono", "L.Plasma", "H.pDC", "O.CD8NaiveT")
# GSEA function
GSEA <- function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  if(dim(fgRes)[1]!=0){
    
    ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
    gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
    
    ups = as.data.frame(gaRes$greater) %>% 
      tibble::rownames_to_column("Pathway") %>% 
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")
    
    downs = as.data.frame(gaRes$less) %>% 
      tibble::rownames_to_column("Pathway") %>% 
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")

    if(dim(rbind(ups,downs))[1]!=0){
      ## Define up / down pathways which are significant in both tests
      keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
      keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
      
      fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                                  c( keepups$pathway, keepdowns$pathway))), ] %>% 
        arrange(desc(NES))
      fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
      
      fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
      filtRes = rbind(head(fgRes, n = 10),
                      tail(fgRes, n = 10 ))
      
      
      upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
      downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
      colos = c(upcols, downcols)
      names(colos) = 1:length(colos)
      filtRes$Index = as.factor(1:nrow(filtRes))
      
      g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
        geom_col( aes(fill = Index )) +
        scale_fill_manual(values = colos ) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="GSEA - Biological Processes") + 
        theme_minimal()
      
      
      output = list("Results" = fgRes, "Plot" = g)
      return(output)
    }else{
      return("Results"=0)
    }
  }else{
    return("Results"=0)
  }
}

###### 
#Run GSEA
# results_df=data.frame(Cell="trash", Sample="trash", Pathway="trash", NES="trash")
# for(j in 1:length(cell_type)){
#   for(i in 1:length(samples_A)){
#     file_path=paste0("CH_ccrem_5BM_0.4/", samples_A[i], "/input/cellHarmony/DifferentialExpression_Fold_1.5_adjp_0.05/GE.", cell_type[j], "_", samples_B[i], "_vs_", cell_type[j], "_cellHarmony-Reference.txt")
#     try({
#       file_1=read.table(file_path, sep="\t", header=T)
#       temp1=file_1$LogFold
#       names(temp1)=file_1$Symbol
#       temp1=sort(temp1, decreasing = T)
#       temp2=GSEA(temp1, GO_file, pval)
#       if(class(temp2$Results)=="data.frame"){
#         assign(paste0(samples_B[i], "_", cell_type[j], "_results"), temp2)
#         x=cbind(cell_type[j], samples_B[i], temp2$Results[,1], temp2$Results[,6])
#         colnames(x)=c("Cell", "Sample", "Pathway", "NES")
#         results_df=rbind(results_df, x)
#         print(paste0(cell_type[j], " (", samples_B[i], "):    ", temp2$Results[,1], ",    NES = ", temp2$Results[,6]))
#       }
#     }, silent=TRUE)
#   }
# } old GSEA code for only significant genes
######

# Run GSEA
results_df=data.frame(Cell="trash", Sample="trash", Pathway="trash", NES="trash")
for(j in 1:length(cell_type)){
  for(i in 1:length(samples_A)){
    file_path=paste0("CH_ccrem_5BM_0.4/", samples_A[i], "/input/cellHarmony/OtherFiles/exp.MarkerFinder-cellHarmony-reference__", samples_B[i], "-AllCells-folds.txt")
    try({
      file_1=read.table(file_path, sep="\t", header=T)
      colnames(file_1)=gsub("_.*", "", colnames(file_1))
      temp1=file_1[,which(colnames(file_1)==cell_type[j])]
      names(temp1)=file_1[,1]
      temp1=sort(temp1, decreasing = T)
      temp2=GSEA(temp1, GO_file, pval)
      if(class(temp2$Results)=="data.frame"){
        assign(paste0(samples_ord[i], "_", cell_type_ord[j], "_results"), temp2)
        x=cbind(cell_type_ord[j], samples_ord[i], temp2$Results[,1], temp2$Results[,6])
        colnames(x)=c("Cell", "Sample", "Pathway", "NES")
        results_df=rbind(results_df, x)
        #print(paste0(cell_type[j], " (", samples_B[i], "):    ", temp2$Results[,1], ",    NES = ", temp2$Results[,6]))
      }
      rm(temp1)
      rm(temp2)
    }, silent=TRUE)
  }
}


# mycolors_master1=c("#9657A3", "#CFE1B6", "#B4D458", "#D7B063", "#AD82BA", "#88D1D2", "#D46483")
mycolors_master1=c("#CFE1B6", "#B4D458", "#D7B063", "#AD82BA", "#88D1D2")
names(mycolors_master1) <- samples_ord

mycolors_master2=c("#C0007C", "#2a4bd7", "#f08080", "#e9debb", "#ffee33", "#BE8A66", "#ad2323", "#ff9233", "#ffcdf3", "#9dafff",
                 "#8126c0", "#000000", "#1d6914", "#29d0d0", "#575757",  "#81c57a", "#708090")
names(mycolors_master2) <- cell_type_ord



# Prep for heatmap
results_df2=results_df
results_df2=results_df2[-1,]
#anno_cols=results_df2[,1:2]


results_df2[,1]=paste0(results_df2[,1], "_", results_df2[,2])
results_df2=results_df2[,-2]

data_wide <- spread(results_df2, key = "Pathway", value = "NES", fill = 0)
anno_cols=cbind(data_wide[,1], data_wide[,1])
row.names(anno_cols)=anno_cols[,1]
anno_cols[,1]=gsub("_.*", "", anno_cols[,1])
anno_cols[,2]=gsub(".*_", "", anno_cols[,2])
colnames(anno_cols)=c("Cell", "Sample")
anno_cols=as.data.frame(anno_cols)

row.names(data_wide)=data_wide[,1]
data_wide=data_wide[,-1]
data_wide2=as.matrix(data_wide)
data_wide3<- matrix(as.numeric(data_wide2),    # Convert to numeric matrix
          ncol = ncol(data_wide2))
row.names(data_wide3)=row.names(data_wide2)
colnames(data_wide3)=colnames(data_wide2)
data_wide3=rbind(data_wide3, colSums(abs(data_wide3)))
data_wide4=data_wide3[,order(-data_wide3[nrow(data_wide3),])]



# Make heatmap
myannocolors1=mycolors_master1[unique(anno_cols$Sample)]
myannocolors2=mycolors_master2[unique(anno_cols$Cell)]
myannocolors <- list(Sample = myannocolors1, Cell = myannocolors2)

myColor=colorRampPalette(c("blue", "white", "red"))(100)
vals=data_wide4[-(nrow(data_wide4)),]
myBreaks <- c(seq(min(vals), 0, length.out=ceiling(100/2) + 1), 0,
              seq(max(vals)/100, max(vals), length.out=floor(100/2)))

myBreaks2 = myBreaks + seq_along(myBreaks) * .Machine$double.eps

pheatmap(mat = t(data_wide4[-(nrow(data_wide4)),]),
         color=myColor,
         cluster_cols = F,
         show_colnames = T,
         fontsize_row = 5,
         width=11,
         annotation_col = anno_cols,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = (nrow(t(data_wide4[-(nrow(data_wide4)),]))/10)+3,
         filename = paste0("GSEA/GSEA_heatmap_C2_noclus_5.pdf"))
pheatmap(mat = t(data_wide4[-(nrow(data_wide4)),]),
         color=myColor,
         cluster_cols = T,
         show_colnames = T,
         fontsize_row = 5,
         width=11,
         annotation_col = anno_cols,
         annotation_colors = myannocolors,
         breaks=myBreaks2,
         height = (nrow(t(data_wide4[-(nrow(data_wide4)),]))/10)+4,
         filename = paste0("GSEA/GSEA_heatmap_C2_clus_5.pdf"))


