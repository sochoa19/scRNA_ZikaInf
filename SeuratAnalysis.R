####Script for analysis of combined Seurat Objects
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Matrix)
library(stringr)
library(future)
library(ggvenn)
library(purrr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

#Read in the normalized Seurat objects
UCC_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_norm_seur.rds")
URN_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/URN_norm_seur.rds")

#Read in the integrated normalized Seurat Object
UCC_int_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_int_seur.rds")

#Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
#it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
Seur_target<-UCC_int_seur
Target_name<-"UCC_Integrated"

#########Make Volcano plots for all Ctrl vs LPs or Ctrl vs PIC isnide of each genetic treatment
#########Also finds the shared significant genes across all 
de.L.list<-list()
de.P.list<-list()
sig.L.names<-list()
sig.P.names<-list()
hotfix<-c("4","7","C","J")

#switched out this : unique(Seur_target$Background) for hotfix to rerun starting from 4

for (i in unique(Seur_target$Background)){
  Cname=paste0("ZSC",i,"C")
  Lname=paste0("ZSC",i,"L")
  Pname=paste0("ZSC",i,"P")
  de.L<-DEVolcano(Seur_target,Lname,Cname,"Sample")
  de.P<-DEVolcano(Seur_target,Pname,Cname,"Sample")
  sig.L.names<-append(sig.L.names,de.L[de.L$expression!="NS",7])
  sig.P.names<-append(sig.P.names,de.P[de.P$expression!="NS",7])
  de.L.list[[i]]<-de.L
  de.P.list[[i]]<-de.P
}
#makes a unqiue lsit of all significant genes across al backgrounds
all.sig.L<-unique(sig.L.names)
all.sig.P<-unique(sig.P.names)
#saves a list fo dataframes, each one with the DE results for each treatment vs its isogenic ctrl
saveRDS(de.L.list,paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/vsL_DE.rds"))
saveRDS(de.P.list,paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/vsP_DE.rds"))

#####Make a list of dataframe of only the upregulated or downregulated genes across samples
##and separate a list of genes 

L.up.gene_sets <- map(de.L.list, ~{
  .x %>%
    filter(expression != "Up") %>%
    pull(gene) %>%
    unique()
})
L.up.WT.genes<-L.up.gene_sets[["C"]]

L.down.gene_sets <- map(de.L.list, ~{
  .x %>%
    filter(expression != "Down") %>%
    pull(gene) %>%
    unique()
})
L.down.WT.genes<-L.down.gene_sets[["C"]]

P.up.gene_sets <- map(de.P.list, ~{
  .x %>%
    filter(expression != "Up") %>%
    pull(gene) %>%
    unique()
})
P.up.WT.genes<-P.up.gene_sets[["C"]]

P.down.gene_sets <- map(de.P.list, ~{
  .x %>%
    filter(expression != "Down") %>%
    pull(gene) %>%
    unique()
})
P.down.WT.genes<-P.down.gene_sets[["C"]]

#Calculate the jaccard index between the control DE significant genes
#and the ones found in other backgrounds 
setname<-P.down.gene_sets
compname<-P.down.WT.genes

jaccard <- map_dbl(setname, function(x){
  
  intersection <- length(intersect(compname, x))
  union <- length(union(compname, x))
  
  intersection / union
  
})

#Turning it into a matrix
J <- matrix(
  jaccard,
  ncol = 1,
  dimnames = list(names(jaccard), "WT")
)

#Visualize jaccard in a heatmap
pheatmap(
  J,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("white", "gold", "red"))(100)
)




# PCA and visualization
Seur_target <- RunPCA(Seur_target, verbose = FALSE)
VizDimLoadings(Seur_target, dims = 1:2, reduction = "pca")
Idents(Seur_target) <- "Treatment"

DimPlot(Seur_target, reduction = "pca") 

#Clustering and UMAP visualization
Seur_target <- FindNeighbors(Seur_target, dims = 1:10, verbose = FALSE)
Seur_target <- FindClusters(Seur_target, verbose = FALSE)


Seur_target <- RunUMAP(Seur_target, dims = 1:10, verbose = FALSE)

#Plotting UMAP grouping by treatments
DimPlot(Seur_target, reduction = "umap",group.by="Treatment",label=FALSE)+
labs(title=paste0(Target_name,": UMAP"))
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Figures/UMAP/TreatmentUmap.png"))

#Plotting UMAP grouping by sample
DimPlot(Seur_target, reduction = "umap",group.by="Sample",label=FALSE)+
  labs(title=paste0(Target_name,": UMAP"))
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Figures/UMAP/SampleUmap.png"))

saveRDS(Seur_target,paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/",Target_name,"_final.rds"))


##FUNCTIONS LIVE BELOW!!!


#DE of WT background Ctrl vs LPS & Ctrl vs PIC
#Focus is the level at which you want to compare expression, by sample, by background etc,
#a string of the Seurat metadata column
#Target1 and Target2 are strings of objects inside of the focus column
#Target1 is the group that you wish to analyse, Target2 is the baseline for comparison
#SeurFile is the Seurat object to be used for DE
#OUtput will be the complete DE matrix for your selected groups and a volcano plot labelling 
#the 10 genes with the lowest pvalue with a fold change above 2
DEVolcano<-function(SeurFile,Target1,Target2,Focus){
  Idents(SeurFile) <- Focus
  
  
  deFile <- FindMarkers(SeurFile, ident.1 = Target1, ident.2 = Target2, verbose = FALSE)
  
  #Volcano Plots for easy vizualization of DE
  deFile <- deFile %>%
    mutate(expression = case_when(
      avg_log2FC >= 1 & p_val_adj <= 0.05 ~ "Up",
      avg_log2FC <= -1 & p_val_adj <= 0.05 ~ "Down",
      TRUE ~ "NS"
    ))
  
  deFile$gene<-rownames(deFile)
  volcano<-ggplot(deFile, aes(x = avg_log2FC, y = -log10(p_val_adj), color = expression)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red")) +
    
    # Threshold lines
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    
    #Labelling the top 10 genes
    geom_text_repel(
      data = head(deFile[deFile$expression != "NS", ][order(deFile[deFile$expression != "NS", ]$p_val_adj), ], 10),
      aes(label = gene),
      color = "black",
      size = 3
    ) +
    
    # Formatting
    labs(
      title = paste0(Target_name," Volcano Plot:",Target1," VS ",Target2),
      x = "Log2 Fold Change",
      y = "-Log10 p-value-adj"
    ) +
    theme_minimal()
  ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Figures/VolcanoPlots/",Target1 ,".png"), width = 6, height = 4, dpi = 300)
  print(volcano)
  invisible(volcano)
  return(deFile)
  
  
}

