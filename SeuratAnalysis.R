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

#Read in the normalized Seurat objects
UCC_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_norm_seur.rds")
URN_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/URN_norm_seur.rds")

#REad in the integrated normalized Seurat Object
UCC_int_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/URN_int_seur.rds")

#Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
#it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
Seur_target<-UCC_int_seur
Target_name<-"UCC_Integrated"


#Making DE lists and volcano plots for the indicated Seurat object
WT.PIC.de<-DEVolcano(Seur_target,"ZSCCP","ZSCCC","Sample")
WT.LPS.de<-DEVolcano(Seur_target,"ZSCCL","ZSCCC","Sample")

J.PIC.de<-DEVolcano(Seur_target,"ZSCJP","ZSCJC","Sample")
J.LPS.de<-DEVolcano(Seur_target,"ZSCJL","ZSCJC","Sample")

Wt.pic.sig<-WT.PIC.de[WT.PIC.de$expression!="NS",c(2,5,7)]
Wt.lps.sig<-WT.LPS.de[WT.LPS.de$expression!="NS",c(2,5,7)]
J.pic.sig<-J.PIC.de[J.PIC.de$expression!="NS",c(2,5,7)]
J.lps.sig<-J.LPS.de[J.LPS.de$expression!="NS",c(2,5,7)]

vennlist<-list(
  WTP=Wt.pic.sig$gene,
  WTL=Wt.lps.sig$gene,
  JP=J.pic.sig$gene,
  JL=J.lps.sig$gene
)

#Venn diagram of all significant genes and all non Zika samples. 
ggvenn(
  vennlist, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title=paste0(Target_name,"All Significant Genes vs their isogenic control"))

ggvenn(
  vennlist[1:2], 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title=paste0(Target_name,"All Significant Genes in WT samples"))

ggvenn(
  vennlist[3:4], 
  fill_color = c( "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title=paste0(Target_name,"All Significant Genes in J samples"))

ggvenn(
  vennlist[c(1,3)], 
  fill_color = c( "#868686FF", "#0073C2FF"),
  stroke_size = 0.5, set_name_size = 4
) +
  labs(title=paste0(Target_name,": All Significant Genes in P samples"))


# PCA and visualization
Seur_target <- RunPCA(Seur_target, verbose = FALSE)
VizDimLoadings(Seur_target, dims = 1:2, reduction = "pca")
Idents(Seur_target) <- "Treatment"

DimPlot(Seur_target, reduction = "pca") 

#Clustering and UMAP visualization
Seur_target <- FindNeighbors(Seur_target, dims = 1:10, verbose = FALSE)
Seur_target <- FindClusters(Seur_target, verbose = FALSE)


Seur_target <- RunUMAP(Seur_target, dims = 1:10, verbose = FALSE)

DimPlot(Seur_target, reduction = "umap",group.by="Treatment",label=FALSE)+
labs(title=paste0(Target_name,": UMAP"))

saveRDS(Seur_target,paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"_final.rds"))


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
  
  print(volcano)
  invisible(volcano)
  return(deFile)
  
  
}

