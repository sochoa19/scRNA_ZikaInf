####Script for analysis of combined Seurat Objects
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Matrix)
library(readr)
library(stringr)
library(future)
library(ggvenn)
library(purrr)
library(pheatmap)
library(fs)


#Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
#it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
Target_name<-"new_UCC"
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/")
dataout<-paste0(datain,Target_name,"/Seurat_Analysis/")

dir_create(dataout)
dir_create(paste0(dataout,"Figures/Dimensional_Reduction"),recurse=TRUE)
dir_create(paste0(dataout,"Figures/VolcanoPlots"),recurse=TRUE)

Seur_target<-readRDS(paste0(datain,Target_name,"_seur.rds"))


########DESEQ for Un-normalized Un integrated Datasets (OPTIONAL)######
Seur_target<- NormalizeData(Seur_target)
Seur_target <- FindVariableFeatures(Seur_target)
Seur_target <- ScaleData(Seur_target)
Seur_target <- RunPCA(Seur_target)



###DIFFERENTIAL EXPRESSION ####
#########Make Volcano plots for all Ctrl vs LPs or Ctrl vs PIC inside of each genetic treatment
#########Also finds the shared significant genes across all 
de.L.list<-list()
de.P.list<-list()
sig.L.names<-list()
sig.P.names<-list()
#Iterate DE Volcano function over each background of the dataset. 
#This way you're getting PIC and LPS DE analysis inside each background
for (i in levels(factor(Seur_target$Background))){
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
saveRDS(de.L.list,paste0(dataout,"vsL_DE.rds"))
saveRDS(de.P.list,paste0(dataout,"vsP_DE.rds"))

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


#Turn each gene list into a .txt file so that It can be input into tfTargets
dir_create(paste0(dataout,"Gene_Lists/"))

for (i in names(P.up.gene_sets)){
  writeLines(unlist(P.up.gene_sets[i]),paste0(dataout,"Gene_Lists/",i,"_up.txt"))
  
}

for (i in names(P.down.gene_sets)){
  writeLines(unlist(P.down.gene_sets[i]),paste0(dataout,"Gene_Lists/",i,"_down.txt"))
  
}


#RUNNING TF TARGETS ON ALL THE TXT files

txt_path<-paste0(dataout,"Gene_Lists")

old_wd<-getwd()

setwd("/home/santi/TF_targets")

for (i in list.files(txt_path,pattern="*.txt",full.names = TRUE)){
  outpath<-gsub(".txt","_tf.csv",i)
  system2(command = "/home/santi/.conda/envs/tf_targets/bin/python",
          args=c("find_TF_regulators.py", paste0("--input=",i)  ,paste0("--output=",outpath)),
          stdout = TRUE,
          stderr = TRUE)
}


setwd(old_wd)

#####Summarizing TF target result CSVs#######
#cycle through all csv files and combine them into a single large dataframe
all_tfs<-data.frame()
for (i in list.files(paste0(dataout,"Gene_Lists/"),pattern="*.csv")){
  tf_csv<-as.data.frame(read_csv(paste0(dataout,"Gene_Lists/",i),show_col_types = FALSE))
  
  #changing colnames to make it easier to process in R. no spaces no hyphens
  colnames(tf_csv)<-gsub("-| ","_",colnames(tf_csv))
  #Removing all samples with a pvalue > 0.05
  tf_csv<-subset(tf_csv,P_Value<0.05)
  
  #add the origin adn riection as extra columns
  tf_csv$Origin<-gsub("_up_tf.csv|_down_tf.csv","",i)
  if (grepl("up",i)){
    tf_csv$Direction<-"Upregulated"
  } else if  (grepl("down",i)){
    tf_csv$Direction<-"Downregulated"
  }
  
  all_tfs<-rbind(all_tfs,tf_csv)
}

all_tfs<-all_tfs %>% group_by(Origin,Direction) %>% arrange(P_Value,by_group=TRUE)

# get top pvalue entries per origin and direction
top_tfs<-all_tfs[,c(2,5,9,10,11,12)] %>% 
  arrange(P_Value) %>%
  group_by(Origin,Direction) %>% 
  slice_head(n=10)

saveRDS(top_tfs,paste0(dataout,"Gene_Lists/Top_TF_Grouped.rds"))  
saveRDS(all_tfs,paste0(dataout,"Gene_Lists/All_TF.rds"))

#Calculate the jaccard index between the control DE significant genes####
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


###########Dimensional Reduction########

# PCA and visualization
Seur_target <- RunPCA(Seur_target, verbose = FALSE)
VizDimLoadings(Seur_target, dims = 1:2, reduction = "pca")
Idents(Seur_target) <- "Treatment"

DimPlot(Seur_target, reduction = "pca") 
ggsave(paste0(dataout,"Figures/Dimensional_Reduction/TreatmentPCA.png"))

#Clustering and UMAP visualization
Seur_target <- FindNeighbors(Seur_target, dims = 1:10, verbose = FALSE)
Seur_target <- FindClusters(Seur_target, verbose = FALSE)


Seur_target <- RunUMAP(Seur_target, dims = 1:10, verbose = FALSE)

#Plotting UMAP grouping by treatments
DimPlot(Seur_target, reduction = "umap",group.by="Treatment",label=FALSE)+
labs(title=paste0(Target_name,": UMAP"))
ggsave(paste0(dataout,"Figures/Dimensional_Reduction/TreatmentUmap.png"))

#Plotting UMAP grouping by sample
DimPlot(Seur_target, reduction = "umap",group.by="Sample",label=FALSE)+
  labs(title=paste0(Target_name,": UMAP"))
ggsave(paste0(dataout,"Figures/Dimensional_Reduction/SampleUmap.png"))

saveRDS(Seur_target,paste0(dataout,Target_name,"_final.rds"))





###COMPARING ALL PIC SAMPLES AGAINST ZSCCP####

#########Make Volcano plots for all Ctrl vs LPs or Ctrl vs PIC inside of each genetic treatment
#########Also finds the shared significant genes across all 

#new data out
dataout<-paste0(dataout,"PIC_ONLY/")
dir_create(paste0(dataout,"Figures/VolcanoPlots"),recursive = TRUE)
de.P.list<-list()
sig.P.names<-list()
#Iterate DE Volcano function over each background of the dataset. 
for (i in setdiff(levels(factor(Seur_target$Background)),"C")){
  Cname=paste0("ZSCCP")
  Pname=paste0("ZSC",i,"P")
  de.P<-DEVolcano(Seur_target,Pname,Cname,"Sample")
  sig.P.names<-append(sig.P.names,de.P[de.P$expression!="NS",7])
  de.P.list[[i]]<-de.P
}
#makes a unqiue lsit of all significant genes across all backgrounds
all.sig.P<-unique(sig.P.names)
#saves a list fo dataframes, each one with the DE results for each treatment vs its isogenic ctrl
saveRDS(de.P.list,paste0(dataout,"vsP_DE.rds"))

#####Make a list of dataframe of only the upregulated or downregulated genes across samples
##and separate a list of genes 


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

#Turn each gene list into a .txt file so that It can be input into tfTargets
for (i in names(P.up.gene_sets)){
  writeLines(unlist(P.up.gene_sets[i]),paste0(dataout,"Gene_Lists/",i,"_up.txt"))
  
}

for (i in names(P.down.gene_sets)){
  writeLines(unlist(P.down.gene_sets[i]),paste0(dataout,"Gene_Lists/",i,"_down.txt"))
  
}

txt_path<-paste0(dataout,"Gene_Lists")

old_wd<-getwd()

setwd("/home/santi/TF_targets")

for (i in list.files(txt_path,pattern="*.txt",full.names = TRUE)){
  outpath<-gsub(".txt","_tf.csv",i)
  system2(command = "/home/santi/.conda/envs/tf_targets/bin/python",
          args=c("find_TF_regulators.py", paste0("--input=",i)  ,paste0("--output=",outpath)),
          stdout = TRUE,
          stderr = TRUE)
}


setwd(old_wd)


#####Summarizing TF target result CSVs
#cycle through all csv files and combine them into a single large dataframe
all_tfs<-data.frame()
for (i in list.files(paste0(dataout,"Gene_Lists/"),pattern="*.csv")){
  tf_csv<-as.data.frame(read_csv(paste0(dataout,"Gene_Lists/",i),show_col_types = FALSE))
  
  #changing colnames to make it easier to process in R. no spaces no hyphens
  colnames(tf_csv)<-gsub("-| ","_",colnames(tf_csv))
  #Removing all samples with a pvalue > 0.05
  tf_csv<-subset(tf_csv,P_Value<0.05)
  
  #add the origin adn riection as extra columns
  tf_csv$Origin<-gsub("_tf_up.csv|_tf_down.csv","",i)
  if (grepl("up",i)){
    tf_csv$Direction<-"Upregulated"
  } else if  (grepl("down",i)){
    tf_csv$Direction<-"Downregulated"
  }
  
  all_tfs<-rbind(all_tfs,tf_csv)
}

# get top pvalue entries per origin and direction
top_tfs<-all_tfs[,c(2,5,9,10,11,12)] %>% 
  arrange(P_Value) %>%
  group_by(Origin,Direction) %>% 
  slice_head(n=15)

saveRDS(top_tfs,paste0(dataout,"Gene_Lists/Top_TF_Grouped.rds"))  
saveRDS(all_tfs,paste0(dataout,"Gene_Lists/All_TF.rds"))

##FUNCTIONS LIVE BELOW!!!####


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
  ggsave(paste0(dataout,"Figures/VolcanoPlots/",Target1 ,".png"), width = 6, height = 4, dpi = 300)
  print(volcano)
  invisible(volcano)
  return(deFile)
  
  
}

