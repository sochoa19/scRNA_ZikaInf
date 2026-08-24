#COMPARING TF ACTIVITY IN DIFFENT BACKGROUNDS AND TREATMENTS
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
library(fs)


#Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
#it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
Target_name<-"URN"
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/")
dataout<-paste0(datain,"Comparison/")

dir_create(dataout)
dir_create(paste0(dataout,"Figures"),recurse=TRUE)


#puling in the Intra_background DE : IB_DE
IB_DE<-readRDS(paste0(datain,"vsP_DE.rds"))
#Pulling in the WT vs ZSC DE across backgrounds: AB_DE
AB_DE<-readRDS(paste0(datain,"PIC_ONLY/vsP_DE.rds"))

#S
#VENN DIAGRAM OF ALL DIFFERENTIALLY EXPRESSED GENES BY BACKGROUND

#defining Boolean categories based on my notes
# 1. Define your 8 desired outcomes
outcomes <- c(
  "NS", # [FALSE, FALSE, FALSE]
  "A", # [TRUE,  FALSE, FALSE]
  "C", # [FALSE, TRUE,  FALSE]
  "B", # [TRUE,  TRUE,  FALSE]
  "F", # [FALSE, FALSE, TRUE]
  "D", # [TRUE,  FALSE, TRUE]
  "E", # [FALSE, TRUE,  TRUE]
  "G"  # [TRUE,  TRUE,  TRUE]
)

# R index arrays using 1 and 2, where FALSE=1 and TRUE=2
lookup_cube <- array(outcomes, dim = c(2, 2, 2))

Venn_list<-list()
for (i in names(AB_DE)){
  #find a shared gene list across all tested sample
  gene_keep<-AB_DE[[i]][AB_DE[[i]]$gene %in% IB_DE[[i]]$gene,]
  gene_keep<-gene_keep[gene_keep$gene %in% IB_DE[["C"]]$gene,"gene"]
  print(paste0(i, " kept ",length(gene_keep)/nrow(AB_DE[[i]])*100,"% of genes"))
  VennDF<-data.frame(
    gene=gene_keep,
    Back_Treat=IB_DE[[i]][IB_DE[[i]]$gene %in% gene_keep,"expression"],
    WT_Treat=IB_DE[["C"]][IB_DE[["C"]]$gene %in% gene_keep,"expression"],
    WT_Back_Treat=AB_DE[[i]][AB_DE[[i]]$gene %in% gene_keep,"expression"]
  )
  
  #mapping the significant expression results onto teh Venn diagram categories
  lookup_idx <- cbind((VennDF$Back_Treat!="NS")+1,(VennDF$WT_Treat!="NS")+1, (VennDF$WT_Back_Treat!="NS")+1)
  VennDF$Classification<-lookup_cube[lookup_idx]
  #adding the dataframe into a summary List
  Venn_list[[i]]<-VennDF
 
}

subset(Venn_list[["J"]],Classification=="G")$gene

saveRDS(Venn_list,paste0(dataout,"All_Venn_List.rds"))

#Making .txt lists of each background and category. To be used in Tf targets later
dir_create(paste0(dataout,"Venn_Categories"))

for (j in names(Venn_list)){
  for (i in levels(factor(VennDF$Classification))) {
   writeLines(unlist(subset(Venn_list[[j]],Classification==i)[,"gene"]),paste0(dataout,"Venn_Categories/",j,"_",i,"_genelist.txt")) 
  }
}

fullpath<-paste0(dataout,"Venn_Categories/")

#Change directory to Tf_targets and tehn change back
old_dir <- getwd()
setwd("/home/santi/TF_targets")


for ( i in list.files(paste0(dataout,"Venn_Categories/"),pattern = "*.txt")){
  j<-gsub("genelist.txt","TF.csv",i)
  result<-system2(command = "/home/santi/.conda/envs/tf_targets/bin/python",
          args=c("/home/santi/TF_targets/find_TF_regulators.py", paste0("--input=",fullpath,i)  ,paste0("--output=",fullpath,j)),
          stdout = TRUE,
          stderr = TRUE)
  cat(result, sep = "\n")
}

setwd(old_dir)

#Reading in all TF target csv and combining them into a single Dataframe
all_tfs<-data.frame()
for (i in list.files(paste0(dataout,"Venn_Categories/"),pattern="*.csv")){
  #Skip all CSVS that are of Non signifcant genes, these are out of univers for our purposes
  if(grepl("NS",i)){
    next
  }
  
  tf_csv<-as.data.frame(read_csv(paste0(dataout,"Venn_Categories/",i),show_col_types = FALSE))
  
  #changing colnames to make it easier to process in R. no spaces no hyphens
  colnames(tf_csv)<-gsub("-| ","_",colnames(tf_csv))
  
  #adding their origin 
  #Removing all samples with a pvalue > 0.05 and 
  tf_csv<-subset(tf_csv,P_Value<0.05)
  
  #add the origin adn riection as extra columns
  tf_csv$Origin<-sub("_.*","",i)
  tf_csv$Category<-str_sub(i,-8,-8)
  
  all_tfs<-rbind(all_tfs,tf_csv)
}


# get top pvalue entries per origin and direction
top_tfs<-all_tfs[,c(2,5,9,10,11,12)] %>% 
  arrange(P_Value) %>%
  group_by(Origin,Category) %>% 
  slice_sample(n=10)

saveRDS(top_tfs,paste0(dataout,"Venn_Categories/Top_TF_Grouped.rds"))  
saveRDS(all_tfs,paste0(dataout,"Venn_Categories/All_TF.rds"))
