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
  
  #appenda ll teh csvs together
  all_tfs<-rbind(all_tfs,tf_csv)
}

#adding columns that capture the predominant direction of all 3 DF comparisons in a captured TF
#to determine hwether the TF is dowregulating or upregulating the gene in question
library(purrr)
library(tidyr)

#Making each TF target gene its own row
tf_long <- all_tfs %>%
  separate_rows(Overlapping, sep = "[ ;]") %>%
  rename(gene = Overlapping)

#Cmbine all vendfs into a single dataframe and add the venn name as the origin
de_combined <- bind_rows(Venn_list, .id = "Origin")

#Combine tf_long and de_combined so that teh DE data is appended to each tf gene pair
tf_long <- tf_long %>%
  left_join(
    de_combined,
    by = c("Origin", "gene")
  )

tf_result <- tf_long %>%
  count(TF, Origin, Back_Treat, Category, name = "n") %>%
  group_by(TF, Origin,Category) %>%
  mutate(
    total = sum(n),
    fraction = n / total
  ) %>%
  summarise(
    Direction = if (any(fraction >= 0.666)) {
      Back_Treat[which(fraction >= 0.666)[1]]
    } else {
      "conflicted"
    },
    total = first(total),
    max_fraction = max(fraction),
    .groups = "drop"
  )

classify_66 <- function(data, class_col) {
  
  data %>%
    group_by(TF, Origin, Category, .data[[class_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(TF, Origin, Category) %>%
    mutate(
      total = sum(n),
      fraction = n / total
    ) %>%
    summarise(
      "{class_col}_result" := if (any(fraction >= 0.666)) {
        .data[[class_col]][which(fraction >= 0.666)[1]]
      } else {
        "conflicted"
      },
      .groups = "drop"
    )
}
back_treat_result <- classify_66(tf_long, "Back_Treat")

wt_treat_result <- classify_66(tf_long, "WT_Treat")

wt_back_treat_result <- classify_66(tf_long, "WT_Back_Treat")
tf_result <- all_tfs %>%
  left_join(back_treat_result, by = c("TF", "Origin","Category")) %>%
  left_join(wt_treat_result, by = c("TF", "Origin","Category")) %>%
  left_join(wt_back_treat_result, by = c("TF", "Origin","Category"))  

tf_result<- tf_result %>%
  group_by(Origin,Category) %>% 
  arrange(P_Value, .by_group = TRUE)

# get top pvalue entries per origin and Category
top_tfs<-tf_result[,-c(1,3,4,6,7,8)] %>% 
  group_by(Origin,Category) %>% 
  arrange(P_Value, .by_group = TRUE)  %>%
  slice_head(n=15)

saveRDS(top_tfs,paste0(dataout,"Venn_Categories/Top_TF_Grouped.rds"))  
saveRDS(tf_result,paste0(dataout,"Venn_Categories/All_TF.rds"))


#Pulling a subset of genes that are upregulated in WT_Control. These are our canonicla PIC response
canon_active<-unique(IB_DE[["C"]][IB_DE[["C"]]$expression=="Up","gene"])
canon_active<-canon_active[!is.na(canon_active)]
canon_suppressed<-unique(IB_DE[["C"]][IB_DE[["C"]]$expression=="Down","gene"])
canon_suppressed<-canon_suppressed[!is.na(canon_suppressed)]

plot_data<-data.frame()
for (i in names(AB_DE)){
  test_up<-subset(AB_DE[[i]],gene %in% canon_active & expression=="Up" )
  test_down<-subset(AB_DE[[i]],gene %in% canon_active & expression=="Down" )
  test_ns<-subset(AB_DE[[i]],gene %in% canon_active & expression=="NS" )
  test_not_found<-length(canon_active)-nrow(test_up)-nrow(test_down)-nrow(test_ns)
  test_df<-data.frame(origin=c(i,i,i,i),
                      gene_num=c(nrow(test_up),nrow(test_down),nrow(test_ns),test_not_found),
                      Wt_Back=c("Up","Down","No Change","Not Found"))
 plot_data<- rbind(plot_data,test_df)
}

plot_data<- plot_data %>%
  mutate(
    Wt_Back = factor(
      Wt_Back,
      levels = c("Up", "Down", "No Change", "Not Found")
    )
  )

ggplot(plot_data, aes(x = factor(origin), y = gene_num, fill = Wt_Back)) +
  geom_col() +
  labs(
    title=("Upregulated genes in WT vs WT-PIC"),
    x = "Origin",
    y = "Number of genes",
    fill = "DE WT-PIC to ZSC-PIC"
  ) +
  theme_classic()

###canonically supressed gnes in PIC

plot_data<-data.frame()
for (i in names(AB_DE)){
  test_up<-subset(AB_DE[[i]],gene %in% canon_suppressed & expression=="Up" )
  test_down<-subset(AB_DE[[i]],gene %in% canon_suppressed & expression=="Down" )
  test_ns<-subset(AB_DE[[i]],gene %in% canon_suppressed & expression=="NS" )
  test_not_found<-length(canon_suppressed)-nrow(test_up)-nrow(test_down)-nrow(test_ns)
  test_df<-data.frame(origin=c(i,i,i,i),
                      gene_num=c(nrow(test_up),nrow(test_down),nrow(test_ns),test_not_found),
                      Wt_Back=c("Up","Down","No Change","Not Found"))
  plot_data<- rbind(plot_data,test_df)
}

plot_data<- plot_data %>%
  mutate(
    Wt_Back = factor(
      Wt_Back,
      levels = c("Up", "Down", "No Change", "Not Found")
    )
  )

ggplot(plot_data, aes(x = factor(origin), y = gene_num, fill = Wt_Back)) +
  geom_col() +
  labs(
    title=("Downregulated genes in WT vs WT-PIC"),
    x = "Origin",
    y = "Number of genes",
    fill = "DE WT-PIC to ZSC-PIC"
  ) +
  theme_classic()



                        