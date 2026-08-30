#OMNIBUS CODE AN ATTEMPT TO PLACE ALL SCRIPTS INTO A SINGLE SCRIPT> EACH PIECE AS ITS OWN FUNCTION

##ALL LIBRARIES NEEDED FOR ALL SCRIPTS

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(scatterpie)
library(ggplot2)
library(ggvenn)
library(pheatmap)
library(ggalluvial)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

#general help
library(fs)
library(purrr)
library(readr)
library(tidyr)
library(ggrepel)
library(tibble)
library(Matrix)
library(stringr)
#GO ANNOTATION
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)



#Set the name of the Seurat object.rds that you wish to work with. don't include _seur.rds

Target_name<-"new_UCC_int"

Attempt_number<-1
#IS THIS SEURAT OBJECT NORMALIZED?
Norm_Check<-FALSE



#Function Used by Seurat Analysis to create a volcano plots and DE Gene lists.
#Target 1 is the sample of interest, Target 2 is teh baseline
#Focus is teh group across which you're comparing. All of our DE tests are across samples so far
DEVolcano<-function(SeurFile,Target1,Target2,Focus,dataout){
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

###SEURAT ANALYSIS########
#Returns teh ANalysed Seurat Object. Saves it as an RDS as well
#Saves all THe DE volcano plot figures alongside teh .txt files and TF_targets results
Seurat_analysis<-function(Target_name,Norm_Check){
  #SETTING DIRECTORIES AND SEURAT OBJECTS
  #Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
  #it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/")
  dataout<-paste0(datain,Target_name,"/Seurat_Analysis/")
  
  dir_create(dataout)
  dir_create(paste0(dataout,"Figures/Dimensional_Reduction"),recurse=TRUE)
  dir_create(paste0(dataout,"Figures/VolcanoPlots"),recurse=TRUE)
  
  Seur_target<-readRDS(paste0(datain,Target_name,"_seur.rds"))
  
  #Normalize the Seurat Object if it is un-Normalized
  if(Norm_Check==TRUE){
    
    ########DESEQ for Un-normalized Un integrated Datasets (OPTIONAL)######
    Seur_target<- NormalizeData(Seur_target)
    Seur_target <- FindVariableFeatures(Seur_target)
    Seur_target <- ScaleData(Seur_target)
    Seur_target <- RunPCA(Seur_target)
    
  }
  
  
  
  
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
    de.L<-DEVolcano(Seur_target,Lname,Cname,"Sample",dataout)
    de.P<-DEVolcano(Seur_target,Pname,Cname,"Sample",dataout)
    sig.L.names<-append(sig.L.names,de.L[de.L$expression!="NS",7])
    sig.P.names<-append(sig.P.names,de.P[de.P$expression!="NS",7])
    de.L.list[[i]]<-de.L
    de.P.list[[i]]<-de.P
  }
  #makes a unqiue list of all significant genes across all backgrounds
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
  
  
  
  
  
  
  
  ###TF_Targets for IB DE####
  
  print('Make sure you have downloaded TF Targets and created a python environment able to run it')
  
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
  
  tf_targets_dir<-"/home/santi/TF_targets"
  tf_target_env_dir<-"/home/santi/.conda/envs/tf_targets/bin/python"
  setwd(tf_targets_dir)
  
  for (i in list.files(txt_path,pattern="*.txt",full.names = TRUE)){
    outpath<-gsub(".txt","_tf.csv",i)
    system2(command =  tf_target_env_dir,
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
  
  
  ###DIMENSIONAL REDUCTION####
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
  
  
  
  
  ###DIFFERENTIAL EXPERSSION FOR PIC ONLY ACROSS BACKGROUNDS####
  
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
    de.P<-DEVolcano(Seur_target,Pname,Cname,"Sample",dataout)
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
  
  dir_create(dataout,"Gene_Lists")
  
  
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
  
  
  
  
 return(Seur_target) 
}

#DE AND TF COMPARISONS
DE_Multi_Comp<-function(Target_name){
  
  #Populate Seur_target with whichever Seurat object you want to run DE and UMAP viz on
  #it will be used downstream for all the analyses, it will also pull the Name of the object to append to plot titles
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/")
  dataout<-paste0(datain,"Comparison/Venn_Categories/")
  
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
  

  saveRDS(Venn_list,paste0(dataout,"All_Venn_List.rds"))
  #Venn_list<-readRDS(paste0(dataout,"All_Venn_List.rds"))
  #Making .txt lists of each background and category. To be used in Tf targets later
  
  for (j in names(Venn_list)){
    for (i in levels(factor(VennDF$Classification))) {
      writeLines(unlist(subset(Venn_list[[j]],Classification==i)[,"gene"]),paste0(dataout,j,"_",i,"_genelist.txt")) 
    }
  }
  
  fullpath<-paste0(dataout)
  
  #Change directory to Tf_targets and tehn change back
  old_dir <- getwd()
  setwd("/home/santi/TF_targets")
  
  
  for ( i in list.files(paste0(dataout),pattern = "*.txt")){
    j<-gsub("genelist.txt","TF.csv",i)
    result<-system2(command = "/home/santi/.conda/envs/tf_targets/bin/python",
                    args=c("/home/santi/TF_targets/find_TF_regulators.py", paste0("--input=",fullpath,i)  ,paste0("--output=",fullpath,j)),
                    stdout = TRUE,
                    stderr = TRUE)
  }
  
  setwd(old_dir)
  
  #Reading in all TF target csv and combining them into a single Dataframe
  all_tfs<-data.frame()
  for (i in list.files(paste0(dataout),pattern="*.csv")){
    #Skip all CSVS that are of Non signifcant genes, these are out of univers for our purposes
    if(grepl("NS",i)){
      next
    }
    
    tf_csv<-as.data.frame(read_csv(paste0(dataout,i),show_col_types = FALSE))
    
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
  #
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
  
  saveRDS(top_tfs,paste0(dataout,"Top_TF_Grouped.rds"))  
  saveRDS(tf_result,paste0(dataout,"All_TF.rds"))
  
  #
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
  
  ggsave(paste0(dataout,"Upregulated_Gene_Fates.png"))
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
  
  ggsave(paste0(dataout,"Downregulated_Gene_Fates.png"))
}

##WT_Only_WGCNA
WT_ONLY_WGCNA<-function(Target_name,Attempt_Number){
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
  dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/","Attempt_",Attempt_Number)
  
  #make sure the data_out folder is still there
  dir_create(dataout)
  #make sure the Figure data
  dir_create(paste0(dataout,"/Figures"))
  
  seurat_obj<-readRDS(paste0(datain,"/",Target_name,"_final.rds"))
  
  
  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())
  
  # set random seed, k nearest neighbors parameter and network type for reproducibility
  set.seed(12345)
  k.parameter<-60
  max_shared<-15
  net.type<-"signed"
  merge_height<-0.20
  deep_split<-3
  mod_size<-60
  
  parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_Number,merge_height,deep_split,mod_size)
  saveRDS(parameter_list,paste0(dataout,"/Parameter_list.rds"))
  
  # optionally enable multithreading
  enableWGCNAThreads(nThreads = 4)
  
  
  ###Create wild type seurat object
  wt_seurat_obj<-subset(seurat_obj, subset= Background=="C")
  
  #Set up subsetted object for WGCNA
  wt_seurat_obj <- SetupForWGCNA(
    wt_seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "WT_isolate" # the name of the hdWGCNA experiment
  )
  
  #Creating metacells using the samples as the biosamples
  #Would it be better to construct them out of the same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
  wt_seurat_obj <- MetacellsByGroups(
    seurat_obj = wt_seurat_obj,
    group.by = c("Background"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = k.parameter, # nearest-neighbors parameter
    max_shared = max_shared, # maximum number of shared cells between two metacells
    ident.group = 'Background' # set the Idents of the metacell seurat object
  )
  
  
  # normalize metacell expression matrix:
  wt_seurat_obj <- NormalizeMetacells(wt_seurat_obj)
  
  
  #frac <- prop.table(table(
  #  wt_seurat_obj$Treatment[cells_in_metacell]
  #))
  
  mc<-GetMetacellObject(wt_seurat_obj)
  
  ###UMAP OF METACELLS ONLY####
  
  dim(mc)
  head(mc@meta.data)
  colnames(mc@meta.data)
  mc <- GetMetacellObject(wt_seurat_obj)
  mc <- FindVariableFeatures(mc)
  mc <- ScaleData(mc)
  mc <- RunPCA(mc)
  mc <- RunUMAP(mc, dims = 1:20)
  
  cell_meta <- wt_seurat_obj@meta.data
  
  frac_df <- lapply(seq_len(nrow(mc)), function(i){
    
    cells <- strsplit(mc$cells_merged[i], ",")[[1]]
    
    tab <- prop.table(table(cell_meta[cells, "Treatment"]))
    
    data.frame(
      metacell = rownames(mc@meta.data)[i],
      Control = ifelse("C" %in% names(tab), tab["C"], 0),
      LPS     = ifelse("L" %in% names(tab), tab["L"], 0),
      PIC     = ifelse("P" %in% names(tab), tab["P"], 0)
    )
  })
  
  frac_df <- bind_rows(frac_df)
  
  frac_df<-frac_df[!is.na(frac_df$metacell),]
  
  umap <- Embeddings(mc, "umap") |>
    as.data.frame() |>
    tibble::rownames_to_column("metacell")
  
  plot_df <- left_join(umap, frac_df, by = "metacell")
  
  ggplot(plot_df) +
    geom_scatterpie(
      aes(x = umap_2, y = umap_1),
      cols = c("Control", "LPS", "PIC"),
      pie_scale = 1.5
    ) +
    xlab("umap_1")+
    ylab("umap_2")+
    labs(title=paste0(Target_name,"_",Attempt_Number,":WT Metacell UMAP"))+
    coord_equal() +
    theme_classic()
  
  ggsave(paste0(dataout,"/Figures/MetacellUmap.png"))
  
  
  DimPlot(wt_seurat_obj,group.by="Treatment")+
    labs(title = paste0(Target_name,"_",Attempt_Number,":WT UMAP"))
  
  ggsave(paste0(dataout,"/Figures/WTcellUmap.png"))
  
  
  
  ##########WGCNA and correlation of modules####
  #Run hdWGCNA on only the WT cells. THis will allow for the discovery of
  #expected module response to LPS and PIC
  wt_seurat_obj <- SetDatExpr(
    wt_seurat_obj,
    assay = 'RNA', # using RNA assay
    layer = 'data' # using normalized data
  )
  
  
  # Test different soft powers:
  wt_seurat_obj <- TestSoftPowers(
    wt_seurat_obj,
    networkType = net.type # you can also use "unsigned" or "signed hybrid"
  )
  
  # plot the results:
  plot_list <- PlotSoftPowers(wt_seurat_obj)
  
  
  # assemble with patchwork
  wrap_plots(plot_list, ncol=2)
  
  #Consturct the WGCNA Network
  wt_seurat_obj <- ConstructNetwork(
    wt_seurat_obj,
    networktype=net.type,
    tom_name = 'WT_isolate' , # name of the topological overlap matrix written to disk
    #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
    deepSplit = deep_split,
    minModuleSize = mod_size,
    mergeCutHeight = merge_height,
    overwrite_tom=TRUE,
    store_tom_in_seurat=TRUE
    
    
  )
  PlotDendrogram(wt_seurat_obj, main='WT isolate hdWGCNA Dendrogram')
  
  
  p1 <- DimPlot(wt_seurat_obj, group.by='Treatment') +
    umap_theme() +
    ggtitle('WT  Isolate') 
  
  
  p2 <- DimPlot(seurat_obj, group.by='Treatment') +
    umap_theme() +
    ggtitle('Full data') 
  
  
  p1 | p2
  
  
  
  # Project modules from query to reference dataset
  seurat_obj <- ProjectModules(
    seurat_obj = seurat_obj,
    seurat_ref = wt_seurat_obj,
    # vars.to.regress = c(), # optionally regress covariates when running ScaleData
    #group.by.vars = "Sample", # column in seurat_query to run harmony on. We already ran harmony on our dataset
    wgcna_name_proj="projected", # name of the new hdWGCNA experiment in the query dataset
    wgcna_name = "WT_isolate" # name of the hdWGCNA experiment in the ref dataset
  )
  
  
  #plot all module MEs across all cells onto the UMAP
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='MEs', # plot the MEs
    order=TRUE # order so the points with highest MEs are on top
  )
  
  ME<-GetMEs(seurat_obj)
  # stitch together with patchwork
  wrap_plots(plot_list, ncol=6)
  ggsave(paste0(dataout,"/Figures/ModuleUMAP.png"),width=10,height=8)
  
  ##correlelogram to check module correlation
  ModuleCorrelogram(seurat_obj,features="MEs")
  ggsave(paste0(dataout,"/Figures/Correlelogram.png"))
  
  
  
  
  ##Correlating modules to samples using a dot plot
  #getting the MEs
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # add MEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  
  
  
  #Make a dot plot correlating Module Eigengenes against all treatments across a single background
  for (i in unique(seurat_obj$Background)){
    z<-DotPlot(subset(x=seurat_obj,subset=Background==i), features=mods, group.by = 'Treatment')
    z <- z +
      RotatedAxis() +
      scale_color_gradient2(high='red', mid='grey95', low='blue') +
      ggtitle(paste0("WT module correlation in ", i, " background"))
    ggsave(paste0(dataout,"/Figures/Background_",i,"Dotplot.png"),width=10,height=4)
    
  }
  
  #Make a dot plot correlating Module Eigengenes against all backgrounds across a single treatment
  for (i in unique(seurat_obj$Treatment)){
    z<-DotPlot(subset(x=seurat_obj,subset=Treatment==i), features=mods, group.by = 'Background')
    z <- z +
      RotatedAxis() +
      scale_color_gradient2(high='red', mid='grey95', low='blue') +
      ggtitle(paste0("WT module correlation in ", i, " background"))
    ggsave(paste0(dataout,"/Figures/Treatment_",i,"Dotplot.png"),width = 10, height=5)
    
  }
  
  #How to 
  #z$data <- z$data %>%
  #  group_by(features.plot) %>%
  #  mutate(
  #    avg.exp.scaled =
  #      avg.exp.scaled -
  #      avg.exp.scaled[id == "C"]
  #  ) %>%
  #  ungroup()
  
  #z
  
  ##############GO TERMS##################
  #reading in the WGCNAed seurat objects. botht the reference and the query
  
  #######RUNNING GO on the modules and on the genes with kme above 0.5
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  #get modules from wild type isolate reference seuarat object
  modules<-GetModules(wt_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  
  #Compiles entrez ids for the universe of possible gene
  universe.symb <- bitr(
    modules$gene_name,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  #takes any genes still in ENSEMBL format and converts them into entrez ids. HIGH ERROR RATE LOOK INTO IT!
  universe.ensembl <- bitr(
    modules[grepl("ENSG",modules$gene_name),1],
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  
  univ.combined<-append(universe.symb$ENTREZID, universe.ensembl$ENTREZID)
  
  #
  
  #loop for creating GO BP terms dot plots
  #initialize a list to populate with modules with no significant GO terms
  noGO<-list()
  
  ontol<-"BP" #This can be changed depending on what ontology database you want to use
  dir_create(paste0(dataout, "/Figures/GO/",ontol),recurse=TRUE)
  
  for (i in mods){
    
    currcolor<-i #sets the name of the module to be put through GO enrichment analysis
    genelist<-modules[modules$color==currcolor,1] #populates genelist with all the genes that make up this module
    
    #Convert all the gene symbols into entrez ids. 
    gene.symb <- bitr(
      genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(genelist[grepl("ENSG",genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      gene.combined<-gene.symb$ENTREZID
      
    }else{
      gene.ensembl<- bitr(
        genelist[grepl("ENSG",genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      gene.combined <-append(gene.symb$ENTREZID,gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
      
    }
    
    
    qcutoff=0.05 #sets the qvalue cutoff to be used in GO enrichment 0.05 is standard
    #GO enrichment 
    ego <- enrichGO(
      gene = gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    #if there are no GO enrichment results that pass the qvalue cutoff add the module name to noGO list
    if(sum(ego@result$qvalue<=qcutoff)==0){
      print(paste0("The ",currcolor, " module does not have any significant terms associated to it"))
      noGO<-append(noGO,currcolor)
      next
    }
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ", currcolor, " module"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",currcolor,"_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    
  }
  capture.output(noGO,file = paste0(dataout,"/Figures/GO/",ontol,"/nonSigModules.txt"))
  
  
  
  ##############HUB GENE DISCOVERY####
  
  
  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'Background', group_name = 'C'
  )
  
  
  p <- PlotKMEs(seurat_obj, ncol=5)
  
  p
  
  hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
  
  head(hub_df)
  
  
  #Saving the final versions of WGCNA analysis. Both the wt_ object containing the metacells and the full dataset.
  saveRDS(seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_WTonly_projected.rds"))
  
  saveRDS(wt_seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_WTonly.rds"))
  
  
  
}


###PIC_ONLY_WGCNA
PIC_ONLY_WGCNA<-function(Target_name,Attempt_Number){
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
  dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/PIC_Only_WGCNA/","Attempt_",Attempt_Number)
  
  #make sure the data_out folder is still there
  dir_create(dataout)
  #make sure the Figure data
  dir_create(paste0(dataout,"/Figures"))
  
  seurat_obj<-readRDS(paste0(datain,"/",Target_name,"_final.rds"))
  
  
  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())
  
  # set random seed, k nearest neighbors parameter and network type for reproducibility
  set.seed(12345)
  k.parameter<-45
  max_shared<-15
  net.type<-"signed"
  merge_height<-0.15
  deep_split<-3
  mod_size<-60
  
  parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_Number,merge_height,deep_split,mod_size)
  saveRDS(parameter_list,paste0(dataout,"/Parameter_list.rds"))
  
  
  # optionally enable multithreading
  enableWGCNAThreads(nThreads = 4)
  
  
  ###Create wild type seurat object
  wt_seurat_obj<-subset(seurat_obj, subset= Background=="C")
  
  ##Create a PIC only wild type seurat object
  wt_seurat_obj<-subset(wt_seurat_obj, subset= Treatment=="P")
  
  
  #Set up subsetted object for WGCNA
  wt_seurat_obj <- SetupForWGCNA(
    wt_seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "WT_isolate" # the name of the hdWGCNA experiment
  )
  
  
  #Creating metacells using the samples as the biosamples
  #Would it be better to construct them out of the same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
  wt_seurat_obj <- MetacellsByGroups(
    seurat_obj = wt_seurat_obj,
    group.by = c("Background"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = k.parameter, # nearest-neighbors parameter
    max_shared = max_shared, # maximum number of shared cells between two metacells
    ident.group = 'Background' # set the Idents of the metacell seurat object
  )
  
  
  # normalize metacell expression matrix:
  wt_seurat_obj <- NormalizeMetacells(wt_seurat_obj)
  
  
  
  ##########WGCNA and correlation of modules####
  #Run hdWGCNA on only the WT cells. THis will allow for the discovery of
  #expected module response to LPS and PIC
  wt_seurat_obj <- SetDatExpr(
    wt_seurat_obj,
    assay = 'RNA', # using RNA assay
    layer = 'data' # using normalized data
  )
  
  
  # Test different soft powers:
  wt_seurat_obj <- TestSoftPowers(
    wt_seurat_obj,
    networkType = net.type # you can also use "unsigned" or "signed hybrid"
  )
  
  # plot the results:
  plot_list <- PlotSoftPowers(wt_seurat_obj)
  
  
  # assemble with patchwork
  wrap_plots(plot_list, ncol=2)
  
  #Consturct the WGCNA Network
  wt_seurat_obj <- ConstructNetwork(
    wt_seurat_obj,
    networktype=net.type,
    tom_name = 'WT_isolate' , # name of the topological overlap matrix written to disk
    #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
    deepSplit = deep_split,
    minModuleSize = mod_size,
    mergeCutHeight = merge_height,
    overwrite_tom=TRUE,
    store_tom_in_seurat=TRUE
  )
  PlotDendrogram(wt_seurat_obj, main='WT PIC hdWGCNA Dendrogram')
  
  #compute all MEs in the full single-cell dataset
  wt_seurat_obj <- ModuleEigengenes(
    wt_seurat_obj
  )
  
  # compute eigengene-based connectivity (kME):
  wt_seurat_obj <- ModuleConnectivity(wt_seurat_obj)
  
  
  p <- PlotKMEs(wt_seurat_obj, ncol=5)
  
  p
  
  #subsetting only the PIC treated samples acrosss all backgrounds
  p_seurat_obj<-subset(seurat_obj, subset= Treatment=="P")
  
  
  p1 <- DimPlot(seurat_obj, group.by='Treatment') +
    umap_theme() +
    ggtitle('Full data') 
  
  
  p2 <- DimPlot(p_seurat_obj, group.by='Background') +
    umap_theme() +
    ggtitle('PIC only') 
  
  
  p1 | p2
  
  
  
  # Project modules from query to reference dataset
  p_seurat_obj <- ProjectModules(
    seurat_obj = p_seurat_obj,
    seurat_ref = wt_seurat_obj,
    # vars.to.regress = c(), # optionally regress covariates when running ScaleData
    #group.by.vars = "Sample", # column in seurat_query to run harmony on. We already ran harmony on our dataset
    wgcna_name_proj="projected", # name of the new hdWGCNA experiment in the query dataset
    wgcna_name = "WT_isolate" # name of the hdWGCNA experiment in the ref dataset
  )
  
  
  #plot all module MEs across all cells onto the UMAP
  plot_list <- ModuleFeaturePlot(
    p_seurat_obj,
    features='MEs', # plot the MEs
    order=TRUE # order so the points with highest MEs are on top
  )
  
  # stitch together with patchwork
  wrap_plots(plot_list, ncol=6)
  ggsave(paste0(dataout,"/Figures/ModuleUMAP.png"),width=10,height=8)
  
  
  ##correlelogram to check module correlation
  ModuleCorrelogram(p_seurat_obj,features="MEs")
  
  
  
  ##Correlating modules to samples using a dot plot
  #getting the MEs
  MEs <- GetMEs(p_seurat_obj, harmonized=FALSE)
  modules <- GetModules(p_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # add MEs to Seurat meta-data:
  p_seurat_obj@meta.data <- cbind(p_seurat_obj@meta.data, MEs)
  
  z<-DotPlot(subset(x=p_seurat_obj,subset=Treatment=="P"), features=mods, group.by = 'Background')
  z <- z +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    ggtitle(paste0("WT module correlation in PIC Treatment"))
  ggsave(paste0(dataout,"/Figures/Background_Dotplot.png"),width = 10, height=5)
  
  
  
  
  
  
  
  #######RUNNING GO on the modules and on the genes with kme above 0.5
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  #get modules from wild type isolate reference seurat object
  modules<-GetModules(wt_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  
  #Compiles entrez ids for the universe of possible gene
  universe.symb <- bitr(
    modules$gene_name,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  #takes any genes still in ENSEMBL format and converts them into entrez ids. HIGH ERROR RATE LOOK INTO IT!
  universe.ensembl <- bitr(
    modules[grepl("ENSG",modules$gene_name),1],
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  
  univ.combined<-append(universe.symb$ENTREZID, universe.ensembl$ENTREZID)
  
  #
  
  #loop for creating GO BP terms dot plots
  #initialize a list to populate with modules with no significant GO terms
  noGO<-list()
  
  topGO<-list()
  
  ontol<-"BP" #This can be changed depending on what ontology database you want to use
  dir_create(paste0(dataout, "/Figures/GO/",ontol),recurse=TRUE)
  
  for (i in mods){
    
    currcolor<-i #sets the name of the module to be put through GO enrichment analysis
    genelist<-modules[modules$color==currcolor,1] #populates genelist with all the genes that make up this module
    
    #Convert all the gene symbols into entrez ids. 
    gene.symb <- bitr(
      genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(genelist[grepl("ENSG",genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      gene.combined<-gene.symb$ENTREZID
      
    }else{
      gene.ensembl<- bitr(
        genelist[grepl("ENSG",genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      gene.combined <-append(gene.symb$ENTREZID,gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
      
    }
    
    
    qcutoff=0.05 #sets the qvalue cutoff to be used in GO enrichment 0.05 is standard
    #GO enrichment 
    ego <- enrichGO(
      gene = gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    #if there are no GO enrichment results that pass the qvalue cutoff add the module name to noGO list
    if(sum(ego@result$qvalue<=qcutoff)==0){
      print(paste0("The ",currcolor, " module does not have any significant terms associated to it"))
      noGO<-append(noGO,currcolor)
      next
    }
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ", currcolor, " module"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",currcolor,"_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    topGO<-append(topGO,paste(head(ego@result$Description,n=5),collapse=","))
  }
  capture.output(noGO,file = paste0(dataout,"/Figures/GO/",ontol,"/nonSigModules.txt"))
  
  
  
  #saving out the fully anlayzed wt_seruat object _in this case the PIC_WT_
  saveRDS(wt_seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
  
  #saving the WGCNA projection onto ahe all PIC subset
  saveRDS(p_seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_PIC_WT_projected_allPIC.rds"))
  
  #####Module preservation across backgrounds####
  
  #Adapted Complete WGCNA function to check module preservation across all backgrounds
  #will output the complete module preservation df pres as an rds for each treatment
  wt_seurat_obj<-readRDS(paste0(dataout,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
  p_seurat_obj<-readRDS(paste0(dataout,"/",Target_name,"WGCNA_final_PIC_WT_projected_allPIC.rds"))
  dir_create(paste0(dataout,"/Module_preservation"))
  
  
  ## WT to BULK MODULE PRESERVATION
  # set expression matrix for reference dataset
  wt_expr <- GetDatExpr(wt_seurat_obj)
  
  
  
  # a loop that iterates over all teh background in the suerat object and cehcks for module preservation.
  #it will save all teh mp objects under the Module preservation directory, to be compared later
  for (i in unique(p_seurat_obj$Background)){
    
    #Set expression matrix for P dataset one background at a time
    # Subset to the cells you want to use as the query
    query <- subset(
      p_seurat_obj,
      subset = Background == i
    )
    
    
    
    #Set up subsetted object for WGCNA
    query <- SetupForWGCNA(
      query,
      gene_select = "fraction", # the gene selection approach
      fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
      wgcna_name = "projected" # the name of the hdWGCNA experiment
    )
    
    
    #Creating metacells for the query dataset as well to compare metcells vs metacells
    query <- MetacellsByGroups(
      seurat_obj = query,
      group.by = c("Background"), # specify the columns in seurat_obj@meta.data to group by
      reduction = 'pca', # select the dimensionality reduction to perform KNN on
      k = k.parameter, # nearest-neighbors parameter
      max_shared = max_shared, # maximum number of shared cells between two metacells
      ident.group = 'Background' # set the Idents of the metacell seurat object
    )
    
    
    # normalize metacell expression matrix:
    query <- NormalizeMetacells(query)
    
    
    
    query <- SetDatExpr(
      query,
      assay = "RNA",
      slot = "data",
      wgcna_name = "projected"
    )
    
    back_expr<-GetDatExpr(query)
    
    # set up the modules
    modules <- GetModules(wt_seurat_obj)
    
    ref_modules <- list(ref = modules$module)
    
    # set up multiExpr:
    setLabels <- c("ref", "query")
    multiExpr <- list(
      ref = list(data=wt_expr),
      query = list(data=back_expr)
    )
    
    
    
    #Calculate module preservation from PIC WT sc data onto other background PIC
    mp <- WGCNA::modulePreservation(
      multiExpr,
      ref_modules,
      referenceNetworks = 1,
      nPermutations = 150 # set this to whatever number is suitable 
    )
    
    saveRDS(mp,paste0(dataout,"/Module_preservation/",i,"mod_pres.rds"))
    
  }
  
  #Pulling al the mp files and summarizing the results in a single graph
  pres_summ<-list()
  median_summ<-list()
  cor_summ<-list()
  
  for (i in unique(p_seurat_obj$Background)){
    df<-readRDS(paste0(dataout,"/Module_preservation/",i,"mod_pres.rds"))
    pres <- df$preservation$Z[[1]][[2]]
    #figure out what to do with median rank. maybe
    pres$MedianRank<-df$quality$observed$ref.ref$inColumnsAlsoPresentIn.query$medianRank.qual
    pres_summ[[i]]<-pres$Zsummary.pres
    median_summ[[i]]<-pres$MedianRank
    cor_summ[[i]]<-df$preservation$observed$ref.ref$inColumnsAlsoPresentIn.query$cor.cor
  }
  
  pres_summ<-as.data.frame(pres_summ)
  rownames(pres_summ)<-rownames(pres)
  colnames(pres_summ)<-gsub("X","",colnames(pres_summ))
  
  plot_df <- pres_summ %>%
    tibble::rownames_to_column("Module") %>%
    pivot_longer(
      cols = -Module,
      names_to = "Background",
      values_to = "Zsummary"
    )
  
  
  plot_df$Preservation <- cut(
    plot_df$Zsummary,
    breaks = c(-Inf, 2, 10, Inf),
    labels = c("Not preserved", "Moderately", "Strongly")
  )
  
  ggplot(plot_df,
         aes(x = Module,
             y = Zsummary,
             color = Background)) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_hline(yintercept = c(2, 10),
               linetype = "dashed",
               color = c("orange", "darkgreen")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Module",
      y = "Zsummary"
    )
  
  ggsave(paste0(dataout,"/Module_preservation/ZSumm.png"))
  
  
  #correlation summary
  cor_summ<-as.data.frame(cor_summ)
  rownames(cor_summ)<-rownames(pres)
  colnames(cor_summ)<-gsub("X","",colnames(cor_summ))
  
  
  
  
  plot_df <- cor_summ %>%
    tibble::rownames_to_column("Module") %>%
    pivot_longer(
      cols = -Module,
      names_to = "Background",
      values_to = "Cor.cor"
    )
  
  
  
  ggplot(plot_df,
         aes(x = Module,
             y = Cor.cor,
             color = Background)) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Module",
      y = "cor.cor score"
    )
  
  ggsave(paste0(dataout,"/Module_preservation/CorrelationSumm.png"))
  #Module Summary
  
  modsumm<-data.frame(
    modules=mods,
    modsize=sapply(mods, function(i) sum(modules$module == i)))
  
}


#COMPARING WT MODULES TO PIC SUBMODULES
Module_Compare<-function(Target_name,Attempt_Number1,Attempt_Number2){
  #READ IN BOTH MODULE GENE LISTS
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/")
  WT_obj<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/",Target_name,"WGCNA_final_WTonly.rds"))
  PIC_obj<-readRDS(paste0(datain,"PIC_Only_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
  
  wtmods<-GetModules(WT_obj)
  picmods<-GetModules(PIC_obj)
  sum(wtmods$gene_name %in%picmods$gene_name)
  merged_mods<-merge(wtmods[,1:2],picmods[,1:2], by="gene_name")
  colnames(merged_mods)<-c("gene_name","WT_Modules","PIC_Modules")
  merged_mods_grey<-merged_mods[! merged_mods$PIC_Modules=="grey",]
  alluvial_data_grey<-merged_mods_grey %>% dplyr::count(WT_Modules,PIC_Modules) 
  merged_mods<-merged_mods[!(merged_mods$WT_Modules=="grey" | merged_mods$PIC_Modules=="grey"),]
  alluvial_data<-merged_mods %>% dplyr::count(WT_Modules,PIC_Modules) 
  
  
  
  module_colors <- c(
    turquoise = "turquoise",
    blue      = "blue",
    brown     = "brown",
    yellow    = "yellow",
    green     = "green",
    red       = "red",
    black     = "grey40",
    pink      = "pink",
    magenta   = "magenta",
    purple    = "purple",
    tan       = "tan",
    salmon    = "salmon",
    grey      = "grey",
    greenyellow="greenyellow"
  )
  
  
  
  ggplot(data=alluvial_data,
    aes(axis1=WT_Modules,axis2=PIC_Modules,y=n))+
    geom_alluvium(aes(fill=WT_Modules))+
    geom_stratum(
      aes(fill = after_stat(stratum)),
      color = "black"
    ) +

    geom_text(stat="stratum",aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("WT_Modules", "PIC_Modules"),
                     expand = c(0.15, 0.05)) +
    scale_fill_manual(
      values = module_colors
    ) +
    scale_x_discrete(
      limits = c("WT_Modules", "PIC_Modules"),
      expand = c(0.15, 0.05)
    ) +
    theme_cowplot()
  ggsave(paste0(datain,"PIC_Only_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny.png"),height=10,width=8)

  
  
  ggplot(data=alluvial_data_grey,
         aes(axis1=WT_Modules,axis2=PIC_Modules,y=n))+
    geom_alluvium(aes(fill=WT_Modules))+
    geom_stratum(
      aes(fill = after_stat(stratum)),
      color = "black"
    ) +
    
    geom_text(stat="stratum",aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("WT_Modules", "PIC_Modules"),
                     expand = c(0.15, 0.05)) +
    scale_fill_manual(
      values = module_colors
    ) +
    scale_x_discrete(
      limits = c("WT_Modules", "PIC_Modules"),
      expand = c(0.15, 0.05)
    ) +
    theme_cowplot()
  ggsave(paste0(datain,"PIC_Only_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny_wGrey.png"),height=10,width=8)
  
  }


##TF TARGETS FOR MODULE
MODULE_TF<-function(Target_name,Attempt_Number){
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seruat_OUT/",Target_name,"/Isolate_WGCNA/",Attempt_Number,"/")
  dataout<-paste0()
}

Seurat_analysis("new_UCC_int",FALSE)
DE_Multi_Comp("new_UCC_int")
WT_ONLY_WGCNA("new_UCC_int",1)
PIC_ONLY_WGCNA("new_UCC_int",1)

print("voila")
Seurat_analysis("UCC_int",FALSE)
DE_Multi_Comp("UCC_int")
WT_ONLY_WGCNA("UCC_int",1)
PIC_ONLY_WGCNA("UCC_int",1)
Module_Compare("UCC_int",1,1)
