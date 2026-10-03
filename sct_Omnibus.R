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
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

#general help
library(fs)
library(purrr)
library(readr)
library(tidyr)
library(tibble)
library(Matrix)
library(stringr)
#GO ANNOTATION
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)



#Set the name of the Seurat object.rds that you wish to work with. don't include _seur.rds

Target_name<-"new_URN_int_sct"

Attempt_number<-1
#IS THIS SEURAT OBJECT UnNORMALIZED?
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
  dir_create(paste0(dataout,"Figures/GO"),recurse=TRUE)
  
  Seur_target<-readRDS(paste0(datain,Target_name,"_seur.rds"))
  
  #Normalize the Seurat Object if it is un-Normalized
  if(Norm_Check==TRUE){
    
    ########DESEQ for Un-normalized Un integrated Datasets (OPTIONAL)######
    Seur_target<- NormalizeData(Seur_target)
    Seur_target <- FindVariableFeatures(Seur_target)
    Seur_target <- ScaleData(Seur_target)
    Seur_target <- RunPCA(Seur_target)
    
  }
  
  Seur_target<-PrepSCTFindMarkers(Seur_target)
  
  
  ###DIFFERENTIAL EXPRESSION ####
  #########Make Volcano plots for all Ctrl vs LPs or Ctrl vs PIC inside of each genetic treatment
  #########Also finds the shared significant genes across all 
  de.P.list<-list()
  sig.P.names<-list()
  #Iterate DE Volcano function over each background of the dataset. 
  #This way you're getting PIC and LPS DE analysis inside each background
  for (i in levels(factor(Seur_target$Background))){
    Cname=paste0("ZSC",i,"C")
    Pname=paste0("ZSC",i,"P")
    de.P<-DEVolcano(Seur_target,Pname,Cname,"Sample",dataout)
    sig.P.names<-append(sig.P.names,de.P[de.P$expression!="NS",7])
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
      filter(expression == "Up") %>%
      pull(gene) %>%
      unique()
  })
  L.up.WT.genes<-L.up.gene_sets[["C"]]
  
  L.down.gene_sets <- map(de.L.list, ~{
    .x %>%
      filter(expression == "Down") %>%
      pull(gene) %>%
      unique()
  })
  L.down.WT.genes<-L.down.gene_sets[["C"]]
  
  
  
  P.up.gene_sets <- map(de.P.list, ~{
    .x %>%
      filter(expression == "Up") %>%
      pull(gene) %>%
      unique()
  })
  P.up.WT.genes<-P.up.gene_sets[["C"]]
  
  P.down.gene_sets <- map(de.P.list, ~{
    .x %>%
      filter(expression == "Down") %>%
      pull(gene) %>%
      unique()
  })
  P.down.WT.genes<-P.down.gene_sets[["C"]]
  
  
  
  ####GO TERMS FOR ALL DE ####
  
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  #getting all genes queried in the DE for each sample
  for (i in names(de.P.list)){
    gene_universe<-de.P.list[[i]]$gene
    
    
    #Compiles entrez ids for the universe of possible gene
    universe.symb <- bitr(
      gene_universe,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #takes any genes still in ENSEMBL format and converts them into entrez ids. HIGH ERROR RATE LOOK INTO IT!
    universe.ensembl <- bitr(
      gene_universe[grepl("ENSG",gene_universe)],
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    univ.combined<-append(universe.symb$ENTREZID, universe.ensembl$ENTREZID)
    
    #converting all genes 
    up_genelist<-P.up.gene_sets[[i]]
    
    #Convert all the gene symbols into entrez ids. 
    up_gene.symb <- bitr(
      up_genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(up_genelist[grepl("ENSG",up_genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      up_gene.combined<-up_gene.symb$ENTREZID
      
    }else{
      up_gene.ensembl<- bitr(
        up_genelist[grepl("ENSG",up_genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      up_gene.combined <-append(up_gene.symb$ENTREZID,up_gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
    }
    
    #converting all downregualted genes
    down_genelist<-P.down.gene_sets[[i]]
    
    #Convert all the gene symbols into entrez ids. 
    down_gene.symb <- bitr(
      down_genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(down_genelist[grepl("ENSG",down_genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      down_gene.combined<-down_gene.symb$ENTREZID
      
    }else{
      down_gene.ensembl<- bitr(
        down_genelist[grepl("ENSG",down_genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      down_gene.combined <-append(down_gene.symb$ENTREZID,down_gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
    }
    
    noGO<-list()
    
    ontol<-"BP" #This can be changed depending on what ontology database you want to use
    dir_create(paste0(dataout, "/Figures/GO/",ontol),recurse=TRUE)
    qcutoff=0.05 #sets the qvalue cutoff to be used in GO enrichment 0.05 is standard
    
    #GO enrichment 
    up_ego <- enrichGO(
      gene = up_gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    
    #GO enrichment 
    down_ego <- enrichGO(
      gene = down_gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(up_ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ZSC", i, " upregulated in PIC response"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",i,"_UP_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    
    
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(down_ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ZSC", i, " downregulated in PIC response"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",i,"_down_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    
    
    
  }
  
  
  
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
  dir_create(paste0(dataout,"Figures/GO"),recurse=TRUE)
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
      filter(expression == "Up") %>%
      pull(gene) %>%
      unique()
  })
  P.up.WT.genes<-P.up.gene_sets[["C"]]
  
  P.down.gene_sets <- map(de.P.list, ~{
    .x %>%
      filter(expression == "Down") %>%
      pull(gene) %>%
      unique()
  })
  P.down.WT.genes<-P.down.gene_sets[["C"]]
  
  
  #RUN GO ON EACH DIFFERNTIAL EXPRESSION LIST
  
  #getting all genes queried in the DE for each sample
  for (i in names(de.P.list)){
    gene_universe<-de.P.list[[i]]$gene
    
    
    #Compiles entrez ids for the universe of possible gene
    universe.symb <- bitr(
      gene_universe,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #takes any genes still in ENSEMBL format and converts them into entrez ids. HIGH ERROR RATE LOOK INTO IT!
    universe.ensembl <- bitr(
      gene_universe[grepl("ENSG",gene_universe)],
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    univ.combined<-append(universe.symb$ENTREZID, universe.ensembl$ENTREZID)
    
    #converting all genes 
    up_genelist<-P.up.gene_sets[[i]]
    
    #Convert all the gene symbols into entrez ids. 
    up_gene.symb <- bitr(
      up_genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(up_genelist[grepl("ENSG",up_genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      up_gene.combined<-up_gene.symb$ENTREZID
      
    }else{
      up_gene.ensembl<- bitr(
        up_genelist[grepl("ENSG",up_genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      up_gene.combined <-append(up_gene.symb$ENTREZID,up_gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
    }
    
    #converting all downregualted genes
    down_genelist<-P.down.gene_sets[[i]]
    
    #Convert all the gene symbols into entrez ids. 
    down_gene.symb <- bitr(
      down_genelist,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    #Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
    if (length(grep(paste(down_genelist[grepl("ENSG",down_genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
    ){
      down_gene.combined<-down_gene.symb$ENTREZID
      
    }else{
      down_gene.ensembl<- bitr(
        down_genelist[grepl("ENSG",down_genelist)],
        fromType= "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
        
      )
      
      down_gene.combined <-append(down_gene.symb$ENTREZID,down_gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
    }
    
    
    ontol<-"BP" #This can be changed depending on what ontology database you want to use
    dir_create(paste0(dataout, "/Figures/GO/",ontol),recurse=TRUE)
    qcutoff=0.05 #sets the qvalue cutoff to be used in GO enrichment 0.05 is standard
    
    #GO enrichment 
    up_ego <- enrichGO(
      gene = up_gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    
    #GO enrichment 
    down_ego <- enrichGO(
      gene = down_gene.combined,
      universe = univ.combined,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontol, #This can be changed depending on what ontology database you want to use
      readable = TRUE,
      qvalueCutoff = qcutoff
      
    )
    
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(up_ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ZSC", i, " upregulated against WT PIC response"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",i,"_UP_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    
    
    #create a dotplot of the 10 most significant enriched terms for that module
    dotplot(down_ego, showCategory = 10)+
      ggtitle(paste0(ontol,"GO enrichment of ZSC", i, " downregulated against WT PIC response"))
    ggsave(paste0(dataout,"/Figures/GO/",ontol,"/",i,"_down_GO.png"),height=1800,width=2200,units= "px", dpi= 300)
    
    
    
  }
  
  
  
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
  
  all_tfs<-all_tfs %>% group_by(Origin,Direction) %>% arrange(P_Value)
  
  # get top pvalue entries per origin and direction
  top_tfs<-all_tfs[,c(2,5,9,10,11,12)] %>% 
    arrange(P_Value) %>%
    group_by(Origin,Direction) %>% 
    slice_head(n=15)
  
  saveRDS(top_tfs,paste0(dataout,"Gene_Lists/Top_TF_Grouped.rds"))  
  saveRDS(all_tfs,paste0(dataout,"Gene_Lists/All_TF.rds"))
  
  
  
  
  return(Seur_target) 
}



# set random seed, k nearest neighbors parameter and network type for reproducibility
set.seed(12345)
k.parameter<-45
max_shared<-15
net.type<-"signed"
merge_height<-0.15
deep_split<-3
mod_size<-60

parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_number,merge_height,deep_split,mod_size)

##WT_Only_WGCNA
WT_ONLY_WGCNA<-function(Target_name,Attempt_number){
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
  dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/","Attempt_",Attempt_number)
  
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
  
  parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_number,merge_height,deep_split,mod_size)
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
    group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = k.parameter, # nearest-neighbors parameter
    max_shared = max_shared, # maximum number of shared cells between two metacells
    ident.group = 'Sample' # set the Idents of the metacell seurat object
  )
  
  
  # normalize metacell expression matrix:
  wt_seurat_obj <- NormalizeMetacells(wt_seurat_obj)
  
  
  
  ###UMAP OF METACELLS ONLY####
  mc<-GetMetacellObject(wt_seurat_obj)
  
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
    labs(title=paste0(Target_name,"_",Attempt_number,":WT Metacell UMAP"))+
    coord_equal() +
    theme_classic()
  
  ggsave(paste0(dataout,"/Figures/MetacellUmap.png"))
  
  
  DimPlot(wt_seurat_obj,group.by="Treatment")+
    labs(title = paste0(Target_name,"_",Attempt_number,":WT UMAP"))
  
  ggsave(paste0(dataout,"/Figures/WTcellUmap.png"))
  
  
  
  ##########WGCNA and correlation of modules####
  #Run hdWGCNA on only the WT cells. THis will allow for the discovery of
  #expected module response to LPS and PIC
  wt_seurat_obj <- SetDatExpr(
    wt_seurat_obj,
    assay = 'SCT', # using RNA assay
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



All_PIC_CTRL_WGCNA<-function(Target_name,Attempt_number,parameter_list){
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
  dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/All_PIC_WGCNA/","Attempt_",Attempt_number)
  
  k.parameter<-as.numeric(parameter_list[1])
  max_shared<-as.numeric(parameter_list[2])
  net.type<-parameter_list[3]
  merge_height<-as.numeric(parameter_list[6])
  deep_split<-as.numeric(parameter_list[7])
  mod_size<-as.numeric(parameter_list[8])
  
  
  
  #make sure the data_out folder is still there
  dir_create(dataout)
  
  saveRDS(parameter_list,paste0(dataout,"/Parameter_list.rds"))
  #make sure the Figure data
  dir_create(paste0(dataout,"/Figures"))
  
  seurat_obj<-readRDS(paste0(datain,"/",Target_name,"_final.rds"))
  
  
  # using the cowplot theme for ggplot
  theme_set(theme_cowplot())
  
  saveRDS(parameter_list,paste0(dataout,"/Parameter_list.rds"))
  
  
  # optionally enable multithreading
  enableWGCNAThreads(nThreads = 4)
  
  
  
  ##Create a PIC and CTRL subset of the seurat object. All backgrounds. 
  sub_seurat_obj<-subset(seurat_obj, subset= Treatment!="L")
  
  
  #Set up subsetted object for WGCNA
  sub_seurat_obj <- SetupForWGCNA(
    sub_seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "sub_isolate" # the name of the hdWGCNA experiment
  )
  
  
  #Creating metacells using the samples as the biosamples
  #Would it be better to construct them out of the same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
  sub_seurat_obj <- MetacellsByGroups(
    seurat_obj = sub_seurat_obj,
    group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = k.parameter, # nearest-neighbors parameter
    max_shared = max_shared, # maximum number of shared cells between two metacells
    ident.group = 'Sample' # set the Idents of the metacell seurat object
  )
  
  
  # normalize metacell expression matrix:
  sub_seurat_obj <- NormalizeMetacells(sub_seurat_obj)
  
  
  ###UMAP OF METACELLS ONLY####
  mc<-GetMetacellObject(sub_seurat_obj)
  
  dim(mc)
  head(mc@meta.data)
  colnames(mc@meta.data)
  mc <- FindVariableFeatures(mc)
  mc <- ScaleData(mc)
  mc <- RunPCA(mc)
  mc <- RunUMAP(mc, dims = 1:20)
  
  
  frac_df <- lapply(seq_len(nrow(mc)), function(i){
    
    cells <- strsplit(mc$cells_merged[i], ",")[[1]]
    
    tab <- prop.table(table(cell_meta[cells, "Treatment"]))
    
    data.frame(
      metacell = rownames(mc@meta.data)[i],
      Control = ifelse("C" %in% names(tab), tab["C"], 0),
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
      cols = c("Control", "PIC"),
      pie_scale = 1.5
    ) +
    xlab("umap_1")+
    ylab("umap_2")+
    labs(title=paste0(Target_name,"_",Attempt_number,":WT Metacell UMAP"))+
    coord_equal() +
    theme_classic()
  
  sum(frac_df$Control>0.7 | frac_df$Control<0.3)
  
  DimPlot(mc)
  
  
  
  ##########WGCNA and correlation of modules####
  #Run hdWGCNA on only the WT cells. THis will allow for the discovery of
  #expected module response to PIC
  sub_seurat_obj <- SetDatExpr(
    sub_seurat_obj,
    assay = 'SCT', # using SCT assay
    layer = 'data' # using normalized data
  )
  
  
  # Test different soft powers:
  sub_seurat_obj <- TestSoftPowers(
    sub_seurat_obj,
    networkType = net.type # you can also use "unsigned" or "signed hybrid"
  )
  
  # plot the results:
  plot_list <- PlotSoftPowers(sub_seurat_obj)
  
  
  # assemble with patchwork
  wrap_plots(plot_list, ncol=2)
  
  #Consturct the WGCNA Network
  sub_seurat_obj <- ConstructNetwork(
    sub_seurat_obj,
    networktype=net.type,
    tom_name = 'sub_isolate' , # name of the topological overlap matrix written to disk
    #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
    deepSplit = deep_split,
    minModuleSize = mod_size,
    mergeCutHeight = merge_height,
    overwrite_tom=TRUE,
    store_tom_in_seurat=TRUE
  )
  PlotDendrogram(sub_seurat_obj, main='WT PIC hdWGCNA Dendrogram')
  
  #compute all MEs in the full single-cell dataset
  sub_seurat_obj <- ModuleEigengenes(
    sub_seurat_obj
  )
  
  # compute eigengene-based connectivity (kME):
  sub_seurat_obj <- ModuleConnectivity(sub_seurat_obj)
  
  
  p <- PlotKMEs(sub_seurat_obj, ncol=5)
  
  p
  
  
  #plot all module MEs across all cells onto the UMAP
  plot_list <- ModuleFeaturePlot(
    sub_seurat_obj,
    features='MEs', # plot the MEs
    order=TRUE # order so the points with highest MEs are on top
  )
  
  # stitch together with patchwork
  wrap_plots(plot_list, ncol=6)
  ggsave(paste0(dataout,"/Figures/ModuleUMAP.png"),width=10,height=8)
  
  
  ##correlelogram to check module correlation
  ModuleCorrelogram(sub_seurat_obj,features="MEs")
  
  
  
  ##Correlating modules to samples using a dot plot
  #getting the MEs
  MEs <- GetMEs(sub_seurat_obj, harmonized=FALSE)
  modules <- GetModules(sub_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # add MEs to Seurat meta-data:
  sub_seurat_obj@meta.data <- cbind(sub_seurat_obj@meta.data, MEs)
  
  z<-DotPlot(subset(x=sub_seurat_obj,subset=Treatment=="P"), features=mods, group.by = 'Background')
  z <- z +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    ggtitle(paste0("PIC&CTRL correlation in PIC Treatment"))
  ggsave(paste0(dataout,"/Figures/Background_Dotplot.png"),width = 10, height=5)
  
  
  
  z<-DotPlot(subset(x=sub_seurat_obj,subset=Treatment=="C"), features=mods, group.by = 'Background')
  z <- z +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    ggtitle(paste0("WT module correlation in CTRL Treatment"))
  
  z
  
  
  #######RUNNING GO on the modules and on the genes with kme above 0.5
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  #get modules from wild type isolate reference seurat object
  modules<-GetModules(sub_seurat_obj)
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
  
  ontol<-"MF" #This can be changed depending on what ontology database you want to use
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
  saveRDS(sub_seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
  
  #Module DE####
  
  
  # set up an empty dataframe for the DMEs
  DMEs <- data.frame()
  
  # loop through the clusters
  for(cur_sample in setdiff(levels(factor(sub_seurat_obj$Sample)),"ZSCJP")){
    
    # identify barcodes for group1 and group2 in eadh cluster
    group1 <- sub_seurat_obj@meta.data %>% subset(Sample == "ZSCJP") %>% rownames
    group2 <- sub_seurat_obj@meta.data %>% subset(Sample == cur_sample) %>% rownames
    
    # run the DME test
    cur_DMEs <- FindDMEs(
      sub_seurat_obj,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use='wilcox',
      pseudocount.use=0.01, # we can also change the pseudocount with this param
      wgcna_name = 'sub_isolate',
      logfc.threshold = 1
      
    )
    
    # add the cluster info to the table
    cur_DMEs$sample <- cur_sample
    
    # append the table
    DMEs <- rbind(DMEs, cur_DMEs)
  }
  
  ###PLOT THE DME as heatmap
  # get the modules table:
  modules <- GetModules(sub_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # make a copy of the DME table for plotting
  plot_df <- DMEs
  
  # set the factor level for the modules so they plot in the right order:
  plot_df$module <- factor(as.character(plot_df$module), levels=mods)
  
  # set a min/max threshold for plotting
  maxval <- 10; minval <- -10
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
  
  # add significance levels
  plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
  
  # change the text color to make it easier to see 
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 5, 'black', 'white')
  
  # make the heatmap with geom_tile
  p <- plot_df %>% 
    ggplot(aes(y=sample, x=module, fill=avg_log2FC)) +
    geom_tile() 
  
  # add the significance levels
  p <- p + 
    geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
  
  # customize the color and theme of the plot
  p <- p + 
    scale_fill_gradient2(low='purple', mid='black', high='yellow') +
    RotatedAxis() +
    theme(
      panel.border = element_rect(fill=NA, color='black', size=1),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.margin=margin(0,0,0,0)
    ) + xlab('') + ylab('') +
    coord_equal()
  p
  
  sub_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/All_PIC_WGCNA/Attempt_1/new_URN_int_sctWGCNA_final_PIC_WT.rds")
  
  #module expression scores
  sub_seurat_obj <- ModuleExprScore(
    sub_seurat_obj,
    n_genes = 25,
    method='UCell'
  )
  
  
  
  # set up an empty dataframe for the DMSs
  DMSs <- data.frame()
  
  # loop through the clusters
  for(cur_sample in setdiff(levels(factor(sub_seurat_obj$Sample)),"ZSCJP")){
    
    # identify barcodes for group1 and group2 in eadh cluster
    group2 <- sub_seurat_obj@meta.data %>% subset(Sample == "ZSCJP") %>% rownames
    group1 <- sub_seurat_obj@meta.data %>% subset(Sample == cur_sample) %>% rownames
    
    # run the DME test
    cur_DMSs <- FindDMEs(
      sub_seurat_obj,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use='wilcox',
      features = "ModuleScores",
      wgcna_name = 'sub_isolate'

    )
    
    # add the cluster info to the table
    cur_DMSs$sample <- cur_sample
    
    # append the table
    DMSs <- rbind(DMSs, cur_DMSs)
  }
  
  ###PLOT THE DME as heatmap
  # get the modules table:
  modules <- GetModules(sub_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  
  # make a copy of the DME table for plotting
  plot_df <- DMSs
  
  # set the factor level for the modules so they plot in the right order:
  plot_df$module <- factor(as.character(plot_df$module), levels=mods)
  
  # set a min/max threshold for plotting
  maxval <- 2; minval <- -2
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
  
  # add significance levels
  plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
  
  # change the text color to make it easier to see 
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.5, 'black', 'white')
  
  # make the heatmap with geom_tile
  p <- plot_df %>% 
    ggplot(aes(y=sample, x=module, fill=avg_log2FC)) +
    geom_tile() 
  
  # add the significance levels
  p <- p + 
    geom_text(label=plot_df$Significance, color=plot_df$textcolor) 
  
  # customize the color and theme of the plot
  p <- p + 
    scale_fill_gradient2(low='purple', mid='black', high='yellow') +
    RotatedAxis() +
    theme(
      panel.border = element_rect(fill=NA, color='black', size=1),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      plot.margin=margin(0,0,0,0)
    ) + xlab('') + ylab('') +
    coord_equal()
  p
  
  
  #hubgenes
  sub_hub_df <- GetHubGenes(sub_seurat_obj, n_hubs = 10)
  

}

###sc_REgNet
#scRegNet tests

library(Seurat)
library(scRegNet)
library(EnsDb.Hsapiens.v86)
library(fs)
library(dplyr)

Target_name<-"new_URN_int_sct"
Attempt_Number1<-1
Attempt_Number2<-1

datain<-(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/"))

seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"_seur.rds"))


#read in the seurat object with the WT module list
temp_seur<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/",Target_name,"WGCNA_final_WTonly.rds"))
wt_mods<-GetModules(temp_seur)

#Read in the seurat object with teh PIC module list
temp_seur<-readRDS(paste0(datain,"All_PIC_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
all_pic_mods<-GetModules(temp_seur)


temp_seur<-NULL

###use all gens from teh seurat, take it before you subset for highly variable genes!
#only do pearson correlation

#How Should I subset my seurat object pull 

re_scale_seurat_obj <- ScaleData(seurat_obj,features=rownames(seurat_obj))


# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "SCT", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "SCT", layer = "scale.data"))


##running for all WT modules
dataout<-paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/scRegNet/")
dir_create(paste0(dataout,"Individual_Module_tfs/"),recurse=TRUE)

#initialize the summary list and the neg and pos associated tfs
all_assoc<-list()
pos_assoc<-data.frame()
neg_assoc<-data.frame()
#iterating over all modules
for ( i in setdiff(levels(factor(wt_mods$module)),"grey")){
  
  mod_list<-wt_mods[wt_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  
  # a message to keep track on what number module you're on
  print(paste0(which(setdiff(levels(factor(wt_mods$module)),"grey")==i)," out of ", length(setdiff(levels(factor(wt_mods$module)),"grey")==i)))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs/",i,".rds"))
  
  all_assoc<-append(all_assoc,assoc_tfs)
  
  neg_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R<0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)
  pos_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R>0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)
  
  #add the module results to teh complete tf dataframe
  neg_assoc<-rbind(neg_assoc,neg_assoc_tfs)
  pos_assoc<-rbind(pos_assoc,pos_assoc_tfs)
  
}

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))
saveRDS(neg_assoc,paste0(dataout,"negative_associated_tfs.rds"))
saveRDS(pos_assoc,paste0(dataout,"positice_associated_tfs.rds"))



####REPEATING FOR PIC SUBMODULES
#SUBSETTING OUT THE LPS TREATED CELLS AND MOST OF THE CONTROL ONES, LEAVING ONLY THE WT
re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment!="L")

#re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment=="P"|Background=="C")

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "SCT", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "SCT", layer = "scale.data"))



##running for all PIC modules
dataout<-paste0(datain,"All_PIC_WGCNA/Attempt_",Attempt_Number2,"/scRegNet/")
dir_create(paste0(dataout,"Individual_Module_tfs/"),recurse=TRUE)

#initialize the summary list and the neg and pos associated tfs
all_assoc<-list()
pos_assoc<-data.frame()
neg_assoc<-data.frame()
#iterating over all modules
for ( i in setdiff(levels(factor(all_pic_mods$module)),"grey")){
  
  mod_list<-all_pic_mods[all_pic_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  # a message to keep track on what number module you're on
  print(paste0(which(setdiff(levels(factor(all_pic_mods$module)),"grey")==i)," out of ", length(setdiff(levels(factor(all_pic_mods$module)),"grey")==i)))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs/",i,".rds"))
  
  
  all_assoc<-rbind(all_assoc,assoc_tfs)
  
  
}

all_assoc<-all_assoc %>% dplyr::filter(Pearson_p_value<0.05) %>% group_by(module) %>% mutate(Direction= case_when(Pearson_R<0 ~ "Down",Pearson_R>0 ~ "Up")) %>% arrange(desc(Pearson_R),Pearson_p_value)

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))

up_associated_tf_summary<- all_assoc %>% dplyr::group_by(module) %>% filter(Pearson_R>0) %>% slice(1:10)
down_associated_tf_summary<- all_assoc %>% dplyr::group_by(module) %>% filter(Pearson_R<0) %>% slice_tail(n=10)
write_csv(up_associated_tf_summary,paste0(dataout,"up_associated_tfs.csv"))
write_csv(down_associated_tf_summary,paste0(dataout,"down_associated_tfs.csv"))

c(unique(down_associated_tfs[down_associated_tfs$module=="blue",3]))





##### MOD COMPARE ####
Attempt_Number1<-1
Attempt_Number2<-1
Target_name<-"new_URN_int_sct"
#PIC object origin,

#READ IN BOTH MODULE GENE LISTS
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/")
WT_obj<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/",Target_name,"WGCNA_final_WTonly.rds"))
PIC_obj<-readRDS(paste0(datain,"All_PIC_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))


wtmods<-GetModules(WT_obj)
picmods<-GetModules(PIC_obj)
sum(wtmods$gene_name %in%picmods$gene_name)
merged_mods<-merge(wtmods[,1:2],picmods[,1:2], by="gene_name")
colnames(merged_mods)<-c("gene_name","WT_Modules","PIC_Modules")
merged_mods_grey<-merged_mods[! merged_mods$PIC_Modules=="grey",]
alluvial_data_grey<-merged_mods_grey %>% dplyr::count(WT_Modules,PIC_Modules) 
merged_mods<-merged_mods[!(merged_mods$WT_Modules=="grey" | merged_mods$PIC_Modules=="grey"),]
alluvial_data<-merged_mods %>% dplyr::count(WT_Modules,PIC_Modules) 


#make a vector of all the module colors and assign their own color to them. 
module_colors<-as.character(unique(c(merged_mods$PIC_Modules,merged_mods$WT_Modules)))
names(module_colors)<-unique(c(merged_mods$PIC_Modules,merged_mods$WT_Modules))





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
ggsave(paste0(datain,"All_PIC_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny.png"),height=10,width=8)



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
ggsave(paste0(datain,"All_PIC_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny_wGrey.png"),height=10,width=8)



##Module network comparisons
tomcheck<-GetTOM(PIC_obj)
picmods<-GetModules(PIC_obj)
royalmod<-picmods[picmods$module=="royalblue",1]
turqmod<-picmods[picmods$module=="turquoise",1]
redmod<-picmods[picmods$module=="red",1]

royalTOM<-(tomcheck[royalmod, royalmod ])
redTOM<-(tomcheck[redmod,redmod ])
turqTOM<-(tomcheck[turqmod,turqmod ])

sumTOM<- colSums(royalTOM, na.rm = TRUE)
tomnames<-names(head(sort(sumTOM, decreasing = TRUE), 25))

# The matrix is mirrored across teh diagonal
#only keep the upper triangular part of the TOM:
royalTOM[upper.tri(royalTOM)] <- NA

# cast the network from wide to long format
cur_network <- royalTOM %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene1 = Var1, gene2 = Var2, weight=value) %>%
  subset(!is.na(weight))


# keep this network before subsetting
cur_network_full <- cur_network 

# keep the top 20% of edges 
edge_percent <- 0.2
cur_network <- cur_network_full %>% 
  dplyr::slice_max(
    order_by = weight, 
    n = round(nrow(cur_network)*edge_percent)
  )

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)


# get the top 25 hub genes for each module
hub_genes <- GetHubGenes(PIC_obj, n_hubs=25) %>% filter(module=="royalblue") %>% .$gene_name
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")

# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight)) + 
  geom_node_point(color='black') +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic') 

p

#Now make a correlation network for these same genes in a sample by sample manner? 

#####TOM GENERATION FOR WT ONLY PIC AND CTRL AND THEN TEH SAME FOR ZSC3####
Attempt_number<-1

datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/","Attempt_",Attempt_number)


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

parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_number,merge_height,deep_split,mod_size)


# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)


###Create wild type seurat object
wildsub_seurat_obj<-subset(seurat_obj, subset= Background=="C" & Treatment %in% c("P","C"))

picsub_seurat_obj<-subset(seurat_obj, subset= Background=="3" & Treatment %in% c("P","C"))

#Set up subsetted object for WGCNA
wildsub_seurat_obj <- SetupForWGCNA(
  wildsub_seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wildsub_isolate" # the name of the hdWGCNA experiment
)

#Creating metacells using the samples as the biosamples
#Would it be better to construct them out of the same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
wildsub_seurat_obj <- MetacellsByGroups(
  seurat_obj = wildsub_seurat_obj,
  group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = k.parameter, # nearest-neighbors parameter
  max_shared = max_shared, # maximum number of shared cells between two metacells
  ident.group = 'Sample' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
wildsub_seurat_obj <- NormalizeMetacells(wildsub_seurat_obj)




##########WGCNA and correlation of modules####
#Run hdWGCNA on only the WT cells. THis will allow for the discovery of
#expected module response to LPS and PIC
wildsub_seurat_obj <- SetDatExpr(
  wildsub_seurat_obj,
  assay = 'SCT', # using RNA assay
  layer = 'data' # using normalized data
)


# Test different soft powers:
wildsub_seurat_obj <- TestSoftPowers(
  wildsub_seurat_obj,
  networkType = net.type # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(wildsub_seurat_obj)


# assemble with patchwork
wrap_plots(plot_list, ncol=2)

#Consturct the WGCNA Network
wildsub_seurat_obj <- ConstructNetwork(
  wildsub_seurat_obj,
  networktype=net.type,
  tom_name = 'wildsub_isolate' , # name of the topological overlap matrix written to disk
  #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
  deepSplit = deep_split,
  minModuleSize = mod_size,
  mergeCutHeight = merge_height,
  overwrite_tom=TRUE,
  soft_power=9,
  store_tom_in_seurat=TRUE
  
  
)

#####FOR PIC ######

#Set up subsetted object for WGCNA
picsub_seurat_obj <- SetupForWGCNA(
  picsub_seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "picsub_isolate" # the name of the hdWGCNA experiment
)

#Creating metacells using the samples as the biosamples
#Would it be better to construct them out of the same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
picsub_seurat_obj <- MetacellsByGroups(
  seurat_obj = picsub_seurat_obj,
  group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = k.parameter, # nearest-neighbors parameter
  max_shared = max_shared, # maximum number of shared cells between two metacells
  ident.group = 'Sample' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
picsub_seurat_obj <- NormalizeMetacells(picsub_seurat_obj)




##########WGCNA and correlation of modules####
#Run hdWGCNA on only the WT cells. THis will allow for the discovery of
#expected module response to LPS and PIC
picsub_seurat_obj <- SetDatExpr(
  picsub_seurat_obj,
  assay = 'SCT', # using RNA assay
  layer = 'data' # using normalized data
)


# Test different soft powers:
picsub_seurat_obj <- TestSoftPowers(
  picsub_seurat_obj,
  networkType = net.type # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(picsub_seurat_obj)


# assemble with patchwork
wrap_plots(plot_list, ncol=2)

#Consturct the WGCNA Network
picsub_seurat_obj <- ConstructNetwork(
  picsub_seurat_obj,
  networktype=net.type,
  tom_name = 'picsub_isolate' , # name of the topological overlap matrix written to disk
  #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
  deepSplit = deep_split,
  minModuleSize = mod_size,
  mergeCutHeight = merge_height,
  overwrite_tom=TRUE,
  soft_power=9,
  store_tom_in_seurat=TRUE
)


#####
wildsubTOM<-GetTOM(wildsub_seurat_obj)

picsubTOM<-GetTOM(picsub_seurat_obj)

wsroyalTOM<-wildsubTOM[royalmod,royalmod]
sumTOM<- colSums(wsroyalTOM, na.rm = TRUE)
wstomnames<-names(head(sort(sumTOM, decreasing = TRUE), 25))

psroyalTOM<-picsubTOM[royalmod,royalmod]
sumTOM<- colSums(psroyalTOM, na.rm = TRUE)
pstomnames<-names(head(sort(sumTOM, decreasing = TRUE), 25))

check<-wsroyalTOM-psroyalTOM
wsroyalTOM[upper.tri(wsroyalTOM)] <- NA
psroyalTOM[upper.tri(psroyalTOM)] <- NA

check<-abs(check)
sumTOM<- colSums(check, na.rm = TRUE)
checktomnames<-names(head(sort(sumTOM, decreasing = TRUE), 25))



# cast the network from wide to long format
cur_network <- psroyalTOM %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene1 = Var1, gene2 = Var2, weight=value) %>%
  subset(!is.na(weight))


# keep this network before subsetting
cur_network_full <- cur_network 

# keep the top 20% of edges 
edge_percent <- 0.2
cur_network <- cur_network_full %>% 
  dplyr::slice_max(
    order_by = weight, 
    n = round(nrow(cur_network)*edge_percent)
  )

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)


# get the top 25 hub genes for each module
hub_genes <- GetHubGenes(PIC_obj, n_hubs=25) %>% filter(module=="royalblue") %>% .$gene_name
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")

# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight)) + 
  geom_node_point(color='black') +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic') 

p


checkhub<-check[hub_genes,hub_genes]
checkhub[upper.tri(checkhub)] <- NA



#####heatmaps of the modules of interest

mc<-GetMetacellObject(PIC_obj) 
mc<-subset(mc,Sample %in% c('ZSCJP', 'ZSC2P'))

expr <- GetAssayData(
  mc,
  assay = "SCT",
  layer = "data"
)

mat <- as.matrix(
  expr[turqmod, , drop = FALSE]
)


mat_scaled <- t(scale(t(mat)))

mat_scaled <- mat_scaled[
  apply(mat_scaled, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]


# Create column annotation showing sample origin
annotation_col <- data.frame(
  Sample = mc$Sample
)

rownames(annotation_col) <- colnames(mat_scaled)

# Optional: control the order of samples in the legend
annotation_col$Sample <- factor(
  annotation_col$Sample,
  levels = c('ZSCJP','ZSC2P')
)



pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  main = paste('turquoise ', "module -", 'WT-PIC and NS2a-PIC'),
  fontsize_row = 7
)
