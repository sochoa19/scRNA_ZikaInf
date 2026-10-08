#SCT METACELL TEST

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)


Target_name<-"new_URN_int_sct"


seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/",Target_name,"_final.rds"))


# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed, k nearest neighbors parameter and network type for reproducibility
set.seed(12345)

net.type<-"signed"

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)


###Create wild type seurat object
sub_seurat_obj<-subset(seurat_obj, subset= Treatment!="L")

#Set up subsetted object for WGCNA
sub_seurat_obj <- SetupForWGCNA(
  sub_seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "sub_isolate" # the name of the hdWGCNA experiment
)



result_df<-data.frame()
mc_tally<-0
#a loop to cycle through all the k parameter sizes
for (i in c(60)){
  #a loop to cycle through all the max_shared values
  for (j in c(15)){
    
    k.parameter<-i
    max_shared<-j
    
    
    sub_seurat_obj <- MetacellsByGroups(
      seurat_obj = sub_seurat_obj,
      group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
      reduction = 'pca', # select the dimensionality reduction to perform KNN on
      k = k.parameter, # nearest-neighbors parameter
      max_shared = max_shared, # maximum number of shared cells between two metacells
      ident.group = 'Sample' # set the Idents of the metacell seurat object
    )
    
    metacell_metadata <- sub_seurat_obj@misc$sub_isolate$wgcna_metacell_obj@meta.data
    
    #count the number of metacells per sample
    mc_count<-metacell_metadata %>% dplyr::count(Sample)
    
    mc_count$k<-k.parameter
    mc_count$shared<-max_shared
    
    result_df<-rbind(result_df,mc_count)
    mc_tally<-mc_tally+1
    print(mc_tally)
  }
}

result_df_change <- result_df %>%
  group_by(k, shared) %>%
  mutate(
    difference =  n[Sample=="ZSCJP"][1]/n
  ) %>%
  ungroup()

result_df_change[result_df_change$Sample=="ZSC3P",]

ggplot(data=result_df_change[result_df_change$Sample=="ZSC3P",],aes(x=k,y=shared,fill=difference))+
  geom_tile(color = "white", linewidth = 0.5) +          # Adds clean white borders between tiles
  scale_fill_viridis_c(option = "viridis") +
  geom_text(aes(label = n))+ # Labels each tile
  labs(
    title = "Ratio of ZSCJP metacells over ZSC3P metacells",
    x = "K parameter",
    y = "Max shared cells between clusters"
  ) +
  theme_minimal() +                                 # Clean, minimal presentation theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Tilts X labels if text overlaps
  )

ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Metacell_Analysis/ZSC3P_Diff.png"))


summ_result<-result_df_change %>% group_by(k,shared) %>% summarise(meandiff=mean(difference))

ggplot(data=summ_result,aes(x=k,y=shared,fill=meandiff))+
  geom_tile(color = "white", linewidth = 0.5) +          # Adds clean white borders between tiles
  scale_fill_viridis_c(option = "viridis") +
labs(
  title = "Average Ratio of ZSCJP metacells over all sample metacells",
  x = "K parameter",
  y = "Max shared cells between clusters"
) +
  theme_minimal() +                                 # Clean, minimal presentation theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Tilts X labels if text overlaps
  )

ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Metacell_Analysis/Mean_Diff.png"))





result_list<-list()
k<-c(20,20,25,25,60)
shared<-c(11,15,15,5,15)

#a loop to cycle through all the k parameter sizes
for (m in 1:5){
 i<-k[m]
 j<-shared[m]
    result_list<-append(result_list,Metacell_analysis(sub_seurat_obj,i,j))
}

chunked_list <- split(result_list, ceiling(seq_along(result_list) / 6))
result_df<- as.data.frame(do.call(rbind, chunked_list))

result_df <- data.frame(lapply(result_df, unlist))
colnames(result_df)<-c("Metacell_number","Avg_grain","Avg_INV","Avg_Separation","K_parameter","max_shared_cells")

result_df<-  result_df %>%
  mutate(
    K_paremeter = as.numeric(K_parameter),
    max_shared_cells = as.numeric(max_shared_cells)
  )

dir.create(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Metacell_Analysis/"))

Metrics<-c("Metacell_number","Avg_INV","Avg_Separation")
for (i in Metrics){
  #Heatmap for the different Metacell quality statistics
  ggplot(result_df, aes(x = K_parameter, y = max_shared_cells, fill = .data[[i]])) +
    geom_tile(color = "white", linewidth = 0.5) +          # Adds clean white borders between tiles
    scale_fill_viridis_c(option = "viridis") +        # Colorblind-friendly continuous color scale
    #geom_text(aes(label = .data[[i]]), color = "white", size = 2) + # Labels each tile
    labs(
      title = paste0(Target_name," ",gsub("_"," ",i), " across metacell creation parameters"),
      x = "K parameter",
      y = "Max shared cells between clusters",
      fill = gsub("_","",i)
    ) +
    theme_minimal() +                                 # Clean, minimal presentation theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1) # Tilts X labels if text overlaps
    )
  ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Metacell_Analysis/",i,".png"))
}

#Heatmap for the different Metacell quality statistics
ggplot(result_df, aes(x = K_parameter, y = max_shared_cells, fill = .data[[Metrics[1]]])) +
  geom_tile(color = "white", linewidth = 0.5) +          # Adds clean white borders between tiles
  scale_fill_viridis_c(option = "viridis") +        # Colorblind-friendly continuous color scale
  geom_text(aes(label = .data[[Metrics[1]]]), color = "white", size = 2) + # Labels each tile
  labs(
    title = paste0(Target_name," ",gsub("_"," ",Metrics[1]), " across metacell creation parameters"),
    x = "K parameter",
    y = "Max shared cells between clusters",
    fill = gsub("_","",Metrics[1])
  ) +
  theme_minimal() +                                 # Clean, minimal presentation theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Tilts X labels if text overlaps
  )
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Metacell_Analysis/Metacell_number.png"))

#Function that calculates INV, and metacell separation from metacell construction based on k and max shared
Metacell_analysis<-function(wt_seurat_obj,k.parameter,max_shared){
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
  
  ####Metacell Analysis #####
  #get metacell details
  metacell_metadata <- wt_seurat_obj@misc$sub_isolate$wgcna_metacell_obj@meta.data
  #get merged cell numbers per metacell
  metacell_metadata$cell_counts <- lengths(strsplit(metacell_metadata$cells_merged, split = ","))
  #take the first 30 PCA componentsof each cell in the wt_seurat_obj
  sc_embeddings <- Embeddings(wt_seurat_obj, reduction = "pca")[, 1:30]
  centroid<-data.frame()
  for (i in rownames(metacell_metadata)){
    #take all the cells that make up metacell i
    cells<-strsplit(metacell_metadata[i,4],split=",")[[1]]
    #calcualte the mean PCA across the cells inside the i metacell, this will be the centroid coordinates
    PCAmean<-t(colMeans(sc_embeddings[cells,]))
    
    rownames(PCAmean) <- i
    #add the centroid coordinates to teh centroid df with teh cluster name as the rowname
    centroid<- rbind(centroid, PCAmean)
  }
  
  # Calculate ALL pairwise Euclidean distances between metacells
  pairwise_distances <- dist(centroid, method = "euclidean")
  
  # Convert to a readable data frame
  dist_df <- as.matrix(pairwise_distances)
  diag(dist_df)<-Inf
  #add teh minimum distance to another metacell to teh metadata. This is the value for separability
  metacell_metadata$mindist <- apply(dist_df, 1, FUN = min)
  
  
  count_matrix <- GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "counts")
  min_umi = 40
  percentile = 0.95
  inv_results<-list()
  for (i in rownames(metacell_metadata)){
    #get a list of all the cells inside metacell i
    cells<-strsplit(metacell_metadata[i,4],split=",")[[1]]
    #subset the matrix count with the cells insidemetacell i
    mc_counts<-count_matrix[,cells]
    
    # Filter for genes with sufficient total UMIs in this metacell
    gene_totals <- rowSums(mc_counts)
    valid_genes <- which(gene_totals >= min_umi)
    mc_counts_filtered <- mc_counts[valid_genes, , drop = FALSE]
    
    # Compute mean and variance per gene
    gene_means <- rowMeans(mc_counts_filtered)
    gene_vars <- rowSums((mc_counts_filtered - gene_means)^2) / (length(cells) - 1)
    
    # Compute index of dispersion (Variance-to-Mean Ratio)
    valid_means <- gene_means > 0
    normalized_variance <- gene_vars[valid_means] / gene_means[valid_means]
    
    # Extract the 95th percentile
    inv_results <- append(inv_results,quantile(normalized_variance, probs = percentile, na.rm = TRUE))
  }
  metacell_metadata$INV<-unlist(inv_results, use.names = FALSE) 
  
  #things I want to plot for: avg INV, avg min centroid distance,average cell per metacell, number of metacells, 
  #x and y are the k parameters and the shared cell number. 
  grain<-mean(metacell_metadata$cell_counts)
  INV<-mean(metacell_metadata$INV)
  separation<-mean(metacell_metadata$mindist)
  metacell_num<-nrow(metacell_metadata)
  
  metacell_analysis_results<-c(metacell_num,grain,INV,separation,k.parameter,max_shared)
  
  
  
  
  return(metacell_analysis_results)
}








