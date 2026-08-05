#PIC Isolate Metacell and WGCNA Diagnostic

#Metacell Parameter Testing
#We want to knwo wha the best k parameter, min cells and max shared numbers might be to produce the highest
#quality metacell while retaining the most real variance in the data

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)

#miscellaneous packages
library(fs)
library(dplyr)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)


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
  metacell_metadata <- wt_seurat_obj@misc$WT_isolate$wgcna_metacell_obj@meta.data
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
  
  
  count_matrix <- GetAssayData(wt_seurat_obj, assay = DefaultAssay(wt_seurat_obj), layer = "counts")
  min_umi = 40
  percentile = 0.95
  inv_results<-list()
  for (i in rownames(metacell_metadata)){
    #get a list of all the cells inside metacell i
    cells<-strsplit(metacell_metadata[i,4],split=",")[[1]]
    #subset the matrix count with the cells inside metacell i
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



Target_name<-"UCC_Integrated"


seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/",Target_name,"_final.rds"))

#Writing a path to output figures and rds files
data_out<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/PIC_tests/")

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed, k nearest neighbors parameter and network type for reproducibility
set.seed(12345)

net.type<-"signed"

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)


###Create wild type seurat object
wt_seurat_obj<-subset(seurat_obj, subset= Background=="C")

wt_p_seurat_obj<-subset(wt_seurat_obj, subset= Treatment=="P")


#Set up subsetted object for WGCNA
wt_p_seurat_obj <- SetupForWGCNA(
  wt_p_seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "WT_isolate" # the name of the hdWGCNA experiment
)

result_list<-list()
#a loop to cycle through all the k parameter sizes
for (i in c(20,25,30,35,40,45,50,55,60,65,70)){
  #a loop to cycle through all the max_shared values
  for (j in c(5,8,11,14,17,20)){
    result_list<-append(result_list,Metacell_analysis(wt_p_seurat_obj,i,j))
  }
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
  ggsave(paste0(data_out,"/Metacell_tests/",i,".png"))
}


#########WGCNA PARAMETERS######


#########FUNCTION#####

Complete_WGCNA<-function(seurat_obj,bulk_obj,param_list){
  #unbox the parameter list
  k_parameter<-as.double(param_list[1])
  net.type<-param_list[2]
  max_shared<-as.double(param_list[3])
  modulesize<-as.double(param_list[4])
  mergeheight<-as.double(param_list[5])
  deepsplit<-as.double(param_list[6])
  ###Create wild type seurat object
  wt_seurat_obj<-subset(seurat_obj, subset= Background=="C")
  wt_seurat_obj<-subset(wt_seurat_obj, subset= Treatment=="P")
  
  #Set up subsetted object for WGCNA
  wt_seurat_obj <- SetupForWGCNA(
    wt_seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "WT_isolate" # the name of the hdWGCNA experiment
  )
  
  #Creating metacells form teh WT group only, using tested parameters
  wt_seurat_obj <- MetacellsByGroups(
    seurat_obj = wt_seurat_obj,
    group.by = c("Background"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = k_parameter, # nearest-neighbors parameter
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
  
  #Construct the WGCNA Network
  wt_seurat_obj <- ConstructNetwork(
    wt_seurat_obj,
    networktype=net.type,
    tom_name = 'WT_isolate' , # name of the topological overlap matrix written to disk
    #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
    overwrite_tom=TRUE,
    mergeCutHeight = mergeheight,
    minModuleSize = modulesize,
    deepSplit = deepsplit
  )
  
  PlotDendrogram(wt_seurat_obj, main='WT isolate hdWGCNA Dendrogram')
  
  
  
  #get modules from wild type isolate reference seurat object
  modules<-GetModules(wt_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  #Create the module DF that will be populated with all module stats and details
  module.df<-data.frame(name=mods)
  rownames(module.df)<-mods
  
  #Calculate Module Eigengenes in the WT subset seurat object
  wt_seurat_obj<-ModuleEigengenes(wt_seurat_obj)
  
  #Caluclate kme of each module gene
  wt_seurat_obj <- ModuleConnectivity(
    wt_seurat_obj
  )
  
  #save out the hub genes into a dataframe
  hb<-GetHubGenes(wt_seurat_obj, n_hubs = 5)
  
  module.df$Hub_Genes<-NA
  
  #Ensembl Ids of all the modules. 
  
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
  
  ontol<-"BP" #This can be changed depending on what ontology database you want to use
  
  #adding a GO column to module df that will be populated inside the loop
  module.df$GO<-NA
  
  #loop for generating GO enriched terms for each module
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
    #if there are no GO enrichment results that pass the qvalue cutoff write No terms in GO column
    if(sum(ego@result$qvalue<=qcutoff)==0){
      module.df[i,"GO"]<-"No terms"
      next
    }
    #add top 5 GO terms to each row
    tophits<-ego@result$Description[1:5]
    module.df[i,"GO"]<-paste(tophits,collapse=",")
    
    #add top 5 hub genes per module
    module.df[i,"Hub_Genes"]<-paste(hb[hb$module==i,1],collapse=",")
    
  }
  
  ## WT to BULK MODULE PRESERVATION
  
  # set expression matrix for reference dataset
  wt_expr <- GetDatExpr(wt_seurat_obj)
  
  
  #subsetting only the genes that are shared across datasets
  # get the genes that are in common: 
  bulk<-t(bulk_obj)
  genes_bulk <- colnames(bulk)
  genes_sc <- colnames(wt_expr)
  genes_keep <- genes_bulk[genes_bulk %in% genes_sc]
  bulk_kept<- bulk[,genes_keep]
  wt_expr_kept<-wt_expr[,genes_keep]
  
  #call good sample genes to amke sure all the genes being sued are 
  gsg_ref <- goodSamplesGenes(wt_expr_kept, verbose = 3)
  gsg_query <- goodSamplesGenes(bulk_kept, verbose = 3)
  
  #adding up both sets of bad genes
  bad_genes <- !gsg_ref$goodGenes | !gsg_query$goodGenes
  
  #Get rid of the bad genes from both expression datasets
  wt_expr2   <- wt_expr_kept[, !bad_genes]
  bulk_kept2 <- bulk_kept[, !bad_genes]
  
  
  # set up the modules
  modules <- GetModules(wt_seurat_obj)
  
  common_genes <- Reduce(intersect, list(
    colnames(wt_expr2),
    colnames(bulk_kept2),
    modules$gene
  ))
  
  wt_expr3   <- wt_expr2[, common_genes]
  bulk_kept3 <- bulk_kept2[, common_genes]
  
  modules3 <- modules[match(common_genes, modules$gene), ]
  #modules3<-modules3[modules3$module!='grey',-5]
  #modules3 <- droplevels(modules3)
  ref_modules <- list(ref = modules3$module)
  
  # set up multiExpr:
  setLabels <- c("ref", "query")
  multiExpr <- list(
    ref = list(data=wt_expr3),
    query = list(data=bulk_kept3)
  )
  
  
  
  #Calculate module preservation from WT sc data onto bulk dataset
  mp <- WGCNA::modulePreservation(
    multiExpr,
    ref_modules,
    referenceNetworks = 1,
    nPermutations = 150 # set this to whatever number is suitable for you
  )
  pres <- mp$preservation$Z[[1]][[2]]
  
  pres<-pres[!rownames(pres) %in% c("grey","gold"),]
  
  module.df$size<-NA
  module.df$Bulk_Preservation<-NA
  module.df$Zsumm<-NA
  for (i in rownames(pres)){
    module.df[i,"size"]<-pres[i,"moduleSize"]
    module.df[i,"Zsumm"]<-pres[i,"Zsummary.pres"]
    
    if(pres[i,"Zsummary.pres"]>10){
      module.df[i,"Bulk_Preservation"]<-"Strongly Preserved"
    }else if(pres[i,"Zsummary.pres"]>2){
      module.df[i,"Bulk_Preservation"]<-"Weakly Preserved"
    }else{
      module.df[i,"Bulk_Preservation"]<-"Not Preserved"
    }
  }
  return(module.df)
}




#Libraries needed fro WGCNA validation, needed for GO step
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

k_parameter<-45
net.type<-"signed"
max_shared<-15


#Create a datapath string for outgoing objects
data_out2<-paste0(data_out,"WGCNA_tests")
#Create out path if not already existing
dir_create(data_out2)
dir_create(paste0(data_out2,"/modsize_merge_deepsplit"))
dir_create(paste0(data_out2,"/modsize_merge"))
dir_create(paste0(data_out2,"/modsize_deepsplit"))


#Read in the TSV bulk RNA seq object
bulk_target<-"HMC3_INF"

bulk_obj<-readRDS(paste0("/data/bulkRNA/",bulk_target,"/DESeq2_results/bulk_df.rds"))




###########MOD SIZE and MERGE HEIGHT#####

modulesizerange <- c(30,40,50,60)
mergerange <- c(0.15,0.20,0.25)
deep_split_range<-4
for(i in modulesizerange){
  for (j in mergerange){
    #Collect all parameters into a single vector to feed into the collected function
    param_list<-c(k_parameter,net.type,max_shared,i,j,deep_split_range)
    combname<-paste0("modsize",i,"_merge",j,".rds")
    saveRDS(object=Complete_WGCNA(seurat_obj,bulk_obj,param_list),file =paste0(data_out2,"/modsize_merge/",combname) )
  }
}


#summarizing results from the parameter search
param_outs<-list()
for (i in list.files(paste0(data_out2,"/modsize_merge"))){
  df<-as.data.frame(readRDS(paste0(data_out2,"/modsize_merge/",i)))
  pass<-sum(df$Zsumm>=2)
  fly<-sum(df$Zsumm>=10)
  passgene<-sum(df[df$Zsumm>2,4])/sum(df$size)
  modnum<-nrow(df)
  
  param_outs[[gsub(".rds","",i)]]<-c(modnum,
                                     as.numeric((pass/modnum)*100),
                                     as.numeric((fly/modnum)*100),
                                     as.numeric(str_sub(i,start=8,end=9)),
                                     as.numeric(str_sub(i,start=16,end=-5)),
                                     as.numeric(passgene))
}
paramdf<-as.data.frame(param_outs)
rownames(paramdf)<-c("size","Percent_passing","percent_strongly_passing","modsize","mergeheight","Pass_Gene")


# Transpose so each row is one parameter combination
plot_df <- as.data.frame(t(paramdf))

# Plot
ggplot(plot_df,
       aes(x = mergeheight,
           y = modsize,
           size = Percent_passing,
           color = Pass_Gene)) +
  geom_point() +
  scale_color_gradient(    low = "magenta",
                           high = "gold")+
  scale_size_continuous(range = c(2, 10)) +
  xlab("Merge Cut Height")+
  ylab("Minimum Module Size")+
  labs(title="Percentage of Modules and Genes Preserved in Bulk",
       color="%Genes Preserved",
       size="%Modules Preserved")+
  theme_bw()

ggsave(paste0(data_out2,"/modsize_merge/parameter_summary.png"))




###Testing all three ranges
modulesizerange <- c(50,60,70)
mergerange <- c(0.15,0.25)
deep_split_range<-c(3,4)
for(i in modulesizerange){
  for (j in mergerange){
    for (h in deep_split_range){
      #Collect all parameters into a single vector to feed into the collected function
      param_list<-c(k_parameter,net.type,max_shared,i,j,h)
      combname<-paste0("modsize",i,"_merge",j,"deepsplit_",h,".rds")
      saveRDS(object=Complete_WGCNA(seurat_obj,bulk_obj,param_list),file =paste0(data_out2,"/modsize_merge_deepsplit/",combname) )
    }
  }
}


#summarizing results from the parameter search
param_outs<-list()
for (i in list.files(paste0(data_out2,"/modsize_merge_deepsplit"))){
  df<-as.data.frame(readRDS(paste0(data_out2,"/modsize_merge_deepsplit/",i)))
  pass<-sum(df$Zsumm>=2)
  fly<-sum(df$Zsumm>=10)
  passgene<-sum(df[df$Zsumm>2,4])/sum(df$size)
  modnum<-nrow(df)
  
  param_outs[[gsub(".rds","",i)]]<-c(modnum,
                                     as.numeric((pass/modnum)*100),
                                     as.numeric((fly/modnum)*100),
                                     as.numeric(str_sub(i,start=8,end=9)),
                                     as.numeric(str_sub(i,start=-5,end=-5)),
                                     as.numeric(str_sub(i,start=16,end=-16)),
                                     as.numeric(passgene))
}
paramdf<-as.data.frame(param_outs)
rownames(paramdf)<-c("size","Percent_passing","percent_strongly_passing","modsize","deep_split","merge_height","Pass_Gene")

# Transpose so each row is one parameter combination
plot_df <- as.data.frame(t(paramdf))

library(ggcube)

# Plot
ggplot(plot_df,
       aes(x = deep_split,
           y = modsize,
           z = merge_height,
           size = Percent_passing,
           color = Pass_Gene)) +
  geom_point() +
  scale_color_gradient(    low = "magenta",
                           high = "gold")+
  scale_size_continuous(range = c(2, 10)) +
  xlab("Deep Split")+
  ylab("Minimum Module Size")+
  labs(title="Percentage of Modules and Genes Preserved in Bulk",
       color="%Genes Preserved",
       size="%Modules Preserved")+
  coord_3d()+
  theme_bw()

ggsave(paste0(data_out2,"/modsize_merge_deepsplit/parameter_summary.png"),width=7, height=7)

test<-readRDS(paste0(data_out2,"/modsize_merge_deepsplit/modsize60_merge0.25deepsplit_3.rds"))

