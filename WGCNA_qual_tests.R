#WGCNA PARAMETER - MODULE QUALITY TEST
#a streamlined pipeline to test Module quality across different WGCNA parameters
#Will return the module preservation across pseudobulk and bulk datasets
#Will also produce GO terms enrichment and a short list of hub genes 

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

#general help
library(fs)
library(dplyr)
library(stringr)
library(readr)


#######RUNNING GO on the modules and on the genes with kme above 0.5
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

theme_set(theme_cowplot())
#set for reproducibility
set.seed(12345)
enableWGCNAThreads(nThreads = 4)


#Name of the starting Seurat object UCC_Integrated or URN_int
Target_name<-"UCC_Integrated"
#Metacell parameters chosen through metacell analysis script
k_parameter<-60
net.type<-"signed"
max_shared<-15

#WGCNA Parameters to be tested


#Create a datapath string for incoming seurat objects
data_in<-"/data/scRNA/HMC3_ZSC/Seurat_OUT"
#Create a datapath string for outgoing objects
data_out<-paste0(data_in,"/",Target_name,"/WGCNA_parameter")
#Create out path if not already existing
dir_create(data_out)
dir_create(paste0(data_out,"/modsize_merge_deepsplit"))
dir_create(paste0(data_out,"/modsize_merge"))
dir_create(paste0(data_out,"/modsize_deepsplit"))



#Read in the Seurat Object
seurat_obj<-readRDS(paste0(data_in,"/",Target_name,"/Seurat_Analysis/",Target_name,"_final.rds"))

#Read in the TSV bulk RNA seq object
bulk_target<-"HMC3_INF"

bulk_obj<-readRDS(paste0("/data/bulkRNA/",bulk_target,"/DESeq2_results/bulk_df.rds"))

###########MOD SIZE and MERGE HEIGHT#####

modulesizerange <- c(20,30,40,50)
mergerange <- c(0.15,0.20,0.25)
deep_split_range<-4
for(i in modulesizerange){
  for (j in mergerange){
    #Collect all parameters into a single vector to feed into the collected function
    param_list<-c(k_parameter,net.type,max_shared,i,j,deep_split_range)
    combname<-paste0("modsize",i,"_merge",j,".rds")
    saveRDS(object=Complete_WGCNA(seurat_obj,bulk_obj,param_list),file =paste0(data_out,"/",combname) )
  }
}


#summarizing results from the parameter search
param_outs<-list()
for (i in list.files(paste0(data_out,"/modsize_merge"))){
  df<-as.data.frame(read_csv(paste0(data_out,"/modsize_merge/",i)))
  pass<-sum(df$Zsumm>=2)
  fly<-sum(df$Zsumm>=10)
  passgene<-sum(df[df$Zsumm>2,4])/sum(df$size)
  modnum<-nrow(df)
  
  param_outs[[gsub(".csv","",i)]]<-c(modnum,
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

ggsave(paste0(data_out,"/modsize_merge/parameter_summary.png"))

#########MOD SIZE AND DEEPSPLIT#####

#witht he results from teh last parameter search we're sticking to 0.25 as the merge cut height
#Now we are iterating over module size and deepsplit
modulesizerange <- c(40,50,60,70)
mergerange <- 0.25
deep_split_range<-c(1,2,3,4)
for(i in modulesizerange){
  for (j in deep_split_range){
    #Collect all parameters into a single vector to feed into the collected function
    param_list<-c(k_parameter,net.type,max_shared,i,mergerange,j)
    combname<-paste0("modsize",i,"_deepsplit",j,".rds")
    saveRDS(object=Complete_WGCNA(seurat_obj,bulk_obj,param_list),file =paste0(data_out,"/modsize_deepsplit/",combname) )
  }
}

#summarizing results from the parameter search
param_outs<-list()
for (i in list.files(paste0(data_out,"/modsize_merge_deepsplit"))){
  df<-as.data.frame(readRDS(paste0(data_out,"/modsize_merge_deepsplit/",i)))
  pass<-sum(df$Zsumm>=2)
  fly<-sum(df$Zsumm>=10)
  passgene<-sum(df[df$Zsumm>2,4])/sum(df$size)
  modnum<-nrow(df)
  
  param_outs[[gsub(".rds","",i)]]<-c(modnum,
                                     as.numeric((pass/modnum)*100),
                                     as.numeric((fly/modnum)*100),
                                     as.numeric(str_sub(i,start=8,end=9)),
                                     as.numeric(str_sub(i,start=-5,end=-5)),
                                     as.numeric(passgene))
}
paramdf<-as.data.frame(param_outs)
rownames(paramdf)<-c("size","Percent_passing","percent_strongly_passing","modsize","deep_split","Pass_Gene")

# Transpose so each row is one parameter combination
plot_df <- as.data.frame(t(paramdf))

# Plot
ggplot(plot_df,
       aes(x = deep_split,
           y = modsize,
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
  theme_bw()

ggsave(paste0(data_out,"/modsize_deepsplit/parameter_summary.png"))

#######MOD,MERGE,DEEPSPLIT#####
dir_create(paste0(data_out,"/modsize_merge_deepsplit"))

modulesizerange <- c(50)
mergerange <- c(0.20,0.25)
deep_split_range<-c(3,4)
for(i in modulesizerange){
  for (j in mergerange){
    for (h in deep_split_range){
      #Collect all parameters into a single vector to feed into the collected function
      param_list<-c(k_parameter,net.type,max_shared,i,j,h)
      combname<-paste0("modsize",i,"_merge",j,"deepsplit_",h,".rds")
      saveRDS(object=Complete_WGCNA(seurat_obj,bulk_obj,param_list),file =paste0(data_out,"/modsize_merge_deepsplit/",combname) )
    }
  }
}


library(ggcube)

#summarizing results from the parameter search
param_outs<-list()
for (i in list.files(paste0(data_out,"/modsize_merge_deepsplit"))){
  df<-as.data.frame(readRDS(paste0(data_out,"/modsize_merge_deepsplit/",i)))
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

ggsave(paste0(data_out,"/modsize_merge_deepsplit/parameter_summary.png"))


####Read in specific parameter rds#######
modsize=60
merge=0.2
deepsplit=4
check<-readRDS(paste0(data_out,"/modsize_merge_deepsplit/modsize",modsize,"_merge",merge,"deepsplit_",deepsplit,".rds"))
#########violinplots of all zsumms

zsumms<-list()

remove_words<-paste(c("merge","modsize","deepsplit",".rds"), collapse = "|")

####violin plots of zsummary of each module exlcuding gray and gold which seemed bugged right now
for (i in list.files(paste0(data_out,"/modsize_merge_deepsplit/"),pattern="\\.rds$")){
  df<-readRDS(paste0(data_out,"/modsize_merge_deepsplit/",i))
  zsumms[[str_remove_all(i, remove_words)]]<-df[!rownames(df) %in% c("grey","gold"),6]
}

plot_df <- stack(zsumms)
colnames(plot_df) <- c("Value", "Group")
ggplot(plot_df, aes(x = Group, y = Value)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  stat_summary(
    fun = median,
    geom= "point",
    shape=3,
    size = 2,
    color = "red"
  ) +
  theme_bw()
ggsave(paste0(data_out,"/modsize_merge_deepsplit/parameter_violin.png"),width=10)



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
  

  

  #get modules from wild type isolate reference seurat object
  modules<-GetModules(wt_seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  #Create the module DF that will be populated with all module stats and details
  module.df<-data.frame(name=mods)
  rownames(module.df)<-mods
  
  #Calculate Module Eigengenes in the WT subset seurat object
  wt_seurat_obj<-ModuleEigengenes(wt_seurat_obj)
  
  #Calucalte kme of each module gene
  wt_seurat_obj <- ModuleConnectivity(
    wt_seurat_obj
  )
  
  #save out the hub gnes into a dataframe
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




