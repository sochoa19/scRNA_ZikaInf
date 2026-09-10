#ALTERNATIVE PIC WGCNA

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

Target_name<-"new_URN_int"

# set random seed, k nearest neighbors parameter and network type for reproducibility
set.seed(12345)
k.parameter<-45
max_shared<-15
net.type<-"signed"
merge_height<-0.15
deep_split<-2
mod_size<-40


Attempt_Number<-1
Attempt_Number1<-1
Attempt_Number2<-1




P_parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_Number,merge_height,deep_split,mod_size)





###PIC_ONLY_WGCNA
J_PIC_ONLY_WGCNA<-function(Target_name,Attempt_Number,parameter_list){
  
  datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis")
  dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Alt_WGCNA/","Attempt_",Attempt_Number)
  
  k.parameter<-as.numeric(parameter_list[1])
  max_shared<-as.numeric(parameter_list[2])
  net.type<-parameter_list[3]
  merge_height<-as.numeric(parameter_list[4])
  deep_split<-as.numeric(parameter_list[5])
  mod_size<-as.numeric(parameter_list[6])
  
 
  
  #make sure the data_out folder is still there
  dir_create(dataout)
  
  saveRDS(parameter_list,paste0(dataout,"/Parameter_list.rds"))
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
  
  
  ###Create a GFP type seurat object
  wt_seurat_obj<-subset(seurat_obj, subset= Background=="J")
  
  ##Create a PIC only GFP seurat object
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
    tom_name = 'J_PIC_isolate' , # name of the topological overlap matrix written to disk
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

#READ IN BOTH MODULE GENE LISTS
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/")
WT_obj<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/",Target_name,"WGCNA_final_WTonly.rds"))
PIC_obj<-readRDS(paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))


wtmods<-GetModules(WT_obj)
picmods<-GetModules(PIC_obj)
sum(wtmods$gene_name %in%picmods$gene_name)
merged_mods<-merge(wtmods[,1:2],picmods[,1:2], by="gene_name")
colnames(merged_mods)<-c("gene_name","WT_Modules","PIC_Modules")
merged_mods_grey<-merged_mods[! merged_mods$PIC_Modules=="grey",]
alluvial_data_grey<-merged_mods_grey %>% dplyr::count(WT_Modules,PIC_Modules) 
merged_mods<-merged_mods[!(merged_mods$WT_Modules=="grey" | merged_mods$PIC_Modules=="grey"),]
alluvial_data<-merged_mods %>% dplyr::count(WT_Modules,PIC_Modules) 


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
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Alt_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny.png"),height=10,width=8)



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
ggsave(paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/Figures/ModuleDestiny_wGrey.png"),height=10,width=8)


##scRegNet on the genes transitioning out of turquoise
library(scRegNet)
library(EnsDb.Hsapiens.v86)

#find out which PIC modules are mostly made up of turquoise WT genes
propmods<-merged_mods_grey %>% group_by(PIC_Modules) %>% dplyr::count(WT_Modules) %>% mutate(prop=n/sum(n))

descendant_mods<-levels(factor(propmods[propmods$WT_Modules=="turquoise" & propmods$prop >0.3 , "PIC_Modules"][[1]]))


turq_fate<-merged_mods %>% dplyr::filter(WT_Modules=="turquoise",PIC_Modules %in% descendant_mods) 

seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"_seur.rds"))

dataout<-paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/Module_Fate/")
dir_create(dataout)
#How Should I subset my seurat object pull 

re_scale_seurat_obj <- ScaleData(seurat_obj,features=rownames(seurat_obj))

##SUBSETTING OUT THE LPS TREATED CELLS AND MOST OF THE CONTROL ONES, LEAVING ONLY THE WT
re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment!="L")

#re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment=="P"|Background=="C")

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "scale.data"))


#initialize the summary list and the neg and pos associated tfs
all_assoc<-data.frame()
pos_assoc<-data.frame()
neg_assoc<-data.frame()

for ( i in levels(factor(turq_fate$PIC_Modules))){
  
  mod_list<-turq_fate[turq_fate$PIC_Modules==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  
  # a message to keep track on what number module you're on
  print(paste0(which(levels(factor(turq_fate$PIC_Modules))==i)," out of ", length(levels(factor(turq_fate$PIC_Modules)))))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs_",i,".rds"))
  
  all_assoc<-rbind(all_assoc,assoc_tfs)
  
  
}

all_assoc<-all_assoc %>% dplyr::filter(Pearson_p_value<0.05) %>% group_by(module) %>% arrange(desc(Pearson_R),Pearson_p_value)

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))


####scRegNet on all modules



##scRegNet on the genes transitioning out of turquoise
library(scRegNet)
library(EnsDb.Hsapiens.v86)


seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"_seur.rds"))

pic_mods<-GetModules(PIC_obj)



dataout<-paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/scRegNet/")
dir_create(dataout)
#How Should I subset my seurat object pull 

re_scale_seurat_obj <- ScaleData(seurat_obj,features=rownames(seurat_obj))

##SUBSETTING OUT THE LPS TREATED CELLS AND MOST OF THE CONTROL ONES, LEAVING ONLY THE WT
re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment!="L")

#re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment=="P"|Background=="C")

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "scale.data"))


#initialize the summary list and the neg and pos associated tfs
all_assoc<-data.frame()
pos_assoc<-data.frame()
neg_assoc<-data.frame()

for ( i in setdiff(levels(factor(pic_mods$module)),"grey")){
  
  mod_list<-pic_mods[pic_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  
  # a message to keep track on what number module you're on
  print(paste0(which(levels(factor(pic_mods$module))==i)," out of ", length(levels(factor(pic_mods$module)))))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(mech_tfs,paste0(dataout,"Mechanistic_Individual_Module_tfs_",i,".rds"))
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs_",i,".rds"))
  
  all_assoc<-rbind(all_assoc,assoc_tfs)
  
  
}

all_assoc<-all_assoc %>% dplyr::filter(Pearson_p_value<0.05) %>% group_by(module) %>% arrange(desc(Pearson_R),Pearson_p_value)

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))



PIC_ONLY_WGCNA(Target_name,Attempt_Number,parameter_list)


test<-GetHubGenes(PIC_obj, n_hubs = 20)
modsumm<-test2 %>% dplyr::count(module)
tests<-All_associated_tfs[,c(3,7,8,10)] %>% 
  mutate(Direction=case_when(Pearson_R>0 ~ "Up",Pearson_R<0 ~ "Down")) %>%
  group_by(module,Direction) %>% 
  arrange(desc(abs(Pearson_R)),Pearson_p_value) %>%
  distinct(TF_symbol, .keep_all = TRUE) %>%
  slice_head(n=10)

