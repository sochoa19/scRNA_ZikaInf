#Repeat WGCNA script but subsetting the WT samples entirely before getting metacells
#hdWGCNA analysis taking SeuratObjects 

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
Target_name<-"URN_int"

#Provide attempt number. This is to keep straight what module results are tied to which metacell construction and hdWGCNA parameter
Attempt_Number<-"1"


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

parameter_list<-c(k.parameter,max_shared,net.type,Target_name,Attempt_Number)
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
  overwrite_tom=TRUE
  
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
saveRDS(seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final.rds"))

saveRDS(wt_seurat_obj,paste0(dataout,"/",Target_name,"WGCNA_final_WTonly.rds"))


wt_seurat_obj <- ScaleMetacells(wt_seurat_obj, features=VariableFeatures(wt_seurat_obj))
wt_seurat_obj <- RunPCAMetacells(wt_seurat_obj, features=VariableFeatures(wt_seurat_obj))
wt_seurat_obj <- RunUMAPMetacells(wt_seurat_obj, reduction='pca', dims=1:15)


p <- DimPlotMetacells(wt_seurat_obj, group.by='Background') + umap_theme() + ggtitle("Sample")

p



