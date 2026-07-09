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

Target_name<-"UCC_Integrated"

seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/",Target_name,"_final.rds"))


# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed, k nearest neighbors parameter and network type for reproducibility
set.seed(12345)
k.parameter<-25
net.type<-"signed"

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)


###Create wild type seurat object
wt_seurat_obj<-subset(seurat_obj, subset= Background=="C")

#Set uo subsetted object for WGCNA
wt_seurat_obj <- SetupForWGCNA(
  wt_seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "WT_isolate" # the name of the hdWGCNA experiment
)

#Creating metacells using the samples as the biosamples
#Would it be better to construct them out of hte same background. so that LPS, PIC ad CTRL cells can comingle into the same emtacell group? In case some cells in the treated samples remain in a homoeostatic state
wt_seurat_obj <- MetacellsByGroups(
  seurat_obj = wt_seurat_obj,
  group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = k.parameter, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Sample' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
wt_seurat_obj <- NormalizeMetacells(wt_seurat_obj)


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
  soft_power=7,
  overwrite_tom=TRUE
  
)
PlotDendrogram(wt_seurat_obj, main='WT isolate hdWGCNA Dendrogram')


p1 <- DimPlot(wt_seurat_obj, group.by='Treatment') +
  umap_theme() +
  ggtitle('WT metacell Isolate') 
  

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

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)


##correlelogram to check module correlation
ModuleCorrelogram(seurat_obj,features="MEs")


##Correlating modules to samples using a dot plot
#getting the MEs
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)



#Make a dot plot correlating Module Eigengenes against all backgrounds across a single treatment
treat.to.plot<-"PIC"
p<- DotPlot(subset(x=seurat_obj,subset=Treatment==substr(treat.to.plot,1,1)), features=mods, group.by = 'Background')

# flip the x/y axes, rotate the axis labels, add titles, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  ggtitle(paste0("WT module correlations in ", treat.to.plot, " samples"))

#save the plot into figures
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Figures/",treat.to.plot,"Dotplot.png"))
# plot output
p


for (i in unique(seurat_obj$Background)){
  z<-DotPlot(subset(x=seurat_obj,subset=Background==i), features=mods, group.by = 'Treatment')
  z <- z +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    ggtitle(paste0("WT module correlation in ", i, " background"))
  ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Figures/",i,"Dotplot.png"))
  
}







#Run this to make sure you get kME values from teh full dataset not the wt_isolate
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = '', group_name = 'INH'
)

seurat_query <- ModuleExprScore(
  seurat_query,
  method='UCell'
)



##Correlating modules to samples using a dot plot
#getting the MEs
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
