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

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "WT_based" # the name of the hdWGCNA experiment
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Background"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = k.parameter, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Background' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


#Run hdWGCNA on only the WT cells. THis will allow for the discovery of
#expected module response to LPS and PIC
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("C"), # the name of the group of interest in the group.by column
  group.by='Background', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = net.type # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

#Consturct the WGCNA Network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  networktype='signed',
  tom_name = 'WTbybackground' , # name of the topological overlap matrix written to disk
  #tom_outdir= paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA"),
  soft_power=7,
  overwrite_tom=TRUE
  
)

pdf(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA/Figures/WTbyBackground_Dendogram.pdf"),width=12, height=9)
PlotDendrogram(seurat_obj, main='WT hdWGCNA Dendrogram')
dev.off()

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=NULL
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
  ggtitle(paste0("Correlations in ", treat.to.plot, " treated samples"))

#save the plot into figures
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA/Figures/",treat.to.plot,"Dotplot.png"))
# plot output
p

for (i in unique(seurat_obj$Background)){
  z<-DotPlot(subset(x=seurat_obj,subset=Background==i), features=mods, group.by = 'Treatment')
  z <- z +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    ggtitle(paste0("Correlations in ", i, " background samples"))
  ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA/Figures/",i,"Dotplot.png"))
  
}
setwd("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_Integrated/WGCNA")
TOM<-(load("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_Integrated/WGCNA/WTbybackground_TOM.rda"))
TOM<-as.matrix(consTomDS)


saveRDS(seurat_obj,paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA/Seurat_WGCNA.rds"))
wg_seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/WGCNA/Seurat_WGCNA.rds"))

# compute eigengene-based connectivity (kME):
wg_seurat_obj <- ModuleConnectivity(
  wg_seurat_obj,
  group.by = 'Background', group_name = 'C',
  
)


# plot genes ranked by kME for each module
p <- PlotKMEs(wg_seurat_obj, ncol=5)


# get hub genes
hub_df <- GetHubGenes(wg_seurat_obj, n_hubs = 10)

head(hub_df)




#Metacells and module eigengene expression of these
metacell_obj <- GetMetacellObject(seurat_obj)

metacell_obj <- NormalizeData(metacell_obj)
metacell_obj <- FindVariableFeatures(metacell_obj)
metacell_obj <- ScaleData(metacell_obj)
metacell_obj <- RunPCA(metacell_obj)


DimPlotMetacells(seurat_obj, group.by='Sample') + umap_theme() + ggtitle("Sample")


