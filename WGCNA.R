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

seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/UCC_Integrated_final.rds")


# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Sample' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


#Run hdWGCNA on only the WT cells. THis will allow for the discovery of
#expected module response to LPS and PIC
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("ZSCCC","ZSCCL","ZSCCP"), # the name of the group of interest in the group.by column
  group.by='Sample', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  networktype='signed',
  tom_name = 'WT' # name of the topological overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='WT hdWGCNA Dendrogram')


# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Sample"
)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Sample', group_name = c('ZSCCC','ZSCCP','ZSCCL')
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "WT"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)



#Metacells and module eigengene expression of these
metacell_obj <- GetMetacellObject(seurat_obj)

metacell_obj <- NormalizeData(metacell_obj)
metacell_obj <- FindVariableFeatures(metacell_obj)
metacell_obj <- ScaleData(metacell_obj)
metacell_obj <- RunPCA(metacell_obj)


DimPlotMetacells(seurat_obj, group.by='Sample') + umap_theme() + ggtitle("Sample")


