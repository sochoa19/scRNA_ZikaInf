#RUnning WGCNA on Pseudobulked WT samples from the ZSC dataset
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(SummarizedExperiment)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

datapath<-"/data/scRNA/HMC3_ZSC/Seurat_OUT/"

Target_name<-"UCC_Integrated"

seurat_obj<-readRDS(paste0(datapath,Target_name,"/", Target_name,"_final.rds"))

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "pseudobulk"
)
length(GetWGCNAGenes(seurat_obj))

# get the counts matrix and the meta-data
X <- GetAssayData(seurat_obj, layer='counts')
meta <- seurat_obj@meta.data

# create a pseudo-bulk SummarizedExperiment object
se <- AggregatePseudobulk(
  X, meta, 
  replicate_col = "Sample", 
  group_col = "Background",
  assay_name = 'counts'
)

# normalize the pseudobulk SummarizedExperiment
se <- NormalizeCounts(
  se, 
  method = 'VST',
  assay_name = 'counts'
)

# Set the pseudobulk matrix using the SummarizedExperiment objectseurat_obj <- SetDatExpr(
seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = se, layer = 'VST'
)

####
# subset the WT samples only for pseudobulking
seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = se[,colData(se)$Background == 'C'], 
  layer = 'VST'
)
##

# select the soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)


# construct the co-expression network and identify gene modules
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name='pseudobulk', 
  overwrite_tom=TRUE,
  mergeCutHeight=0.15,
  soft_power = 20
)

#plotting and saving the dendogram
pdf(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/PseudobulkWGCNA/Figures/Pseodubulk_Dendogram.pdf"),width=12, height=9)
PlotDendrogram(seurat_obj, main='Pseudobulk WGCNA Dendrogram')
dev.off()

# compute the MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

# get MEs from seurat object
MEs <- GetMEs(seurat_obj)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add MEs to Seurat meta-data for plotting:
meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)

#Make a dot plot correlating Module Eigengenes against all backgrounds across a single treatment
treat.to.plot<-"P"
p<- DotPlot(subset(x=seurat_obj,subset=Treatment==substr(treat.to.plot,1,1)), features=mods, group.by = 'Background')

# flip the x/y axes, rotate the axis labels, add titles, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  ggtitle(paste0("Correlations in ", treat.to.plot, " treated samples"))

#output plot
p

#save the plot into figures
ggsave(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/PseudobulkWGCNA/Figures/",treat.to.plot,"Dotplot.png"),height=1800,width=5000,units= "px", dpi= 300)



# compute the co-expression network umap 
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=10,
  min_dist=0.4,
  spread=3,
  supervised=TRUE,
  target_weight=0.3
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) + 
  umap_theme() 

# add the module names to the plot by taking the mean coordinates
centroid_df <- umap_df %>% 
  dplyr::group_by(module) %>%
  dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p <- p + geom_label(
  data = centroid_df, 
  label=as.character(centroid_df$module), 
  fontface='bold', size=2) + 
  theme(panel.background = element_rect(fill='black'))

p


