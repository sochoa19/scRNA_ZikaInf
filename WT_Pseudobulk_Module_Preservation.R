#Script to compare WT scRNA seq modules to pseudobulked data

#All needed libraries

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

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 4)

#Change these values to read in different versions of the seurat objects
Attempt_number<-"2"
Target_name<-"UCC_Integrated"

#read in wt seurat object post WGCNA 
wtdatain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Attempt_",Attempt_number,"/")
wt_seurat_obj<-readRDS(paste0(wtdatain,Target_name,"WGCNA_final_WTonly.rds"))

#read in paramters used to create the WT modules, to re use in the 
parameter_list<-readRDS(paste0(wtdatain,"Parameter_list.rds"))
k_parameter<-as.numeric(parameter_list[1])
max_shared<-as.numeric(parameter_list[2])
net.type<-parameter_list[3]

#read in full seurat object before WGCNA
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/",Target_name,"_final.rds")
seurat_obj<-readRDS(datain)

#Set up seurat object as pseudobulk object and WGCNA
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

expr <- assay(se, "VST")

# transpose to samples x genes
expr <- t(expr)

seurat_obj <- SetDatExpr(
  seurat_obj,
  mat = expr
)

#Project wt modules onto Pseudobulkdata
seurat_obj <- ProjectModules(
  seurat_obj = seurat_obj,
  seurat_ref = wt_seurat_obj,
  wgcna_name = "WT_isolate",
  wgcna_name_proj="projected",
  assay="RNA" # assay for query dataset
)

# set expression matrix for reference dataset
wt_seurat_obj <- SetDatExpr(
  wt_seurat_obj,
  group_name = "C",
  group.by = "Background"
)


#Run Module prservation of wt modules on complete pseudobulked dataset.. 
seurat_obj <- ModulePreservation(
  seurat_obj,
  seurat_ref = wt_seurat_obj,
  name="WT to Pseudobulk",
  verbose=3,
  n_permutations=150 # can be lowered
)

# get the module preservation table
mod_pres <- GetModulePreservation(seurat_obj, "WT to Pseudobulk")$Z
obs_df <- GetModulePreservation(seurat_obj, "WT to Pseudobulk")$obs

dir_create(paste0(wtdatain,"WT_to_Pseudo/Figures"),recurse = TRUE)

saveRDS(mod_pres,paste0(wtdatain,"WT_to_Pseudo/mod_pres.rds"))
saveRDS(obs_df,paste0(wtdatain,"WT_to_Pseudo/obs_df.rds"))

#plot summary stats
plot_list <- PlotModulePreservation(
  seurat_obj,
  name="WT to Pseudobulk",
  statistics = "summary"
)

wrap_plots(plot_list, ncol=2)
ggsave(paste0(wtdatain,"WT_to_Pseudo/Figures/Summary.png"),width=8,height=5,bg="white")

#plot ranking stats
plot_list <- PlotModulePreservation(
  seurat_obj,
  name="WT to Pseudobulk",
  statistics = "rank"
)

wrap_plots(plot_list, ncol=2)
ggsave(paste0(wtdatain,"WT_to_Pseudo/Figures/Rank.png"),width=8,height=8)


#plot all stats
plot_list <- PlotModulePreservation(
  seurat_obj,
  name="WT to Pseudobulk",
  statistics = "all",
  plot_labels = FALSE,
  #label_size = 0.4,
  mod_point_size=3
)

wrap_plots(plot_list, ncol=6)
ggsave(paste0(wtdatain,"WT_to_Pseudo/Figures/ALL.png"),width=21,height=14,bg="white")




