#Script to look at module preservation across backgrounds of our ZSC scRNA dataset

#takes in the subset wt seurat object with the complete module data included, taken from the wt_isolate_WGCNA script

#takes in the complete Seurat object from Seurat Analysis. This will be subset by backgrounds and then the wt modules will be tested across all backgrounds

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
Target_name<-"UCC_Integrated"

#Provide attempt number. This is to keep straight what module results are tied to which metacell construction and hdWGCNA parameter
Attempt_Number<-"2"


datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Attempt_",Attempt_Number)
dataout<-paste0(dataout,"/Module_Preservation")

dir_create(dataout)


##Reading in seurat objects and parameter list used in ATTEMPT_NUMBER
seurat_obj<-readRDS(paste0(datain,"/",Target_name,"WGCNA_final.rds"))
wt_seurat_obj<-readRDS(past(datain,"/",Target_name,"WGCNA_final_WTonly.rds"))

parameter_list<-readRDS(paste0(datain,"/Parameter_list.rds"))
k.parameter<-as.numeric(parameter_list[1])
max_shared<-as.numeric(parameter_list[2])
net.type<-parameter_list[3]

#Using the native module preservation protocol from WGCNA
##############MODULE PRESERVATION#####

for (i in unique(seurat_obj$Background)){
  #Iterate the module preservation test over all the 
  Backtest<-i
  
  #subset the background to check for module preservation
  pres_seurat_obj<-subset(seurat_obj, subset= Background==Backtest)
  
  
  #Set up subsetted object for WGCNA
  pres_seurat_obj <- SetupForWGCNA(
    pres_seurat_obj,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "Back_isolate" # the name of the hdWGCNA experiment
  )
  
  # construct metacells for query dataset:
  pres_seurat_obj <- MetacellsByGroups(
    seurat_obj = pres_seurat_obj,
    group.by ="Background",
    k = k.parameter,
    max_shared = max_shared,
    reduction = 'pca',
    ident.group = 'Background'
  )
  
  pres_seurat_obj <- NormalizeMetacells(pres_seurat_obj)
  
  
  wt_seurat_obj <- SetDatExpr(
    wt_seurat_obj,
    assay = 'RNA', # using RNA assay
    layer = 'data', # using normalized data
    use_metacells=TRUE
  )
  
  
  pres_seurat_obj <- SetDatExpr(
    pres_seurat_obj,
    assay = 'RNA', # using RNA assay
    layer = 'data', # using normalized data
    use_metacells=TRUE
  )
  
  #Run Module Preservation. 
  pres_seurat_obj <- ModulePreservation(
    pres_seurat_obj,
    seurat_ref = wt_seurat_obj,
    name=paste0("WT to ",Backtest),
    verbose=3,
    n_permutations=150 # can be lowered
  )
  
  # get the module preservation table
  mod_pres <- GetModulePreservation(pres_seurat_obj, paste0("WT to ",Backtest))$Z
  obs_df <- GetModulePreservation(pres_seurat_obj, paste0("WT to ",Backtest))$obs
  saveRDS(obs_df,paste0(dataout,"/obs_DF_",Backtest,".rds"))
  saveRDS(mods_pres,paste0(dataout,"/Z_DF_",Backtest,".rds"))
  
}



#Using the NetRep protocol 

#