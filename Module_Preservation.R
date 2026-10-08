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
library(dplyr)

Target_name<-"new_URN_int_sct"

#Provide attempt number. This is to keep straight what module results are tied to which metacell construction and hdWGCNA parameter
Attempt_Number<-"2"


datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Attempt_",Attempt_Number)
dataout<-paste0(datain,"/Module_Preservation")

dir_create(dataout)


##Reading in seurat objects and parameter list used in ATTEMPT_NUMBER
seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Seurat_Analysis/",Target_name,"_final.rds"))
wt_seurat_obj<-readRDS(paste0(datain,"/",Target_name,"WGCNA_final_WTonly.rds"))

parameter_list<-readRDS(paste0(datain,"/Parameter_list.rds"))
k.parameter<-as.numeric(parameter_list[1])
max_shared<-as.numeric(parameter_list[2])
net.type<-parameter_list[3]

#Using the native module preservation protocol from WGCNA
##############MODULE PRESERVATION#####

for (i in setdiff(levels(factor(seurat_obj$Background)),"12")){
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
    group.by ="Sample",
    k = k.parameter,
    max_shared = max_shared,
    reduction = 'pca',
    ident.group = 'Sample'
  )
  
  pres_seurat_obj <- NormalizeMetacells(pres_seurat_obj)
  
  
  wt_seurat_obj <- SetDatExpr(
    wt_seurat_obj,
    assay = 'SCT', # using SCT assay
    layer = 'data', # using normalized data
    use_metacells=TRUE
  )
  
  
  pres_seurat_obj <- SetDatExpr(
    pres_seurat_obj,
    assay = 'SCT', # using RNA assay
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
  saveRDS(mod_pres,paste0(dataout,"/Z_DF_",Backtest,".rds"))
  
}

#Plotting all the results in a single plot

wt_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/new_URN_int_sctWGCNA_final_WTonly_projected.rds")
#Pulling al the mp files and summarizing the results in a single graph
pres_summ<-list()
cor_summ<-list()


for (i in unique(wt_seurat_obj$Background)){
  df<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/Module_Preservation/Z_DF_",i,".rds"))
  
  #figure out what to do with median rank. maybe
  pres_summ[[i]]<-df$Zsummary.pres
  cor_summ[[i]]<-df$Z.cor.cor
}

pres_summ<-as.data.frame(pres_summ)
rownames(pres_summ)<-rownames(df)
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


#change background values to the right thing
plot_df$Background<-replace_values(plot_df$Background, 
                                   "C" ~ "WT-PIC",
                                   "1"~"NS1-PIC",
                                   "2"~"NS2a-PIC",
                                   "3"~"NS2b-PIC",
                                   "4"~"NS3-PIC",
                                   "7"~"NS4a-PIC",
                                   "12"~"NS5-PIC",
                                   "J"~"GFP-PIC",)


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
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8, color = "black"),
    legend.title=element_blank(),
    axis.title.x= element_text(size=11,color="black"),
    axis.title.y= element_text(size=11,color="black"),
    plot.title = element_text(size=11,color="black",hjust=0.5,face="plain"),
    axis.text = element_text(size=10)
    
  ) +
  labs(
    x = "Module",
    y = "Zsummary",
    title="Module preservation from WT to ZIKV NS Backgrounds"
  )


