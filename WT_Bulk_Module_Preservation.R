#WT to Bulk Module PReservation

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
Attempt_number<-"1"
Target_name<-"URN_int"

#read in wt seurat object post WGCNA 
wtdatain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Isolate_WGCNA/Attempt_",Attempt_number,"/")
wt_seurat_obj<-readRDS(paste0(wtdatain,Target_name,"WGCNA_final_WTonly.rds"))

data_out<-paste0(wtdatain,"WT_to_Bulk")
#Make a directory to save the module preservation outputs
dir_create(data_out)

#read in paramters used to create the WT modules, to re use in the 
parameter_list<-readRDS(paste0(wtdatain,"Parameter_list.rds"))
k_parameter<-as.numeric(parameter_list[1])
max_shared<-as.numeric(parameter_list[2])
net.type<-parameter_list[3]

# set expression matrix for reference dataset
wt_expr <- GetDatExpr(wt_seurat_obj)

#read in the bulk RNA seq dataframe
bulk_target<-"HMC3_INF"
bulk_filtered<-readRDS(paste0("/data/bulkRNA/",bulk_target,"/DESeq2_results/filt_bulk_df.rds"))
bulk<-readRDS(paste0("/data/bulkRNA/",bulk_target,"/DESeq2_results/bulk_df.rds"))


#subsetting only the genes that are shared across datasets
# get the genes that are in common: 
bulk_filtered<-t(bulk_filtered)
genes_bulk_filtered <- colnames(bulk_filtered)
genes_sc <- colnames(wt_expr)
genes_keep <- genes_bulk_filtered[genes_bulk_filtered %in% genes_sc]
bulk_kept_filtered <- bulk_filtered[,genes_keep]

# set up multiExpr:
setLabels <- c("ref", "query")
multiExpr <- list(
  ref = list(data=wt_expr),
  query = list(data=bulk_kept_filtered)
)

# set up the modules
ref_modules <- list(ref = GetModules(wt_seurat_obj)$module)

mp <- WGCNA::modulePreservation(
  multiExpr,
  ref_modules,
  referenceNetworks = 1,
  nPermutations = 150 # set this to whatever number is suitable for you
)

#Saving the module preservation file 
saveRDS(mp,paste0(data_out,"/Module_Preservation.rds"))


pres <- mp$preservation$Z[[1]][[2]]

plot(
  pres$moduleSize,
  pres$Zsummary.pres,
  pch=19,
  xlab="Module size",
  ylab="Zsummary"
)

abline(h=2, lty=2, col="red")
abline(h=10, lty=2, col="blue")

text(
  pres$moduleSize,
  pres$Zsummary.pres,
  labels=rownames(pres),
  pos=3,
  cex=0.8
)


plot(
  pres$moduleSize,
  mp$preservation$observed$ref.ref$inColumnsAlsoPresentIn.query$medianRank.pres,
  pch=19,
  xlab="Module size",
  ylab="medianRank"
)

text(
  pres$moduleSize,
  mp$preservation$observed$ref.ref$inColumnsAlsoPresentIn.query$medianRank.pres,
  labels=rownames(pres),
  pos=3
)


plot.df <- data.frame(
  Module = rownames(pres),
  ModuleSize = pres$moduleSize,
  Zsummary = pres$Zsummary.pres,
  medianRank = mp$preservation$observed$ref.ref$inColumnsAlsoPresentIn.query$medianRank.pres
)

ggplot(plot.df,
       aes(ModuleSize,
           Zsummary,
           color=medianRank,
           )) +
  geom_point(alpha=.8,size=4) +
  scale_colour_gradient(low="blue",high="yellow")+
  geom_text(aes(label=Module),
            color='black',
            nudge_y=1) +
  geom_hline(yintercept=2,colour="red",linetype="dashed")+
  geom_hline(yintercept=10,colour="black",linetype="dashed")
  theme_classic()

  ggsave(paste0(data_out,"/Zsummary.png"))

