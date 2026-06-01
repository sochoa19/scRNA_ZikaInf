####Scrip to make seurat objects out of DRAGEN single cell analysis output files
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(tibble)
library(Matrix)
library(stringr)

datapath<-"/data/scRNA/HMC3_ZKV/DRAGEN/Realign"
filenames<-list.dirs(datapath, recursive=FALSE)
filenames<-gsub("/data/scRNA/HMC3_ZKV/DRAGEN/Realign/","",filenames[grep("ZSC",filenames)])

obj_list<-list()
#Loop to make a Suerat object from each alignemnt dataset directory
#Also makes a list of the Seurat object for downstream filtering and joining
for (i in filenames){
  specpath<-paste0(datapath,"/",i,"/",i)
  expression_matrix <- ReadMtx(
    mtx = paste0(specpath,".scRNA.matrix.mtx.gz"),
    cells = paste0(specpath,".scRNA.barcodes.tsv.gz"),
    features = paste0(specpath,".scRNA.features.tsv.gz")
  )
  seurat_obj <- CreateSeuratObject(
    counts = expression_matrix,
    project = "ZSC",
    min.cells = 3,
    min.features = 200
  )
  saveRDS(seurat_obj,paste0(datapath,"/Seurat_Out/Raw_Obj/",i,".rds"))
  obj_list<-append(obj_list,seurat_obj)
}

#Filtering to cull any cells with a mt percent above 5% and features <200 >15000
#also creates plots of percent mt vs nCount and nFeature vs nCount. both filtered obj
#and plot are saved into an directory named filtered

filt_obj_list<-list()
for (i in seq_along(obj_list)){
  seur<-obj_list[[i]]
  #Add metadata column of mitochondrial gene percentage
  seur[["percent.mt"]]<-PercentageFeatureSet(seur,pattern="^MT-")
  #add Treatment C=Control, L=LPS, P=PIC
  seur[["Treatment"]]<-str_sub(filenames[i],-1)
  #add ZSC plasmid that was being transduced
  seur[["Background"]]<-str_sub(filenames[i],-2,-2)
  seur[["Sample"]]<-filenames[i]
  ###save plots to a png
  plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  ggsave(paste0(datapath,"/Seurat_Out/Filtered/plot",filenames[i],".png"))
  #
  seur <- subset(seur, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 5)
  filt_obj_list<-append(filt_obj_list,seur)
  saveRDS(seur,paste0(datapath,"/Seurat_Out/Filtered/",filenames[i],".rds"))
  
}


merge(x=filt_obj_list[1],y=filt_obj_list[2:24])
