####Scrip to make seurat objects out of DRAGEN single cell analysis output files
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(tibble)
library(Matrix)
library(stringr)

#Change datapaht and the gsub path to reflect directory where the DRAGEN outputs in question are stored

datapath<-"/data/scRNA/HMC3_ZKV/DRAGEN/Realign_again/Output_ds.d6e2c4c6825d46fda615ccfc230f0d78"
#lists all paths tot he dircetories stored under DRAGEN outputs
filenames<-list.dirs(datapath, recursive=FALSE)
#Mkaes a lsit of all biosamples directories in teh DRAGEN outputs. In this case taking the biosample suffix ZSC as the basis of identifying them
filenames<-gsub("/data/scRNA/HMC3_ZKV/DRAGEN/Realign_again/Output_ds.d6e2c4c6825d46fda615ccfc230f0d78/","",filenames[grep("ZSC",filenames)])

obj_list<-list()
#Loop to make a Seurat object from each alignemnt dataset directory
#Also makes a list of the Seurat object for downstream filtering and joining

dir.create(paste0(datapath,"/Seurat_Out/Raw_Obj"),recursive=TRUE)
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

###########Make matrix for expression of TRANSCDS my trangenic gene by sample
transframe<-data.frame(row_id=1:4)
rownames(transframe)<-c("counts","cells","GAPDHcount","Norm")
for (i in seq_along(filt_obj_list)){
  seur<-filt_obj_list[[i]]
  samp<-seur$Sample[[1]]
  transexp<-FetchData(object=seur,vars=c("GAPDH","TRANSCDS"),layer="counts")
  counts<-sum(transexp$TRANSCDS)
  cells<-sum(transexp$TRANSCDS>0)
  hkg<-sum(transexp$GAPDH)
  normTRANS<-counts/hkg
  transframe[[samp]]<-c(counts,cells,hkg,normTRANS)
}

#####UNIFORM CELL COUNT and UNIFORM READ NUMBER merged Seurat objects
####UCC and URN
##subsetting each filtered seurat subject to teh top 512 cells by feature count
UCClist<-list()
for (i in seq_along(filt_obj_list)){
  seurat_obj<-filt_obj_list[[i]]
  #Get the 512 top cells 
  top_cells<-seurat_obj@meta.data %>% arrange(desc(nCount_RNA)) %>% slice(1:512)%>% rownames()
  seurat_top512 <- subset(seurat_obj, cells = top_cells)
  UCClist<-append(UCClist,seurat_top512)
}

UCC_seur<-merge(x=UCClist[[1]],y=UCClist[2:24])
UCC_seur<-JoinLayers(UCC_seur)
#Normalize using SCTransform and save RDs into merged objects directory
UCC_seur<-SCTransform(UCC_seur, vars.to.regress="percent.mt", verbose =FALSE)
saveRDS(UCC_seur,paste0(datapath,"/Seurat_Out/MergedObjects/URN_seur.rds"))

#Subsetting all seurat objects to only include cells with RNA counts above 10,000 and then emrging them into a single Seurat object

URNlist<-list()

for (i in seq_along(filt_obj_list)){
  seurat_obj<-filt_obj_list[[i]]
  #Get any cells above 10,000 feature count 
  seurat_counthresh <- subset(seurat_obj, subset = nCount_RNA>=10000)
  URNlist<-append(URNlist,seurat_counthresh) 
  print(seurat_obj$Sample[[1]])
  print(ncol(seurat_counthresh))
}

URN_seur<-merge(x=URNlist[[1]],y=URNlist[2:24])
URN_seur<-JoinLayers(URN_seur)
URN_seur<-SCTransform(URN_seur, vars.to.regress="percent.mt", verbose =FALSE)
saveRDS(URN_seur,paste0(datapath,"/Seurat_Out/MergedObjects/URN_seur.rds"))
