#scRegNet tests

library(Seurat)
library(scRegNet)
library(EnsDb.Hsapiens.v86)
library(fs)
library(dplyr)

Target_name<-"new_URN_int"
Attempt_Number1<-1
Attempt_Number2<-1

datain<-(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/"))

seurat_obj<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"_seur.rds"))


#read in the seurat object with the WT module list
temp_seur<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/",Target_name,"WGCNA_final_WTonly.rds"))
wt_mods<-GetModules(temp_seur)

#Read in the seurat object with teh PIC module list
temp_seur<-readRDS(paste0(datain,"PIC_Only_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
pic_mods<-GetModules(temp_seur)

temp_seur<-NULL

#Read in the seurat object with teh PIC module list
temp_seur<-readRDS(paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/",Target_name,"WGCNA_final_PIC_WT.rds"))
j_mods<-GetModules(temp_seur)

temp_seur<-NULL

###use all gens from teh seurat, take it before you subset for highly variable genes!
#only do pearson correlation

#How Should I subset my seurat object pull 

re_scale_seurat_obj <- ScaleData(seurat_obj,features=rownames(seurat_obj))

#

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "scale.data"))


##running for all WT modules
dataout<-paste0(datain,"Isolate_WGCNA/Attempt_",Attempt_Number1,"/scRegNet/")
dir_create(paste0(dataout,"Individual_Module_tfs/"),recurse=TRUE)

#initialize the summary list and the neg and pos associated tfs
all_assoc<-list()
pos_assoc<-data.frame()
neg_assoc<-data.frame()
#iterating over all modules
for ( i in setdiff(levels(factor(wt_mods$module)),"grey")){
  
  mod_list<-wt_mods[wt_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  
  # a message to keep track on what number module you're on
  print(paste0(which(setdiff(levels(factor(wt_mods$module)),"grey")==i)," out of ", length(setdiff(levels(factor(wt_mods$module)),"grey")==i)))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs/",i,".rds"))
  
  all_assoc<-append(all_assoc,assoc_tfs)
  
  neg_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R<0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)
  pos_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R>0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)
  
  #add the module results to teh complete tf dataframe
  neg_assoc<-rbind(neg_assoc,neg_assoc_tfs)
  pos_assoc<-rbind(pos_assoc,pos_assoc_tfs)
  
}

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))
saveRDS(neg_assoc,paste0(dataout,"negative_associated_tfs.rds"))
saveRDS(pos_assoc,paste0(dataout,"positice_associated_tfs.rds"))



####REPEATING FOR PIC SUBMODULES
#SUBSETTING OUT THE LPS TREATED CELLS AND MOST OF THE CONTROL ONES, LEAVING ONLY THE WT
re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment!="L")

#re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment=="P"|Background=="C")

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "scale.data"))



##running for all PIC modules
dataout<-paste0(datain,"PIC_Only_WGCNA/Attempt_",Attempt_Number2,"/scRegNet/")
dir_create(paste0(dataout,"Individual_Module_tfs/"),recurse=TRUE)

#initialize the summary list and the neg and pos associated tfs
all_assoc<-list()
pos_assoc<-data.frame()
neg_assoc<-data.frame()
#iterating over all modules
for ( i in setdiff(levels(factor(pic_mods$module)),"grey")){
  
  mod_list<-pic_mods[pic_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  # a message to keep track on what number module you're on
  print(paste0(which(setdiff(levels(factor(pic_mods$module)),"grey")==i)," out of ", length(setdiff(levels(factor(pic_mods$module)),"grey")==i)))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs/",i,".rds"))
  
  
  all_assoc<-rbind(all_assoc,assoc_tfs)
  

}

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))
saveRDS(neg_assoc,paste0(dataout,"negative_associated_tfs.rds"))
saveRDS(pos_assoc,paste0(dataout,"positice_associated_tfs.rds"))

neg_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R<0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)
pos_assoc_tfs<-assoc_tfs[assoc_tfs$Pearson_R>0,] %>% arrange(Pearson_p_value) %>% slice_head(n=20)

#add the module results to teh complete tf dataframe
neg_assoc<-rbind(neg_assoc,neg_assoc_tfs)
pos_assoc<-rbind(pos_assoc,pos_assoc_tfs)



all_assoc<-readRDS(paste0(dataout,"All_associated_tfs.rds"))
neg_assoc<-readRDS(paste0(dataout,"negative_associated_tfs.rds"))
pos_assoc<-readRDS(paste0(dataout,"positice_associated_tfs.rds"))


# 3. Create a grouping factor to split the indices evenly
groups <- cut(seq_along(all_assoc), breaks = 12, labels = FALSE)

# 4. Split the list into X sub-lists
split_lists <- split(all_assoc, groups)

all_tfs<-map(split_lists,bind_cols) |> bind_rows()

#SUMMARIZE THE TF RESULTS
#FILTER OUT ANY TFS WITH A PEARSON Q VALUE HIGHER THAN 0.05
count_summ<-dplyr::filter(all_tfs,Pearson_q_value < 0.05) %>%  
  #ADDA COLUMN STATING WHETHER THE PEARSON CORRELATION IS NEGATIVE OR POSITIVE
  mutate(Direction=case_when( Pearson_R >0 ~ "Positive", Pearson_R <0 ~ "Negative" )) %>%
  #ARRANGE BY THE Q VALUE AS WELL AS TEH ABSOLUTE PEARSON CORRELATION, TO KEEP THE MOST CORRELATED NEGATIVE AND POSITIVE TFS
  arrange(Pearson_q_value,abs(Pearson_R))%>%
  #SEPARATE TEH DF BY MODULE AND DIRECTION AND TAKE TEH 50 HIGHEST PERFORMING TFS
  group_by(module, Direction) %>% 
  slice_head(n=50) %>%
  dplyr::count(TF_symbol) 




neg_unique<-neg_assoc %>% group_by(module) %>% dplyr::count(TF_symbol) %>% arrange(desc(n))
pos_unique<-pos_assoc %>% group_by(module) %>% dplyr::count(TF_symbol) %>% arrange(desc(n))

neg_unique<-neg_unique[!is.na(neg_unique$module),]
pos_unique<-pos_unique[!is.na(pos_unique$module),]






#############J PIC ONLY ######


####REPEATING FOR PIC SUBMODULES
#SUBSETTING OUT THE LPS TREATED CELLS AND MOST OF THE CONTROL ONES, LEAVING ONLY THE WT
re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment!="L")

#re_scale_seurat_obj<-subset(re_scale_seurat_obj,subset= Treatment=="P"|Background=="C")

# all genes from the assay
expressed_genes <- rownames(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "data"))

# a scaled expression matrix
m <- as.matrix(GetAssayData(re_scale_seurat_obj, assay = "RNA", layer = "scale.data"))



##running for all PIC modules
dataout<-paste0(datain,"Alt_WGCNA/Attempt_",Attempt_Number2,"/scRegNet/")
dir_create(paste0(dataout,"Individual_Module_tfs/"),recurse=TRUE)

#initialize the summary list and the neg and pos associated tfs
all_assoc<-list()
pos_assoc<-data.frame()
neg_assoc<-data.frame()
#iterating over all modules
for ( i in setdiff(levels(factor(j_mods$module)),"grey")){
  
  mod_list<-j_mods[j_mods$module==i,1]
  
  #get rid of any NA values inside of the module gene list
  mod_list<-mod_list[!is.na(mod_list)]
  
  mech_tfs<-mechanistic(mod_list,expressed_genes,organism="human",gene_id="symbol")
  
  
  #change qvalue depending on teh lowest Q value found in the mech tf list
  q_value=0.005
  if (min(mech_tfs$q_value)>q_value){
    q_value=min(mech_tfs$q_value)
  }
  
  # a message to keep track on what number module you're on
  print(paste0(which(setdiff(levels(factor(j_mods$module)),"grey")==i)," out of ", length(setdiff(levels(factor(j_mods$module)),"grey")==i)))
  assoc_tfs<-associate_with_expression(m,mech_tfs,organism="human",gene_id="symbol",method=c("pearson"),q_value = q_value)
  
  assoc_tfs$module<-i
  
  saveRDS(assoc_tfs,paste0(dataout,"Individual_Module_tfs/",i,".rds"))
  
  
  all_assoc<-rbind(all_assoc,assoc_tfs)
  
  
}


all_assoc<-all_assoc %>% dplyr::filter(Pearson_p_value<0.05) %>% group_by(module) %>% mutate(Direction= case_when(Pearson_R<0 ~ "Down",Pearson_R>0 ~ "Up")) %>% arrange(desc(Pearson_R),Pearson_p_value)

saveRDS(all_assoc,paste0(dataout,"All_associated_tfs.rds"))


summ<-all_assoc %>% arrange(desc(abs(Pearson_R)),Pearson_p_value) %>% group_by(module,Direction) %>% distinct(TF_symbol,.keep_all=TRUE) %>% arrange(desc(Pearson_R),Pearson_p_value)


all_assoc




  