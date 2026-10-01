#KLF hunting
#looking to find the effects of NS2b on the KLF targets

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(scatterpie)
library(ggplot2)
library(ggvenn)
library(pheatmap)
library(ggalluvial)
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

#general help
library(fs)
library(purrr)
library(readr)
library(tidyr)
library(tibble)
library(Matrix)
library(stringr)

#1. Finding KLF targets in the brown yellow and blue modules
all_tfs<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/ALL_PIC_WGCNA/Attempt_2/scRegNet/All_associated_tfs.rds")
all_tfs<-all_tfs %>% dplyr::filter(Pearson_q_value < 0.05)
#susbetting the scREGNET results just to the module desired
yellow_tf<-all_tfs %>% dplyr::filter(module=="yellow")
#making a subset of the module TFs to only include TFs in the KLF family
yellow_klf<-yellow_tf[grepl("KLF",yellow_tf$TF_symbol),]
yellow_klf_target<-unique(unlist(unlist(str_split(yellow_klf$TF_targets_symbol , " "))))

yellow_non<-yellow_tf[!grepl("KLF",yellow_tf$TF_symbol),]
yellow_non_target<-unique(unlist(unlist(str_split(yellow_non$TF_targets_symbol , " "))))

#remove any shared targets form non target. There is a complete overlap. between non targets and klf targets. 
yellow_non_target<-yellow_non_target[!(yellow_non_target %in% yellow_klf_target)]
#total number of genes in module is 272, yell targets is 99, non targets 147, the rest are not targetted


blue_tf<-all_tfs %>% dplyr::filter(module=="blue")

#making a subset of the module TFs to only include TFs in the KLF family
blue_klf<-blue_tf[grepl("KLF",blue_tf$TF_symbol),]
blue_klf_target<-unique(unlist(unlist(str_split(blue_klf$TF_targets_symbol , " "))))

blue_non<-blue_tf[!grepl("KLF",blue_tf$TF_symbol),]
blue_non_target<-unique(unlist(unlist(str_split(blue_non$TF_targets_symbol , " "))))

#remove any shared targets form non target. There is a complete overlap. between non targets and klf targets. 
blue_non_target<-blue_non_target[!(blue_non_target %in% blue_klf_target)]
#total number of genes in module is 1276, klf targets is 646, non targets 405,  the rest are not targeted.



brown_tf<-all_tfs %>% dplyr::filter(module=="brown")

brown_klf<-brown_tf[grepl("KLF",brown_tf$TF_symbol),]
brown_klf_target<-unique(unlist(unlist(str_split(brown_klf$TF_targets_symbol , " "))))

brown_non<-brown_tf[!grepl("KLF",brown_tf$TF_symbol),]
brown_non_target<-unique(unlist(unlist(str_split(brown_non$TF_targets_symbol , " "))))

#remove any shared targets form non target. between non targets and klf targets. 
brown_non_target<-brown_non_target[!(brown_non_target %in% brown_klf_target)]
#total number of genes in module is 319, klf targets is 81, non targets 190, of those 7 are unique to KLF, teh rest are not targeted.
# Unique genes "ATP8B3" "BEX2"   "CDC25C" "LDOC1"  "MLF2"   "OSGEP"  "WDR45B"


#calculate metagenes and append to expression data
seurat<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_seur.rds")

cells<-subset(seurat, Background %in% c("J","3","1") & Treatment %in% c("P","C") )
meta<-cells@meta.data
expr<-GetAssayData(cells, assay="RNA", layer = "data")

genes <- intersect(brown_klf_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$brown_klf<-pc1



genes <- intersect(brown_non_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$brown_non<-pc1


genes <- intersect(blue_klf_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$blue_klf<-pc1



genes <- intersect(blue_non_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$blue_non<-pc1



genes <- intersect(yellow_klf_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$yellow_klf<-pc1



genes <- intersect(yellow_non_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]

# Mean expression of the gene set
mean_expr <- rowMeans(x)

# Orient PC1 consistently
if (cor(pc1, mean_expr) < 0) {
  pc1 <- -pc1
}

meta$yellow_non<-pc1


#adding all the target genes together and getting PCA 1 and two of those to comapre PCa
yellow_all_target<-append(yellow_klf_target, yellow_non_target)
blue_all_target<-append(blue_klf_target, blue_non_target)
brown_all_target<-append(brown_klf_target, brown_non_target)



genes <- intersect(yellow_all_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]
# Mean expression of the gene set
mean_expr <- rowMeans(x)
pc1_loadings <- pca$rotation[, 2]
top_pc1_genes <- sort(abs(pc1_loadings), decreasing = TRUE)
yellow_movers<-names(head(top_pc1_genes, 20)) # View top 20 driving genes


meta$yellow_pc1<-pc1
meta$yellow_pc2<-pc2


genes <- intersect(blue_all_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]
# Mean expression of the gene set
mean_expr <- rowMeans(x)
pc1_loadings <- pca$rotation[, 2]
top_pc1_genes <- sort(abs(pc1_loadings), decreasing = TRUE)
blue_movers<-names(head(top_pc1_genes, 20)) # View top 20 driving genes

meta$blue_pc1<-pc1
meta$blue_pc2<-pc2


genes <- intersect(brown_all_target, rownames(expr))
x <- t(as.matrix(expr[genes, , drop = FALSE]))
# PCA
pca <- prcomp(
  x,
  center = TRUE,
  scale. = TRUE
)
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]
# Mean expression of the gene set
mean_expr <- rowMeans(x)

pc1_loadings <- pca$rotation[, 2]
top_pc1_genes <- sort(abs(pc1_loadings), decreasing = TRUE)
brown_movers<-names(head(top_pc1_genes, 20)) # View top 20 driving genes


meta$brown_pc1<-pc1
meta$brown_pc2<-pc2

movers<-list(blue=blue_movers,yellow=yellow_movers,brown=brown_movers)

saveRDS(meta,"/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/TF_metagene_analysis/metagene_meta.rds")
saveRDS(movers,"/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/TF_metagene_analysis/movers.rds")
meta<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/TF_metagene_analysis/metagene_meta.rds")








ggplot(data=meta %>% filter(Treatment=="P"), aes(x=brown_pc1,y=brown_pc2,color=Sample))+
  geom_point(size=0.5)+
  labs(title="Brown PCA")

ggplot(data=meta %>% filter(Treatment=="P"), aes(x=blue_pc1,y=blue_pc2,color=Sample))+
  geom_point(size=0.5)+
  labs(title="Blue PCA")


ggplot(data=meta %>% filter(Treatment=="P"), aes(x=yellow_pc1,y=yellow_pc2,color=Sample))+
  geom_point(size=0.5)+
  labs(title="Yellow PCA")



meta_long<-pivot_longer(meta,cols=brown_klf:yellow_non,names_to = "gene_list")

ggplot(meta_long %>% dplyr::filter(Treatment=="P"),aes(x=gene_list, y=value,color = Sample))+
  geom_violin()




