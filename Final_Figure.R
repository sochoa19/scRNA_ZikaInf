
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
#GO ANNOTATION
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)


#Recreating Figures
dataout<-"/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Paper_Figures"
dir_create(dataout)
base=2


#Figure 1 B ####
Seur_target<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/new_URN_int_sct_final.rds")
Seur_target$Background<-replace_values(Seur_target$Background, 
                                       "C" ~ "WT",
                                       "1"~"NS1",
                                       "2"~"NS2a",
                                       "3"~"NS2b",
                                       "4"~"NS3",
                                       "7"~"NS4a",
                                       "12"~"NS5",
                                       "J"~"GFP",)


DimPlot(Seur_target, reduction = "umap",group.by="Background",label=FALSE,stroke.size=base/15)+
  labs(title="")+
  
  theme_cowplot( )+
  guides(color = guide_legend(override.aes = list(size = base/2)))+
  theme(legend.text = element_text(size = base*2, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.title.x= element_text(size=base*4,color="black"),
        legend.key.size = unit(0.2,"line"),
        legend.box.spacing = unit(0, "pt"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_blank(),
        axis.text = element_text(size=base*2)
  )
ggsave(paste0(dataout,"/fig1B.png"),width=base,height=base,unit="in",dpi=300)


#Figure 1 C ####
Seur_target<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/new_URN_int_sct_final.rds")

DimPlot(Seur_target, reduction = "umap",group.by="Treatment",label=FALSE,stroke.size=base/15)+
  labs(title="")+
  
  scale_color_discrete(labels =c("C"= "CTRL", "L"="LPS", "P"="PIC"),palette = c("C"="#619CFF","L"="#00BA38","P"="#F8766D") )+
  theme_cowplot( )+
  guides(color = guide_legend(override.aes = list(size = base/2)))+
  theme(legend.text = element_text(size = base*2, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.title.x= element_text(size=base*4,color="black"),
        legend.key.size = unit(0.2,"line"),
        legend.box.spacing = unit(0, "pt"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_blank(),
        axis.text = element_text(size=base*2)
  )
ggsave(paste0(dataout,"/fig1C.png"),width=base,height=base,unit="in",dpi=300)




#Figure 1 D ####

#volcanoplots for new_URN_int with teh fixed titles and the same list of ISGs, IFIT

vsP_DE<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/vsP_DE.rds")
WT<-vsP_DE[["C"]]

#Capping pval adjusted at 10^-302 to allow for better plotting
WT[WT$p_val_adj==0,"p_val_adj"]<-1*10^-302

#Pulling out the significant genes with teh lowest pvalue, that are known PIC associated genes, ISG IFIT, and oas
WT_PIC<-grepl("ISG|IFIT|OAS|CCL",WT[WT$expression != "NS", ][order(WT[WT$expression != "NS", ]$p_val_adj), "gene"])
WT_PIC_GENES<-head(WT[WT_PIC,],10)

WT$isg<-ifelse(grepl("ISG|IFIT|OAS|CCL|IFI",WT$gene), "Yes", "No")

sigWT<-WT[WT$expression != "NS", ]
highlight_data<-sigWT[grepl("ISG|IFIT|OAS|CCL|IFI",sigWT$gene),]

volcano<-ggplot(WT, aes(x = avg_log2FC, y = -log10(p_val_adj) , color = expression)) +
  geom_point( size = base/30) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red"),
                     labels =c("Downregulated", "Not Significant", "Upregulated") ) +
  # Layer 2: Plot only the selected points (bright and prominent)
  geom_point(data = highlight_data, 
             color="orange",
             size = base/10) +

  # Threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black",linewidth=base/7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",linewidth=base/7) +
  expand_limits(y = -log10(320))+
  
  

  scale_x_continuous(breaks = seq(-4,8, by = 2), limits=c(-4,8)) +

  scale_y_continuous(breaks = seq(0, 320, by = 100), limits=c(0,330)) +
  # Formatting
  labs(
    title = "WT-CTRL vs WT-PIC",
    x = "Log2 Fold Change",
    y = "-Log10 p-value-adj",
    color = ""
  ) +
  theme_cowplot( )+
  theme(legend.text = element_text(size = base, color = "black"),
        legend.position = "none",
        axis.title.x= element_text(size=base*4,color="black"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_text(size=base*5,color="black",hjust=1,face="plain"),
        axis.text = element_text(size=base*2)
        )
volcano
ggsave(paste0(dataout,"/fig1D.png"), width = base, height = base, units="in",dpi = 300)



#Figure 1 E ####

#volcanoplots for new_URN_int with teh fixed titles and the same list of ISGs, IFIT

WT<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/vsL_DE.rds")
WT$gene<-rownames(WT)

#Capping pval adjusted at 10^-302 to allow for better plotting
WT[WT$p_val_adj==0,"p_val_adj"]<-1*10^-302

#Pulling out the significant genes with teh lowest pvalue, that are known PIC associated genes, ISG IFIT, and oas
WT_LPS<-grepl("ISG|IFIT|OAS|CCL",WT[WT$expression != "NS", ][order(WT[WT$expression != "NS", ]$p_val_adj), "gene"])
WT_LPS_GENES<-head(WT[WT_PIC,],10)


sigWT<-WT[WT$expression != "NS", ]
highlight_data<-sigWT[grepl("ISG|IFIT|OAS|CCL|IFI",sigWT$gene),]

volcano<-ggplot(WT, aes(x = avg_log2FC, y = -log10(p_val_adj) , color = expression)) +
  geom_point( size = base/30) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red"),
                     labels =c("Downregulated", "Not Significant", "Upregulated") ) +
  # Layer 2: Plot only the selected points (bright and prominent)
  geom_point(data = highlight_data, 
             color="orange",
             size = base/10) +
  
  # Threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black",linewidth=base/7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",linewidth=base/7) +
  expand_limits(y = -log10(320))+
  
  
  
  scale_x_continuous(breaks = seq(-4,8, by = 2), limits=c(-4,8)) +
  
  scale_y_continuous(breaks = seq(0, 320, by = 100), limits=c(0,330)) +
  # Formatting
  labs(
    title = "WT-CTRL against WT-LPS",
    x = "Log2 Fold Change",
    y = "-Log10 p-value-adj",
    color = ""
  ) +
  theme_cowplot( )+
  theme(legend.text = element_text(size = base, color = "black"),
        legend.position = "none",
        axis.title.x= element_text(size=base*4,color="black"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_text(size=base*5,color="black",hjust=1,face="plain"),
        axis.text = element_text(size=base*2)
  )
volcano
ggsave(paste0(dataout,"/fig1E.png"), width = base, height = base, units="in",dpi = 300)



#Figure1F

#Plotting the cluster I keep getting different cluster each time I try to do this. Might ahve to redo at some point
cluster_plot<-Seur_target@meta.data %>% group_by(Treatment) %>% dplyr::count(seurat_clusters)
ggplot(cluster_plot,aes(x=Treatment,y=n,fill=seurat_clusters))+
  geom_bar(position="stack",stat = "identity")+
  theme_cowplot( )+
  labs(
   fill=" Cluster",
    x = "Treatment",
    y = "Number of Cells"
  ) +
  guides(fill = guide_legend(ncol = 2))+
  theme(legend.text = element_text(size = base*2, color = "black"),
        legend.key.size = unit(base/20, "in"),
        legend.title = element_text(size=base*4,color="black"),
        axis.title.x= element_text(size=base*4,color="black"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_text(size=base*5,color="black",hjust=1,face="plain"),
        axis.text = element_text(size=base*2)
  )
ggsave(paste0(dataout,"/fig1F.png"), width = base, height = base, units="in",dpi = 300)






#Figure 2 A####
wt_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/new_URN_int_sctWGCNA_final_WTonly.rds")
png(paste0(dataout,"/fig2A.png"), width = base*2, height = base, units="in", res = 300)
PlotDendrogram(wt_seurat_obj, main='WT isolate hdWGCNA Dendrogram',return_plot=TRUE)
dev.off()



#Figure 2B ####
WT_WGCNA<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/Isolate_WGCNA/Attempt_1/new_URN_intWGCNA_final_WTonly_projected.rds")
MEs<-WT_WGCNA@meta.data %>% arrange(Treatment) %>% dplyr::filter(Background=="C")
MEs$Cell<-seq_along(rownames(MEs))

MEs_long <- MEs %>%
  pivot_longer(
    cols = c(turquoise, green, brown),  # your modules
    names_to = "Module",
    values_to = "Eigengene"
  )

#adding y offset values so that I can plota ll three module tracks in teh same graph
module_offsets <- c(
  turquoise = 0,
  green      = 50,
  brown     = 100
)

MEs_long <- MEs_long %>%
  mutate(
    Offset = module_offsets[Module],
    Y = Eigengene + Offset
  )

treatment_ranges <- MEs_long %>%
  group_by(Treatment) %>%
  summarise(
    xmin = min(Cell) - 0.5,
    xmax = max(Cell) + 0.5,
    .groups = "drop"
  )

module_axes <- MEs_long %>%
  distinct(Module, Offset) %>%
  mutate(
    ymin = Offset - 20,
    ymax = Offset + 20
  )

ggplot(MEs_long, aes(
  x = Cell,
  y = Y,
  color = Treatment
)) +
  geom_segment(
    aes(
      xend = Cell,
      y = Offset,
      yend = Y
    ),
    linewidth = 0.5
  ) +
  
  #Creating the colored treatemnt rectangle labels
  geom_rect(
    data = treatment_ranges,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = -25
    ),     
    inherit.aes = FALSE,
    fill = c("#F8766D","#00BA38","#619CFF"),
    color = "transparent",
    linewidth = 0.5
  ) +
  #adding the label text to the middle of the rectangles
  geom_text(data=treatment_ranges, aes(x = (xmin + xmax) / 2, 
                y = -31, 
                label = c("CTRL","LPS","PIC")),
            inherit.aes=FALSE,
            size = 3)+
  
  # Module baselines
  geom_segment(
    aes(y = Offset,
        yend=Offset,
        x=0,
        xend=max(Cell)+100),
    linetype = "dotted",
    linewidth = 0.4,
    color = "black"
    
  ) +
  
  # Module labels
  geom_text(
    data = module_axes,
    aes(
      x = min(MEs_long$Cell) - 20,
      y = Offset,
      label = c("Turquoise","Green","Brown")
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3
  ) +
  
  #adding right hand side axis
  geom_segment(
    data = module_axes,
    aes(
      x = max(MEs_long$Cell) + 100,
      xend = max(MEs_long$Cell) + 100,
      y = ymin,
      yend = ymax
    ),
    inherit.aes = FALSE
  ) +
  
  #adding right hand axis labels, min
  geom_text(
    data = module_axes,
    aes(
      x = max(MEs_long$Cell) + 150,
      y = ymin,
      label = "-20"
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust=0,
    size = 3
  ) +
  #adding right hand axis labels, min
  geom_text(
    data = module_axes,
    aes(
      x = max(MEs_long$Cell) + 150,
      y = ymax,
      label = "20"
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust=1,
    size = 3
  ) +
  
  
  coord_cartesian(clip = "off")+
  #facet_wrap(~ turquoise, scales = "free_y") +
  theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y= element_blank(),
    plot.title = element_text(size=11,color="black",hjust=0.5,face="plain"),
    axis.text = element_text(size=10),
    # Legend at bottom
    legend.position = "none",
    
    # Make room for module labels and right-side axes
    plot.margin = margin(
      5, 10, 5, 30)
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL,
    title = "Module Activity in WT cells"
  )
ggsave(paste0(dataout,"/fig2B.png"), width = base*2, height = base, units="in", dpi = 300)








#Figure 2C####
#GO terms

#######RUNNING GO on the modules and on the genes with kme above 0.5
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

#read in teh seurat object
wt_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/new_URN_int_sctWGCNA_final_WTonly.rds")

#get modules from wild type isolate reference seuarat object
modules<-GetModules(wt_seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']


#Compiles entrez ids for the universe of possible gene
universe.symb <- bitr(
  modules$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#takes any genes still in ENSEMBL format and converts them into entrez ids. HIGH ERROR RATE LOOK INTO IT!
universe.ensembl <- bitr(
  modules[grepl("ENSG",modules$gene_name),1],
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


univ.combined<-append(universe.symb$ENTREZID, universe.ensembl$ENTREZID)

#


ontol<-"BP" #This can be changed depending on what ontology database you want to use


currcolor<-"turquoise" #sets the name of the module to be put through GO enrichment analysis
genelist<-modules[modules$color==currcolor,1] #populates genelist with all the genes that make up this module

#Convert all the gene symbols into entrez ids. 
gene.symb <- bitr(
  genelist,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#Checks to see if any of the ENSEMBL ids in the module can be converted. If not skip conversion 
if (length(grep(paste(genelist[grepl("ENSG",genelist)],collapse="|"),keys(org.Hs.eg.db, keytype = "ENSEMBL")))==0
){
  gene.combined<-gene.symb$ENTREZID
  
}else{
  gene.ensembl<- bitr(
    genelist[grepl("ENSG",genelist)],
    fromType= "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
    
  )
  
  gene.combined <-append(gene.symb$ENTREZID,gene.ensembl$ENTREZID) #combine symbol->entrez and ensembl->entrez
  
}


qcutoff=0.05 #sets the qvalue cutoff to be used in GO enrichment 0.05 is standard
#GO enrichment 
ego <- enrichGO(
  gene = gene.combined,
  universe = univ.combined,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = ontol, #This can be changed depending on what ontology database you want to use
  readable = TRUE,
  qvalueCutoff = qcutoff
  
)

#create a dotplot of the 10 most significant enriched terms for that module
dotplot(ego, showCategory = 10)+
  labs(title="BP GO Enrichment of Turquoise Module",
       x="Ratio of GO Genes in Module")+
  ggtitle(paste0(ontol,"GO enrichment of ", currcolor, " module"))+
  scale_size(range = c(base/2, base*2))+
  theme_cowplot()+
  theme(legend.text = element_text(size = base*2, color = "black"),
        legend.key.size = unit(base/20, "in"),
        legend.title = element_text(size=base*4,color="black"),
        axis.title.x= element_text(size=base*4,color="black"),
        axis.line = element_line(colour = 'black', size = base/7),
        axis.ticks = element_line(size=base/7),
        axis.title.y= element_text(size=base*4,color="black"),
        plot.title = element_text(size=base*5,color="black",hjust=1,face="plain"),
        axis.text = element_text(size=base*2)
  )

ggsave(paste0(dataout,"/fig2C.png"),height=base,width=2*base,units= "in", dpi= 300)



#Figure 2D ####
#Venn diagram
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggvenn)
wt_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/new_URN_int_sctWGCNA_final_WTonly.rds")
turq_mod<-GetModules(wt_seurat_obj) %>% dplyr::filter(module=="turquoise")

#subsetting the module genes
turq_genes<-turq_mod[,1]

#retrieving the WT-PIC DE genes

vsP_DE<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/vsP_DE.rds")
WT<-vsP_DE[["C"]]

#Capping pval adjusted at 10^-302 to allow for better plotting
WT[WT$p_val_adj==0,"p_val_adj"]<-1*10^-302

#subsetting upregulated WT-PIC genes
DE_genes<-WT[WT$expression=="Up","gene"]


#Retreiving GO genes
genes_df <- select(org.Hs.eg.db, 
                   keys = "GO:0051607", 
                   columns = c("SYMBOL"), 
                   keytype = "GOALL")
GO_symbols<-unique(genes_df$SYMBOL)
GO_genes<-GO_symbols[GO_symbols %in% rownames(wt_seurat_obj@assays$SCT$data)]

universe<-unique(rownames(wt_seurat_obj@assays$SCT$data))

venn_lists <- list(
  "Turquoise Module"=turq_genes,
  "DE"=DE_genes,
  "GO"=GO_genes
)

ggvenn(
  venn_lists, 
  fill_color = c("#27dbe2", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, 
  set_name_size = 4,
  show_percentage=FALSE
)



#PIC object origin,

#READ IN BOTH MODULE GENE LISTS
datain<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/")
WT_obj<-readRDS(paste0(datain,"Isolate_WGCNA/Attempt_1/new_URN_intWGCNA_final_WTonly.rds"))
PIC_obj<-readRDS(paste0(datain,"Alt_WGCNA/Attempt_1/new_URN_intWGCNA_final_PIC_WT.rds"))


wtmods<-GetModules(WT_obj)
picmods<-GetModules(PIC_obj)
sum(wtmods$gene_name %in%picmods$gene_name)
merged_mods<-merge(wtmods[,1:2],picmods[,1:2], by="gene_name")
colnames(merged_mods)<-c("gene_name","WT_Modules","PIC_Modules")
merged_mods_grey<-merged_mods[! merged_mods$PIC_Modules=="grey",]
alluvial_data_grey<-merged_mods_grey %>% dplyr::count(WT_Modules,PIC_Modules) 
merged_mods<-merged_mods[!(merged_mods$WT_Modules=="grey" | merged_mods$PIC_Modules=="grey"),]
alluvial_data<-merged_mods %>% dplyr::count(WT_Modules,PIC_Modules) 


#make a vector of all the module colors and assign their own color to them. 
module_colors<-as.character(unique(c(merged_mods$PIC_Modules,merged_mods$WT_Modules)))
names(module_colors)<-unique(c(merged_mods$PIC_Modules,merged_mods$WT_Modules))

module_colors<-gsub("black","grey40",module_colors)



ggplot(data=alluvial_data_grey,
       aes(axis1=WT_Modules,axis2=PIC_Modules,y=n))+
  geom_alluvium(aes(fill=WT_Modules),width=1/5)+
  geom_stratum(
    aes(fill = after_stat(stratum)),
    color = "transparent",
    width=1/5
  ) +
  
  geom_text(stat="stratum",aes(label = after_stat(stratum)),
            color="black",
            size=3,
            width=1/5)+
  scale_x_discrete(limits = c("WT_Modules", "PIC_Modules"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(
    values = module_colors
  ) +
  scale_x_discrete(
    limits = c("WT_All_Modules", "GFP_PIC_Modules"),
    expand = c(0.15, 0.05)
  ) +
  labs(
    y="Number of genes",
    title="Gene Module Identity"
  )+
  theme_cowplot()+
  theme(legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size=11,color="black"),
        axis.title.x= element_text(size=11,color="black"),
        axis.title.y= element_text(size=11,color="black"),
        plot.title = element_text(size=11,color="black",hjust=0.5,face="plain"),
        axis.text = element_text(size=10),
        legend.position = "none",
  )

ggsave(paste0(dataout,"/fig2D.png"),height=8,width=6,dpi=300)



#Figure 2E
#Module preservation

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




#Figure 3A####
pic_all_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/All_PIC_WGCNA/Attempt_1/new_URN_intWGCNA_final_PIC_WT.rds")


MEs <- GetMEs(pic_all_obj, harmonized=FALSE)
modules <- GetModules(pic_all_obj)
pic_all_obj<-AddMetaData(pic_all_obj,MEs)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
meta<-pic_all_obj@meta.data
j_only<-subset(pic_all_obj,Sample=="ZSCJP")
#j_only<-subset(pic_all_obj,Background=="J")


distances<-list()
for (i in mods){
  #finding the mean and standard deviation form the 
  mod_sd<-sd(j_only@meta.data[,i])
  mod_mean<-mean(j_only@meta.data[,i])
  
 # distances[[i]]<-(as.data.frame(meta %>% group_by(Background) %>% summarise(mean(.data[[i]])))[,2]-mod_mean)/mod_sd
  distances[[i]]<-(as.data.frame(meta %>% group_by(Sample) %>% summarise(mean(.data[[i]])))[,2]-mod_mean)/mod_sd
  
}
plot_values<-as.data.frame(distances)
#plot_values$origin<-levels(factor(pic_all_obj$Background))
plot_values$origin<-levels(factor(pic_all_obj$Sample))

plot_long<-plot_values %>% pivot_longer(cols=-origin, names_to = "Modules", values_to = "Eigengene")

plot_long$Treatment<-str_sub(gsub("ZSC","",plot_long$origin),-1,-1) 
plot_long$Background<-str_sub(gsub("ZSC","",plot_long$origin),1,-2) 

#change background values to the right thing
plot_long$origin<-replace_values(plot_long$origin, 
                                       "C" ~ "WT-PIC",
                                       "1"~"NS1-PIC",
                                       "2"~"NS2a-PIC",
                                       "3"~"NS2b-PIC",
                                       "4"~"NS3-PIC",
                                       "7"~"NS4a-PIC",
                                       "12"~"NS5-PIC",
                                       "J"~"GFP-PIC",)

#remove the the GFP sample
#plot_long<-plot_long %>% dplyr::filter(origin!="GFP-PIC")
plot_long<-plot_long %>% dplyr::filter(Background!="J")


ggplot(plot_long, aes(y = Modules, x = Eigengene)) +
  geom_vline(
    aes(xintercept = 0),
    linetype = "dotted",
    linewidth = 0.8
  ) +
  geom_vline(
    aes(xintercept = -2),
    linetype = "dotted",
    linewidth = 0.3,
    color="red"
  ) +
  geom_vline(
    aes(xintercept = 2),
    linetype = "dotted",
    linewidth = 0.3,
    color="red"
  ) +
  geom_point(
    aes(color = Background,shape=Treatment),
    #aes(color = origin)
    size = 3,
    position = position_jitter(width = 0, height = 0.08)
  ) +
  theme_cowplot() +
  labs(
    x = "Z score of eigengene expression, relative to GFP-PIC",
    y = "Module",
    color = "Background"
  )+
  theme(legend.text = element_text(size = 8, color = "black"),
        legend.title=element_blank(),
        axis.title.x= element_text(size=11,color="black"),
        axis.title.y= element_text(size=11,color="black"),
        plot.title = element_text(size=11,color="black",hjust=0.5,face="plain"),
        axis.text = element_text(size=10))

ggsave(paste0(dataout,"/fig3A.png"),height=5, width=6, dpi=300)




#Figure 3B ####

p_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/Alt_WGCNA/Attempt_1/new_URN_intWGCNA_final_PIC_WT_projected_allPIC.rds")
#Pulling al the mp files and summarizing the results in a single graph
pres_summ<-list()
median_summ<-list()
cor_summ<-list()


for (i in unique(p_seurat_obj$Background)){
  df<-readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int/Alt_WGCNA/Attempt_1/Module_preservation/",i,"mod_pres.rds"))
  pres <- df$preservation$Z[[1]][[2]]
  #figure out what to do with median rank. maybe
  pres$MedianRank<-df$quality$observed$ref.ref$inColumnsAlsoPresentIn.query$medianRank.qual
  pres_summ[[i]]<-pres$Zsummary.pres
  median_summ[[i]]<-pres$MedianRank
  cor_summ[[i]]<-df$preservation$observed$ref.ref$inColumnsAlsoPresentIn.query$cor.cor
}

pres_summ<-as.data.frame(pres_summ)
rownames(pres_summ)<-rownames(pres)
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
    y = "Zsummary"
  )

ggsave(paste0(dataout,"/fig3b.png"),width=6, height=5, dpi=300)



###VENN DIAGRAM DE GO and TURQMODULE#######
library(AnnotationDbi)
library(org.Hs.eg.db)

# Retrieve Ensembl IDs for a specific GO term
gene_list <- select(org.Hs.eg.db, 
                    keys = "GO:0051607", 
                    keytype = "GOALL", 
                    columns = c("SYMBOL", "ENSEMBL", "ENTREZID"))


DE_data<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Seurat_Analysis/vsP_DE.rds")
C_DE<-dplyr::filter(DE_data[["C"]], expression=="Up")
WT_seur<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/Isolate_WGCNA/Attempt_2/new_URN_int_sctWGCNA_final_WTonly.rds")
turq_mod<-dplyr::filter(GetModules(WT_seur), module=="turquoise")

mod_venn<-turq_mod$gene_name
de_venn<-C_DE$gene
go_venn<-unique(gene_list$SYMBOL)

#Have to find which of the GO terms are even present in our analysis. Subsetting them for comparison
gene_universe<-DE_data[["C"]]$gene
go_venn_filtered<-go_venn[go_venn %in% gene_universe]

library(ggvenn)
venn<-list(Turquoise_Module=mod_venn, PIC_DE=de_venn, Antiviral_GO=go_venn_filtered)
ggvenn(venn, fill_color = c("red", "blue", "green"))





####Network highlighting####
pic_seurat_obj<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/All_PIC_WGCNA/Attempt_3/new_URN_int_sctWGCNA_final_PIC_WT.rds")

#Pull the topological overlap matrix form teh seuart object
TOM<-GetTOM(pic_seurat_obj)

#Pull the list of gene module membership from teh seuart obejct
pic_mods<-GetModules(pic_seurat_obj)

#Select which module you want to plot
mod_name<-"red"
#pull the gene list in teh selected module
current_mod<-dplyr::filter(pic_mods,module==mod_name)

#Subset the TOM with the genes inside the selected module
currentTOM<-TOM[current_mod$gene_name,current_mod$gene_name]

#Zero out the upper half to the TOM, which is symetrical to the bottom half. across a diagonal
currentTOM[upper.tri(currentTOM)] <- NA


# cast the network from wide to long format
cur_network <- currentTOM %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene1 = Var1, gene2 = Var2, weight=value) %>%
  subset(!is.na(weight))


# keep this network before subsetting
cur_network_full <- cur_network 

# keep the top 10% of edges 
edge_percent <- 0.1
cur_network <- cur_network_full %>% 
  dplyr::slice_max(
    order_by = weight, 
    n = round(nrow(cur_network)*edge_percent)
  )

# make the graph object with tidygraph
graph <- cur_network %>% 
  igraph::graph_from_data_frame() %>%
  tidygraph::as_tbl_graph(directed=FALSE) %>% 
  tidygraph::activate(nodes)


# get the top 25 hub genes for the selected module
hub_genes<-head(current_mod %>% arrange(desc(get(paste0("kME_",mod_name)))),n=10)
hub_genes<-hub_genes$gene_name                
V(graph)$hub <- ifelse(V(graph)$name %in% hub_genes, V(graph)$name, "")


#Change the color of the nodes depending if they are Differntially expressed in ZSC3P or not
DE_df<-readRDS("/data/scRNA/HMC3_ZSC/Seurat_OUT/new_URN_int_sct/All_PIC_WGCNA/Attempt_3/vsP_DE.rds")
curr_background<-dplyr::filter(DE_df[["3"]],expression=="Down")[,"gene"]
test<-DE_df[["3"]]


graph <- graph %>%
  activate(nodes) %>%
  mutate(
    in_gene_set = name %in% curr_background
  )


# make the plot with gggraph
p <- ggraph(graph) + 
  geom_edge_link(aes(alpha=weight)) + 
  geom_node_point(
    aes(color = name %in% curr_background)
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "red")
  ) +
  geom_node_label(aes(label=hub), repel=TRUE, max.overlaps=Inf, fontface='italic')+
  theme_cowplot() +
  labs(title="Red Module Network, NS2b-PIC DE")+
  theme(
    legend.text = element_text(size = 8, color = "black"),
    legend.title=element_blank(),
    axis.title.x= element_blank(),
    axis.title.y= element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size=11,color="black",hjust=0.5,face="plain"),
    axis.text = element_blank()
    
  ) +
  labs(
    x = "Module",
    y = "Zsummary"
  )


p



