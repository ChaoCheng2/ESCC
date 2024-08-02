##FigureS1——related to Figure 1######
setwd(dir = workdir)
cb_pattern <- c("#4dbbd5ff","#E64b35ff","#00a087ff",'#fdbf6f',"#3c5488ff","#f39b7fff","#8491b4ff","#a65628","#91d1c2ff",
                "#7876b1ff",'#fb9a99',"#377eb8","#4daf4a","#984ea3","#ff7f00","#7e6148ff","#f781bf","#999999",'#1f78b4',
                '#b2df8a','#33a02c', '#F4D160','#9EB384','#CDC2AE',"#b09c85ff")

library(ggplot2)
library(dplyr)
library(Seurat)
library(lisi)
library(dplyr)
library(tidyr)
library(magrittr)
library(reshape2)
library(ggpubr)
library(harmony)
library(ggsci)
library(vegan)

######### funtion ##############
process_seurat <-function(obj, dim=30, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  vf_gene = VariableFeatures(obj)
  obj <- ScaleData(obj, features = vf_gene)
  obj <- RunPCA(obj, features = vf_gene, npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "pca", n.components=n.components)
  return(obj)
}

process_seurat_harmony <-function(obj, dim=20, n.components=2,resolution=0.3){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  vf_gene = VariableFeatures(obj)
  obj <- ScaleData(obj, features = vf_gene)
  obj <- RunPCA(obj, features = vf_gene, npcs = 50, verbose = TRUE)
  obj = RunHarmony(obj, "batch", plot_convergence = TRUE)
  obj <- FindNeighbors(obj, reduction = "harmony",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "harmony", n.components=n.components)
  return(obj)
}

# Figure S1 

All_cell_pre_harmony <- readRDS("./Final_RDS/All_cells_before_harmony.rds")
All_cell_post_harmony <- readRDS("./Final_RDS/All_cells_after_harmony.rds")
All_cell_post_harmony$cell.type[is.na(All_cell_post_harmony$cell.type)]='Tumor or epithelial cells'
All_cell_post_harmony@meta.data[All_cell_post_harmony@meta.data[, 'cell.type'] == 'NA',  'cell.type'] = 'Tumor or epithelial cells'
All_cell_pre_harmony$cell.type[is.na(All_cell_pre_harmony$cell.type)]='Tumor or epithelial cells'
All_cell_pre_harmony@meta.data[All_cell_pre_harmony@meta.data[, 'cell.type'] == 'NA',  'cell.type'] = 'Tumor or epithelial cells'

save(All_cell_pre_harmony,All_cell_post_harmony,file = 'Final_RDS/supplement_Figure_S1A_All_cell.Rdata')


load('Final_RDS/supplement_Figure_S1A_All_cell.Rdata')

#before
P1 = DimPlot(All_cell_pre_harmony, reduction = "umap",label = F, group.by = "batch",raster=FALSE) + 
  scale_color_manual(values = c("#E64b35ff","#8491b4ff","#91d1c2ff"))+
  ggtitle('Before batch') +
  NoLegend()

#after
P2 = DimPlot(All_cell_post_harmony, reduction = "umap",label = F, group.by = "batch",raster=FALSE) + 
  scale_color_manual(values = c("#E64b35ff","#8491b4ff","#91d1c2ff"))+
  ggtitle('After batch') +
  NoLegend()



##Tcell batch
T_cell <- readRDS("./Final_RDS/T_cell.rds")
T_cell_before = process_seurat(T_cell, resolution = 0.3)
save(T_cell,T_cell_before,file = 'Final_RDS/supplement_Figure_S1A_Tcell.Rdata')

load('Final_RDS/supplement_Figure_S1A_Tcell.Rdata')

#before-Tcell
P3 = DimPlot(T_cell_before, reduction = "umap",label = F, group.by = "batch", raster=FALSE) + 
  scale_color_manual(values = c("#E64b35ff","#8491b4ff","#91d1c2ff"))+
  ggtitle('Before batch') +
  NoLegend()

#after-Tcell
P4 = DimPlot(T_cell, reduction = "umap",label = F, group.by = "batch",raster=FALSE) + 
  scale_color_manual(values = c("#E64b35ff","#8491b4ff","#91d1c2ff"))+
  ggtitle('After batch')



pdf(file = 'Supple_Figure/Supp_Figure1/Fig_S1A.pdf',width = 20, height = 5) 
ggarrange(P1,P2,P3,P4,nrow = 1,widths = c(1,1,1,1.1))
dev.off()




# ------------------------------------------------------------------------------
# Figure S1B
rm(list = ls())
gc()

# Function
tmp_diversity_func <- function(x){
  Shannon_diversity = diversity(table(x), index = 'shannon') # 
}

diversity_index = function(before_srt, after_srt){
  after_batch = after_srt@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(diversity=tmp_diversity_func(orig.ident)) %>% as.data.frame()
  after_batch$type='after'
  before_batch = before_srt@meta.data %>%
    group_by(seurat_clusters) %>%
    dplyr::summarise(diversity=tmp_diversity_func(orig.ident)) %>% as.data.frame()
  before_batch$type='before'
  a = rbind(after_batch, before_batch) %>%
    mutate(type = factor(type, levels = c('before', 'after')))
  return(a)
}

process_seurat <-function(obj, dim=30, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  vf_gene = VariableFeatures(obj)
  obj <- ScaleData(obj, features = vf_gene)
  obj <- RunPCA(obj, features = vf_gene, npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "pca", n.components=n.components)
  return(obj)
}

process_seurat_harmony <-function(obj, dim=30, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  vf_gene = VariableFeatures(obj)
  obj <- ScaleData(obj, features = vf_gene)
  obj <- RunPCA(obj, features = vf_gene, npcs = 50, verbose = TRUE)
  obj = RunHarmony(obj, "batch", plot_convergence = TRUE)
  obj <- FindNeighbors(obj, reduction = "harmony",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "harmony", n.components=n.components)
  return(obj)
}

# # Tcell diversity

load('Final_RDS/supplement_Figure_S1A_Tcell.Rdata')
Fig.S1B_Tcell <- diversity_index(before_srt = T_srt_before, after_srt = T_srt)
Fig.S1B_Tcell$CellType = 'T cell';head(Fig.S1B_Tcell)


load('Final_RDS/supplement_Figure_S1B_Bcell.Rdata')
Fig.S1B_Bcell = diversity_index(srt_B_before, srt_B_after)
Fig.S1B_Bcell$CellType = 'B cell';head(Fig.S1B_Bcell)

load('Final_RDS/supplement_Figure_S1B_Myeloidcell.Rdata')
Fig.S1B_Myeloid = diversity_index(srt_Myeloid_before, srt_Myeloid_after)
Fig.S1B_Myeloid$CellType = 'Myeloid';head(Fig.S1B_Myeloid)

load('Final_RDS/supplement_Figure_S1B_Fibroblast.Rdata')
Fig.S1B_Fibroblast = diversity_index(srt_Fibroblast_before, srt_Fibroblast_after)
Fig.S1B_Fibroblast$CellType = 'Fibroblast';head(Fig.S1B_Fibroblast)

load('Final_RDS/supplement_Figure_S1B_Endothelial.Rdata')
Fig.S1B_Endothelial = diversity_index(srt_Endothelial_before, srt_Endothelial_after)
Fig.S1B_Endothelial$CellType = 'Endothelial';head(Fig.S1B_Endothelial)


Fig.S1B <- rbind(Fig.S1B_Tcell,Fig.S1B_Bcell,Fig.S1B_Myeloid,Fig.S1B_Fibroblast,Fig.S1B_Endothelial);head(Fig.S1B)
Fig.S1B$CellType <- factor(Fig.S1B$CellType,levels = c('T cell','B cell','Myeloid','Fibroblast','Endothelial'))


pdf(file = 'Supple_Figure/Supp_Figure1/Fig_S1B.pdf',width = 6,height = 3)
ggplot(data = Fig.S1B,aes(x = CellType,y = diversity, color = type))+
  geom_boxplot(position = position_dodge(0.8),
               width= 0.6, size=0.4,outlier.shape = NA) +
  geom_jitter(aes(fill = type),position = position_jitterdodge(jitter.width = 0.2),shape = 21)+
  stat_compare_means(aes(group = type), label = "p.format",label.y = 3.2)+
  scale_color_nejm() +
  scale_fill_nejm(alpha = 0.6)+
  xlab('Cell Type') +
  ylab('Diversity(shannon)') +
  theme_classic() +
  theme(axis.text = element_markdown(color = 'black'))    
dev.off()

# Figure S1C 
# Before
load('Final_RDS/supplement_Figure_S1A_Tcell.Rdata')
Before_cell.prop <- as.data.frame(prop.table(table(T_cell_before$seurat_clusters, T_cell_before$orig.ident)))
colnames(Before_cell.prop)<-c("seurat_clusters","orig.ident","proportion")

P1 = ggplot(Before_cell.prop,aes(y= seurat_clusters,x = proportion,fill=orig.ident))+
  geom_bar(stat="identity",position="fill")+ 
  scale_fill_manual(values = cb_pattern )+
  theme_bw() +
  theme(axis.text = element_markdown(colour = 'black'),
        legend.position = 'no')

after_cell.prop <- as.data.frame(prop.table(table(T_cell$seurat_clusters, T_cell$orig.ident)))
colnames(after_cell.prop)<-c("seurat_clusters","orig.ident","proportion")
after_cell.prop <- subset(after_cell.prop,seurat_clusters!='12'&seurat_clusters !='13')

P2 = ggplot(after_cell.prop,aes(y = seurat_clusters,x = proportion,fill=orig.ident))+
  geom_bar(stat="identity",position="fill",na.rm = T)+ 
  scale_fill_manual(values = cb_pattern )+
  theme_bw() +
  theme(axis.text = element_markdown(colour = 'black'),
        legend.position = 'right')

pdf(file = 'Supple_Figure/Supp_Figure1/Fig_S1C.pdf',width = 8,height = 3)
P1 + P2
dev.off()



# ------------------------------------------------------------------------------
rm(list = ls())
gc()

# Figure S1D 
# cell marker experssion
load('Final_RDS/supplement_Figure_S1A_All_cell.Rdata')
pdf('Supple_Figure/Supp_Figure1/Fig_S1D.pdf',width = 18,height = 12)
FeaturePlot(All_cell_post_harmony,features = c("KRT5","CD3D","CD3E","MS4A1","MZB1","CD68","LYZ","TPSAB1","VWF","DCN","CXCL8","S100A8"),raster = FALSE)
dev.off()


# ------------------------------------------------------------------------------
# Figure S1E
# cell type
pdf(paste0(output_path,'/Figure S1E.pdf'), width =10, height=10)
DimPlot(All_cell_post_harmony, reduction = "umap",label = T,pt.size = 0.01,label.size = 4, group.by = "cell.type", raster=FALSE) + 
  scale_color_manual(values = color_cell_type)+
  NoLegend()
dev.off()
# Figure S1F
# Response/Treatment
P1 = DimPlot(All_cell_post_harmony, reduction = "umap", label = F, group.by = "response", raster=FALSE) + 
  scale_color_manual( values = c("#0099CCFF","#FFC20AFF"))  

P2 = DimPlot(All_cell_post_harmony, reduction = "umap",label = F, group.by = "time", raster=FALSE) +
  scale_color_manual(values = c('#91d1c2ff','#f768a1'))  

pdf(file = 'Supple_Figure/Supp_Figure1/Fig_S1F.pdf',width = 12,height = 5)
P1 + P2
dev.off()

