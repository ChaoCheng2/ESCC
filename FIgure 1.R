#' ------------------------------------------------------------------------------------
#' Figure1 

rm(list = ls())

# library
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(dplyr)

library(magrittr)
library(cowplot)
library(egg)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
 
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

zylcolor40 <- ggsci::pal_igv()(40)
color_palette_Response <- c('Non_responder' = "#0099CCFF",'Responder' = "#FFC20AFF")
color_cell_type = c("B cell" = "#E64B35FF", "Endothelial cell"= "#4DBBD5FF", "Fibroblast" = "#00A087FF", "Macrophage/Monocyte" = "#3C5488FF", 
                     "Mast cell" = "#F39B7FFF", "Neutrophils"= "#8491B4FF", "Plasma" = "#91D1C2FF",  "T cell" = "#fdbf6f", "Tumor or epithelial cells" = "#7E6148FF")
color_palette_Treatment = c('Post_T'='#91d1c2ff', 'Pre_T'= '#f768a1')

# --------------------------------------------------------------------------------
# Setwd path
workdir <- "/Volumes/Camellia/Project_ESCC" ; setwd(workdir)
output_path = 'Figure_output/Figure1_original'
dir.create(file.path(workdir,output_path))

# --------------------------------------------------------------------------------

All_cell_post_harmony <- readRDS("./Final_RDS/All_cells_after_harmony_10.13.rds")
Tcell_srt <- readRDS("./Final_RDS/T_cell.rds")
T_cell_no_NK <- readRDS("./Final_RDS/T_cell_no_NK_updated_8.8.rds")

# -------------------------------------------------------------------------------
gene_test = c("IGFBP5","ACKR1","FABP5","SLCO2A1","KCNIP4","CD74","HLA-DRA","PLAT","TFF3","HLA-DRB1","FABP4","SELE","MT2A","LYZ","FLRT2","ANO2","ADIRF","ZNF385D","NR2F2","IFI30","S100A10","S100A6","GPM6A","NRG3","VCAM1","HLA-DQB1","NNMT","CEBPD")
png('gene.test.png',width = 1200,height = 1500, res = 100)
FeaturePlot(All_cell_post_harmony, features = gene_test,raster=FALSE)
dev.off()

All_cell_post_harmony$gene_test_mean = colMeans(as.matrix(All_cell_post_harmony@assays[['RNA']]@data[gene_test,]))

png('gene.test2.png',width = 500,height = 400, res = 100)
FeaturePlot(All_cell_post_harmony, features = 'gene_test_mean',raster=FALSE)
dev.off()

##----------------------------------------------------------

dir.create('Test')
gene_test = read.csv('/Volumes/Camellia/Project_ESCC/DATA/endo_c5_sigs_1.csv')
head(gene_test)
gene =gene_test$c5.sigs
length(unique(gene))
png('./Test/gene.test_11.19_1.png',width = 4000,height = 3000, res = 200)
FeaturePlot(All_cell_post_harmony, features = gene[1:30],raster=FALSE,ncol = 6)
dev.off()
png('./Test/gene.test_11.19_2.png',width = 4000,height = 3000, res = 200)
FeaturePlot(All_cell_post_harmony, features = gene[31:60],raster=FALSE,ncol = 6)
dev.off()
png('./Test/gene.test_11.19_3.png',width = 4000,height = 3000, res = 200)
FeaturePlot(All_cell_post_harmony, features = gene[61:90],raster=FALSE,ncol = 6)
dev.off()

All_cell_post_harmony$gene_test_mean = colMeans(as.matrix(All_cell_post_harmony@assays[['RNA']]@data[gene,]))

png('./Test/Mean_2.png',width = 500,height = 400, res = 100)
FeaturePlot(All_cell_post_harmony, features = 'gene_test_mean',raster=FALSE)
dev.off()

# ------------------------------------
endo = subset(All_cell_post_harmony,cell.type == 'Endothelial cell')

pdf('./Test/cellType.pdf', width =10, height=10)
DimPlot(endo, reduction = "umap",label = T,pt.size = 0.01,label.size = 4, group.by = "cell.type", raster=FALSE) + 
    scale_color_manual(values = color_cell_type)+
    NoLegend()
dev.off()
 

endo = process_seurat(endo)

pdf('./Test/cellType2.pdf', width =5, height=5)
DimPlot(endo, reduction = "umap",label = T,pt.size = 0.01,label.size = 4, group.by = "cell.type", raster=FALSE) + 
    scale_color_manual(values = color_cell_type)+
    NoLegend()
dev.off()
aa = pal_igv()(15)
aa

pdf('./Test/cellType4.pdf', width =5, height=4)
DimPlot(endo, reduction = "umap",label = T,pt.size = 0.01,label.size = 4, group.by = c( "response"), raster=FALSE) +
 scale_color_manual(values = color_palette_Response)
dev.off()

png('./Test/Mean_endo.png',width = 500,height = 400, res = 100)
FeaturePlot(endo, features = 'gene_test_mean',raster=FALSE)
dev.off()


pdf('./Test/vln.pdf', width =5, height=5)
VlnPlot(endo,features = 'gene_test_mean',group.by = 'response')
dev.off()

png('./Test/gene.test_11.19_endo1.png',width = 4000,height = 3000, res = 200)
FeaturePlot(endo, features = gene[1:30],raster=FALSE,ncol = 6)
dev.off()
png('./Test/gene.test_11.19_endo2.png',width = 4000,height = 3000, res = 200)
FeaturePlot(endo, features = gene[31:60],raster=FALSE,ncol = 6)
dev.off()
png('./Test/gene.test_11.19_endo3.png',width = 4000,height = 3000, res = 200)
FeaturePlot(endo, features = gene[61:90],raster=FALSE,ncol = 6)
dev.off()

#-----------------------------------------

library(ggplot2)
 
# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
data = endo@meta.data[,c('seurat_clusters','response')]
data = table(data$seurat_clusters,data$response)
head(data)
data = as.data.frame(data)
table(data$Var2)

library(tidyr)

data2= gather(data)
head(data2)
melt()
head(data)
library(viridis)
library(hrbrthemes)

# Stacked
pdf('./Test/ggplot_barplot.pdf',width =4,height = 3)
ggplot(data, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = color_palette_Response) +
   
    geom_text(aes(label = Freq),position=position_stack(vjust =0.5),size =1.8)+ 
    #geom_text(aes(label = Freq,y = Freq-10))+
    theme_few()+    
    xlab("")
dev.off()

pdf('./Test/ggplot_VlnPlot.pdf',width = 6,height = 3)
VlnPlot(endo,features = 'gene_test_mean',split.by ='response',pt.size = 0.001 )+
    scale_fill_manual(values = color_palette_Response) 
dev.off()

# -----------------------------------------
# Figure 1A
#study design


# --------------------------------------------------------------------------------
# Figure 1B
pdf(paste0(output_path,'/Figure1B.pdf'), width =10, height=10)
DimPlot(All_cell_post_harmony, reduction = "umap",label = T,pt.size = 0.01,label.size = 4, group.by = "cell.type", raster=FALSE) + 
    scale_color_manual(values = color_cell_type)+
    NoLegend()
dev.off()

# --------------------------------------------------------------------------------
# Figure 1C
# FeaturePlot(All_cell_post_harmony, features = 'PDCD1',pt.size = 0.2,raster=FALSE) 
Figure1C <- data.frame(All_cell_post_harmony@reductions$umap@cell.embeddings)
Figure1C$PDCD1 <-  All_cell_post_harmony@assays$RNA@data[c('PDCD1'),]
Figure1C_gray <- subset(Figure1C,PDCD1==0)
Figure1C_red <- subset(Figure1C,PDCD1!=0)

pdf(paste0(output_path,'/Figure1C.pdf'), width =11, height=10)
ggplot()+
    geom_point(data = Figure1C_gray,mapping = aes(x = UMAP_1, y = UMAP_2,color = PDCD1),size = 0.01,inherit.aes = F) +
    geom_point(data = Figure1C_red,mapping = aes(x = UMAP_1, y = UMAP_2,color = PDCD1),size = 0.01,inherit.aes = F) +
    scale_color_gradientn(colours = c("lightgrey", "blue")) +
    theme_classic()
dev.off()

# --------------------------------------------------------------------------------
# Figure 1D
T_pre <- subset(Tcell_srt,time=="Pre_T")
IL7R_exp = FetchData(T_pre, vars="PDCD1")
IL7R_exp$sample = T_pre@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = T_pre@meta.data[rownames(IL7R_exp), "response"]
IL7R_exp$time = T_pre@meta.data[rownames(IL7R_exp), "time"]

IL7R_exp_mean = IL7R_exp %>%
    group_by(sample,response,time) %>%
    dplyr::summarise_all(list(mean)) %>% as.data.frame()

P1 = ggplot(data = IL7R_exp_mean, aes(x=response, y=PDCD1))+
    geom_boxplot(color = 'black',width=0.7, size=0.4,outlier.shape = NA)+
    geom_jitter(aes(color=response),width = 0.2,size=2) +
    stat_compare_means(comparisons = list(c('Responder', 'Non_responder')), label = 'p.format',face="bold",method = "t.test")+
    scale_color_manual(values = color_palette_Response) + 
    xlab("") +
    ylab("PD1 expression level in T cells pre-treatment") +
    theme_classic( )+
    theme(
        legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black"), 
        axis.line.y = element_line(),
        axis.line = element_line(colour = "black"))

seurat_meta = T_pre@meta.data
seurat_meta$Exp = FetchData(T_pre, vars = "PDCD1")
# counts of exp > 0
countsExp <- function(x){
    sum(x>0)
}

count_exp_data = seurat_meta %>% 
    group_by(orig.ident,time,response) %>%
    dplyr::summarise(ce=countsExp(Exp)) %>% as.data.frame()

sample_cell_counts = seurat_meta %>% 
    group_by(orig.ident) %>%
    dplyr::summarise(n=n()) %>% as.data.frame()

plot_counts_PD1 = merge(count_exp_data, sample_cell_counts, by='orig.ident')
plot_counts_PD1$p = plot_counts_PD1$ce / plot_counts_PD1$n
plot_counts_PD1<-select(plot_counts_PD1,-c("ce","n"))
colnames(plot_counts_PD1)<-c("sample","time","response","PD1_cell_prop")

P2 = ggplot(data = plot_counts_PD1, aes(x=response, y=PD1_cell_prop))+
    geom_boxplot(color = 'black',width=0.7, size=0.4,outlier.shape = NA)+
    geom_jitter(aes(color=response),width = 0.2,size=2) +
    stat_compare_means(comparisons = list(c('Responder', 'Non_responder')), label = 'p.format',face="bold",method = "t.test")+
    scale_color_manual(values = color_palette_Response) + 
     
    xlab("") +
    ylab("% PD1+ T cells pre-treatment") +
    scale_y_continuous(position = c('left'))+
    theme_classic(base_line_size = 0.3)+
    theme(
        legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black"), 
        axis.line.y = element_line(),
        axis.line = element_line(colour = "black"))

pdf(paste0(output_path,'/Figure1D.pdf'),width = 4,height = 4)
P1+P2
dev.off()


# ---------------------------------------------------------------------------------
# Figure 1E
P1 = DimPlot(Tcell_srt, reduction = "umap",label = F, pt.size = 0.4, group.by = "cluster_name", raster=FALSE,cols = zylcolor40) + NoLegend()
P2 = DimPlot(Tcell_srt, reduction = "umap",label = F, pt.size = 0.4, group.by = "response",raster=FALSE) +  scale_color_manual(values = color_palette_Response) +  NoLegend()
P3 = DimPlot(Tcell_srt, reduction = "umap",label = F, group.by = "time",pt.size = 0.4, raster=FALSE) + scale_color_manual(values =color_palette_Treatment) + NoLegend()

pdf(paste0(output_path,'/Figure1E.pdf'), width = 18, height=7)
P1+P2+P3
dev.off()

# --------------------------------------------------------------------------------
# Figure 1F
source('Code/Function/ccGene.R')
tmeta = Tcell_srt@meta.data[,c('orig.ident', 'response', 'time')]
tmeta = tmeta[!duplicated(tmeta),]
rownames(tmeta) = tmeta$orig.ident

# load cytosig
cytosig = t(read.table('./Final_RDS/Cytoseq_output/output_only_sample_expand.Zscore', sep='\t'))#
new_Tcell_srt = AddModuleScore(Tcell_srt, features = list('cc_score'=ccgene), name='cc_score')
cc_exp_score = new_Tcell_srt@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(score = median(cc_score1))
cc_exp = cc_exp_score$score
names(cc_exp) = cc_exp_score$orig.ident

plodata = as.data.frame(cytosig)
plodata$Proliferation = cc_exp[rownames(cytosig)]
plodata[, c('orig.ident','response', 'time')] = tmeta[rownames(cytosig), c('orig.ident','response', 'time')]
plodata = plodata %>% 
  filter(orig.ident!='EC_07_post_T') #%>%


cor_res = c()
for(cyto in colnames(cytosig)){
  cor_res = rbind(cor_res, c(cyto,cor(plodata[, cyto], plodata$Proliferation)))
}
cor_res = as.data.frame(cor_res)
colnames(cor_res) = c('cytosig', 'cor')
cor_res$cor = as.numeric(cor_res$cor)
cor_res = cor_res[order(cor_res$cor),]
 
# add zzm score
new_Tcell_srt = AddModuleScore(Tcell_srt, features = list('cc_score'=ccgene_zzm), name='cc_score')
cc_exp_score = new_Tcell_srt@meta.data %>%
  group_by(orig.ident) %>%
  summarise(score = median(cc_score1))
cc_exp = cc_exp_score$score
names(cc_exp) = cc_exp_score$orig.ident

plodata$cor_ZZM = cc_exp[rownames(plodata)]

cor_ttest_res = c()
for(cyto in colnames(cytosig)){
  #cyto = 'TRAIL'
  tmp_cor =  cor(plodata[, cyto], plodata$Proliferation)
  tmp_cor_zzm = cor(plodata[, cyto], plodata$cor_ZZM)
  
  tmp_ttest = wilcox.test(plodata[plodata$response=='Responder', cyto], plodata[plodata$response=='Non_responder', cyto])
  tmp_fc = median(plodata[plodata$response=='Non_responder', cyto]) - median(plodata[plodata$response=='Responder', cyto])
  cor_ttest_res = rbind(cor_ttest_res, c(cyto,tmp_cor,tmp_ttest$p.value, tmp_fc, tmp_cor_zzm))
}

cor_ttest_res = as.data.frame(cor_ttest_res)
colnames(cor_ttest_res) = c('cytosig', 'cor_Tres', 'Diff_pvalue', 'NR_sub_R', 'cor_ZZM')
cor_ttest_res$cor_Tres = as.numeric(cor_ttest_res$cor_Tres)
cor_ttest_res$Diff_pvalue = as.numeric(cor_ttest_res$Diff_pvalue)
cor_ttest_res$NR_sub_R = as.numeric(cor_ttest_res$NR_sub_R)
cor_ttest_res$cor_ZZM = as.numeric(cor_ttest_res$cor_ZZM)
cor_ttest_res$size = 1#cor_ttest_res$Diff_pvalue<0.05
cor_ttest_res$size[cor_ttest_res$Diff_pvalue<0.05] = 1.5

library(forcats)
library(ggrepel)
library(ggtext)
pdata = cor_ttest_res %>%
  dplyr::arrange(Diff_pvalue) %>%
  dplyr::mutate(cytosig=fct_reorder(cytosig, -NR_sub_R)) %>%
  dplyr::mutate(x=1:51)

pdata$label <- ifelse(pdata$Diff_pvalue<=0.05, as.character(pdata$cytosig), NA)

# ------------------------------------------------------------------------------

pdf(paste0(output_path,'/Figure1F.pdf'),width = 6,height = 4.5)
ggplot(data = pdata, aes(x=cor_Tres, y = -log10(Diff_pvalue))) +
    geom_point(aes(size = -log10(Diff_pvalue), fill = cor_Tres), shape=21) +
    geom_text_repel(aes(label = label),color = 'black',  size= 3, max.overlaps = 1000,seed = 7654) +
    
    #scale_fill_gradient2(midpoint = 0, low = '#2171b5',mid = 'white',high =  '#cc4c02',name='Correlation') +
    scale_fill_gradientn(limits = c(-0.8,0.8),colours = c("#053061",'#2171b5', 'white','#cc4c02',"#67000d"),breaks =  c(-0.5,0,0.5),name='Correlation')+
    scale_size(range = c(0.5, 6)) +
    ylim(0,1.5)+
    geom_hline(yintercept = -log10(0.05),linetype = 'dashed')+
    xlab('Correlation') +
    ylab('-log10(P-Value)')+
    theme_few() +
    theme(
        panel.border = element_rect(color = "black", fill=NA),
        axis.text = element_markdown(color = 'black'),
        legend.text = element_text(size = 8),
        panel.grid = element_blank())
dev.off()

# ------------------------------------------------------------------------------
# Supplement 
pdata$Name = 'name'
data_NR_polygon = data.frame(x=c(0.5, 0.5, nrow(cor_ttest_res), 0.5),y=c(1.2,1.5, 1.2, 1.2))

p1 <- ggplot() +
    geom_point(data=pdata, aes(x = cytosig, y = Name, size = -log10(Diff_pvalue),fill = cor_Tres), shape=21,inherit.aes = F) +
    #geom_text_repel(aes(label = cytosig),color = 'black',  size= 3, max.overlaps = 1000,seed = 7654) +
    scale_fill_gradientn(limits = c(-0.8,0.8),colours = c("#053061",'#2171b5', 'white','#cc4c02',"#67000d"),breaks =  c(-0.5,0,0.5),name='Correlation')+
    scale_size(range = c(0.5, 8)) +
    theme_few() +
    xlab('Correlation') +
    ylab('-log10(P-Value)')+
    theme(
        legend.position = 'bottom',
        panel.border = element_rect(color = "black", fill=NA),
        axis.text.x = element_markdown(color = 'black',angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(colour = NA),
        legend.text = element_text(size = 8),
        panel.grid = element_blank())
p2 <- ggplot()+
    geom_polygon(aes(x=x, y=y), data=data_NR_polygon, fill='#deebf7',alpha=0.5, color='gray') +
    scale_x_discrete(expand = expansion()) +
    #scale_y_discrete(expand = expansion(add = 0.5)) +
    theme_void()

library(aplot)
pdf('Supple_Figure/Supp_Figure2/FigS2G_hh.pdf',width = 18,height = 2)
p1%>%insert_top(p2,height = 0.5)
dev.off()


# ----------------------------------------------------------------------------- 
####Figure1G #######
cyto_res = t(read.table('./Final_RDS/Cytoseq_output/new_updated/output_cell_expand_noNK.Zscore', sep='\t'))
cyto_res_sample = t(read.table('./Final_RDS/Cytoseq_output/new_updated/output_sample_expand_noNK.Zscore', sep='\t'))
cyto_res_TRAIL = cyto_res_sample[, c('TRAIL')]
cyto_res = cyto_res[, c('TGFB1', 'PGE2', 'TRAIL')]
rownames(cyto_res) = gsub('\\.', '-',rownames(cyto_res))

T_cell_no_NK <- readRDS("./Final_RDS/T_cell_no_NK_updated_8.8.rds")
Tres_score = T_cell_no_NK$Tres_proliferation_score1
#
cyto_res = as.data.frame(cyto_res)
cyto_res$proliferation = Tres_score[rownames(cyto_res)]
plodata = cyto_res
plodata$response = T_cell_no_NK@meta.data[rownames(plodata), 'response']
plodata$sample = T_cell_no_NK@meta.data[rownames(plodata), 'orig.ident']

PD1_exp = FetchData(T_cell_no_NK, vars="PDCD1")
plodata$PDCD1 = PD1_exp[rownames(plodata), 'PDCD1']

##
save = plodata %>%
  filter(sample !="EC_07_post_T") %>%filter(sample !="EC_16_pre_T") %>%
  group_by(sample) %>%
  dplyr::summarise(TGFB1=cor(TGFB1, proliferation), 
            PDCD1=cor(PDCD1, proliferation), 
            PGE2=cor(PGE2, proliferation),
            TRAIL=cor(TRAIL, proliferation))
save[is.na(save)] = 0

### zhangzemin
ccgene = c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
           'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
           'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1')

T_cell_no_NK = AddModuleScore(T_cell_no_NK, features = list('cc_score'=ccgene), name='Proliferation_score_ZZM')
cc_exp = T_cell_no_NK$Proliferation_score_ZZM1
#
cyto_res = as.data.frame(cyto_res)
cyto_res$proliferation = cc_exp[rownames(cyto_res)]
plodata = cyto_res
plodata$response = T_cell_no_NK@meta.data[rownames(plodata), 'response']
plodata$sample = T_cell_no_NK@meta.data[rownames(plodata), 'orig.ident']

PD1_exp = FetchData(T_cell_no_NK, vars="PDCD1")
plodata$PDCD1 = PD1_exp[rownames(plodata), 'PDCD1']

##
save2 = plodata %>%
  filter(sample !="EC_07_post_T") %>%filter(sample !="EC_16_pre_T") %>%
  group_by(sample) %>%
  dplyr::summarise(TGFB1=cor(TGFB1, proliferation), 
            PDCD1=cor(PDCD1, proliferation), 
            PGE2=cor(PGE2, proliferation),
            TRAIL=cor(TRAIL, proliferation))
save2[is.na(save2)] = 0

#
library(tidyr)
saveData_long <- gather(save2[,c('sample','PDCD1','TRAIL')],key = 'gene',value = 'sample')

pdf(paste0(output_path,'/Figure1G.pdf'),width = 5,height = 4)
ggplot(saveData_long, aes(sample, fill = gene)) + 
    geom_density(adjust = 1.5,alpha = 0.6,size = 0.8) + 
    xlim(-0.7,0.3)+
    scale_fill_manual(values = c('#8A94DB','#F4E0B9'))+
    theme_few()+
    xlab('Correlation with proliferation')+
    ylab('Density')+
    theme(
        panel.border = element_rect(color = "black", fill=NA),
        axis.text.x = element_markdown(colour = 'black'),
        axis.text.y = element_markdown(colour = 'black'),
        legend.justification = c(1,0),
        legend.text = element_text(size = 8),panel.grid = element_blank())
dev.off()


# ------------------------------------------------------------------------------
# Figure1H####
# cytoseq
output_dir = './Final_RDS/Cytoseq_output/new_updated/' # output path

# draw

cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample <- cyto_res[,"TRAIL"]
TRAIL_sample <- as.data.frame(TRAIL_sample)


TRAIL_exp = FetchData(T_cell_no_NK, vars="TNFSF10")
TRAIL_exp$sample = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "orig.ident"]
TRAIL_exp$time = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "time"]
TRAIL_exp$response = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "response"]
TRAIL_exp_mean = TRAIL_exp %>%
  group_by(sample,time,response) %>%
  summarise_all(list(mean)) %>% as.data.frame()


TRAIL_sample_all <- cbind(TRAIL_exp_mean,TRAIL_sample)
TRAIL_sample_all <- TRAIL_sample_all[,c(1,2,3,5)]
colnames(TRAIL_sample_all) <- c("orig.ident","time","response","TRAIL_Signaling")
TRAIL_sample_all_no_16pre <- subset(TRAIL_sample_all,orig.ident!="EC_16_pre_T")
TRAIL_sample_all_no_16pre <- subset(TRAIL_sample_all_no_16pre,orig.ident!="EC_07_post_T")
TRAIL_sample_pre <- subset(TRAIL_sample_all_no_16pre,time=="Pre_T")

Figure1F <- TRAIL_sample_pre
Figure1F$Method = 'Single Cell'

###Figure1G-BULK####
library(clusterProfiler)
library(fgsea)
library(forcats)
library(enrichplot)
library(Seurat)
dataMatrix = read.table('./Final_RDS/EC_bulk_count.txt', sep='\t')
sample_info = read.table('./Final_RDS/sample_info.txt', header = F, sep='\t')
colnames(sample_info) = c('sample', 'time', 'response')
sample_info$response2 = sapply(sample_info$response, function(x)ifelse(x%in%c("CR","MPR","CR "), 'Response', 'Non_response'))
sample_info$response = sample_info$response2
sample_info$sample = gsub('-', '_', sample_info$sample)
rownames(sample_info) = sample_info$sample
sample_info = sample_info[sample_info$time%in%c('Pre-T'),]
share_name = intersect(sample_info$sample, colnames(dataMatrix))
dataMatrix = dataMatrix[,share_name]
sample_info = sample_info[share_name, ]
dataMatrix = dataMatrix[rowSums(dataMatrix>10)  > (ncol(dataMatrix)/4),]

###
library(limma)
library(edgeR)
dataMatrix <- DGEList(counts=dataMatrix)
group <- factor(sample_info$response, levels = c( 'Non_response','Response'))
dataMatrix$samples$group = sample_info[rownames(dataMatrix$samples), 'response']
keep <- filterByExpr(dataMatrix, group=group)
dataMatrix <- dataMatrix[keep, , keep.lib.sizes=FALSE]

gene_length = read.table('./Final_RDS/hg19_gene.txt')
gene_vector = gene_length$width
names(gene_vector) = gene_length$gene_id
logCPM <- rpkm(dataMatrix,gene_vector, log = T,prior.count=1)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) = colnames(dataMatrix)
# 
fit <- lmFit(logCPM, design)

df.matrix <- makeContrasts('Non_response-Response', levels = design)
fit <- contrasts.fit(fit, df.matrix)
fit2 <- eBayes(fit)

allDiff=topTable(fit2, coef=1, n=Inf)
allDiff = allDiff[order(allDiff$logFC, decreasing = T), ]
genelist = allDiff$logFC
names(genelist) = rownames(allDiff)

####
# Cytosig TRAIL gene
TRAILmark <- read.gmt('./Final_RDS/TRAIL_pathway.gmt')
TRAILmark = TRAILmark %>% split(.$term) %>% lapply('[[',2)
TRAILmark_gene = c(#TRAILmark$GOBP_TRAIL_ACTIVATED_APOPTOTIC_SIGNALING_PATHWAY,
  TRAILmark$PID_TRAIL_PATHWAY,
  TRAILmark$GOMF_TRAIL_BINDING,
  TRAILmark$REACTOME_TRAIL_SIGNALING
)
TRAIL_marker = read.table('./Final_RDS/signature.centroid.expand', sep='\t')
TRAIL_marker = TRAIL_marker[,'TRAIL', drop=F]
TRAIL_marker = TRAIL_marker[order(TRAIL_marker$TRAIL, decreasing = T), , drop=F]
# hallmark = list('TRAIL' = c(rownames(head(TRAIL_marker, 100))))
#target_gene = c('CXCL2', 'CCL20', 'CXCL8', 'TNFAIP3', 'BIRC3', 'IRF1', 'CXCL1', 'NFKBIA', 'FOS', 'JUNB', 'PPP1R15A', 'EGR1', 'PLK3', 'RHOB')
#term2gene = data.frame(name='TRAIL', gene=c(target_gene),TRAILmark_gene)#,  'PTPRC', 'CD3D', 'CD3E', 'CD4', 'CD8A'))
term2gene = data.frame(name='TRAIL', gene=c(rownames(head(TRAIL_marker, 20)),TRAILmark_gene,'CD3D', 'CD3E'))
# term2gene = data.frame(name='TRAIL', gene=c(rownames(TRAIL_marker[TRAIL_marker>0.5,,drop=F]),TRAILmark_gene,'CD3D', 'CD3E'))
term2gene = subset(term2gene, gene%in%names(genelist))
term2gene = subset(term2gene, gene!='IKBKG')
term2gene = subset(term2gene, gene!='IKBKB')
term2gene = subset(term2gene, gene!='NFKB2')
term2gene = subset(term2gene, gene!='CXCL8')
term2gene = subset(term2gene, gene!='TRIB1')

pre_bulk_data = logCPM[term2gene$gene, ]

library(reshape2)
library(ggpubr)
#pre_bulk_data = as.data.frame(t(apply(pre_bulk_data, 1, scale)))
#pheatmap::pheatmap(pre_bulk_data)
#colnames(pre_bulk_data) = colnames(logCPM)
pre_bulk_data = melt(as.matrix(pre_bulk_data))

pre_bulk_data$response = sample_info[pre_bulk_data$Var2, 'response']
pre_bulk_data$response[pre_bulk_data$response=='Response'] = 'Responder'
pre_bulk_data$response[pre_bulk_data$response=='Non_response'] = 'Non_responder'
pre_bulk_data$response = factor(pre_bulk_data$response, levels=c('Responder', 'Non_responder'))

Figure1G <- pre_bulk_data %>%
    #filter(Var1%in%names(sort(genelist[term2gene$gene], decreasing = T)[1:50])) %>%
    group_by(Var2, response) %>%
    summarise(TRAIL_Signaling=mean(value))

Figure1G$Method <- 'Bulk'
Figure1F_1G <- rbind(Figure1F[,c('response','TRAIL_Signaling','Method')],Figure1G[,c( "response","TRAIL_Signaling","Method" )])
Figure1F_1G$Method <- factor(Figure1F_1G$Method,levels = c('Single Cell','Bulk'))
Figure1F_1G$response <- factor(Figure1F_1G$response,levels = c('Non_responder','Responder'))
save(Figure1F_1G,file = 'DATA/Figure1G.RData')


pdf(paste0(output_path,'/Figure1H.pdf'),width = 3.5,height = 4)
ggplot(data = Figure1F_1G,aes(x = response,y = TRAIL_Signaling))+
    geom_boxplot(width = 0.7,color = '#4d4d4d', outlier.size = 0.0,outlier.color = 'white',size = 0.6)+
    geom_jitter(aes(color = response),width = 0.15,size = 2.5)+
    scale_color_manual(values = color_palette_Response)+
    stat_compare_means(comparisons = list(c('Responder', 'Non_responder')),label = "p.format",face="bold")+
    facet_wrap(~Method,scales = 'free')+
    theme_bw()+
    xlab('Response')+
    ylab('TRAIL Signaling')+
    NoLegend()+
    theme(
        strip.background = element_rect(color="black", fill="white",  linetype="solid" ),
        panel.border = element_rect(color = "black" ,),
        axis.text.x = element_markdown(colour = 'black'),
        axis.text.y = element_markdown(colour = 'black'),
        legend.text = element_text(),
        panel.grid = element_blank())
dev.off()

 

# ------------------------------------------------------------------------------------
# Figure 1I
mygseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                         rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                         ES_geom = "line",ymax = 0.75) {
    
    gseaScores<-function (geneList, geneSet, exponent = 1, fortify = FALSE) {
        geneSet <- intersect(geneSet, names(geneList))
        N <- length(geneList)
        Nh <- length(geneSet)
        Phit <- Pmiss <- numeric(N)
        hits <- names(geneList) %in% geneSet
        Phit[hits] <- abs(geneList[hits])^exponent
        NR <- sum(Phit)
        Phit <- cumsum(Phit/NR)
        Pmiss[!hits] <- 1/(N - Nh)
        Pmiss <- cumsum(Pmiss)
        runningES <- Phit - Pmiss
        max.ES <- max(runningES)
        min.ES <- min(runningES)
        if (abs(max.ES) > abs(min.ES)) {
            ES <- max.ES
        }
        else {
            ES <- min.ES
        }
        df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
                         position = as.integer(hits))
        if (fortify == TRUE) {
            return(df)
        }
        df$gene = names(geneList)
        res <- list(ES = ES, runningES = df)
        return(res)
    }
    
    gsInfo = function (object, geneSetID) 
    {
        geneList <- object@geneList
        if (is.numeric(geneSetID)) 
            geneSetID <- object@result[geneSetID, "ID"]
        geneSet <- object@geneSets[[geneSetID]]
        exponent <- object@params[["exponent"]]
        df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
        df$ymin <- 0
        df$ymax <- 0
        pos <- df$position == 1
        h <- diff(range(df$runningScore))/20
        df$ymin[pos] <- -h
        df$ymax[pos] <- h
        df$geneList <- geneList
        df$Description <- object@result[geneSetID, "Description"]
        return(df)
    }
    
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    }
    else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) + 
        xlab(NULL) + 
        theme_classic(base_size) + 
        geom_hline(yintercept = 0, linetype='dashed')+
        theme(panel.grid.major = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(colour = "grey92"), 
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
        scale_x_continuous(expand = c(0, 0))
    
    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                              size = 1)
    }
    else {
        es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                               size = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p +
        es_layer + 
        theme(
            legend.position = c(0.8, 0.8), 
            legend.title = element_blank(), 
            legend.background = element_rect(fill = "transparent"))
    
    p.res <- p.res + 
        ylab("Running Enrichment Score") + 
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.line.x = element_blank(),  
            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"))
    
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description ==  term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    
    p2 <- ggplot(gsdata, aes_(x = ~x)) + 
        geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
        ylab(NULL) + 
        theme_classic(base_size) + 
        theme(legend.position = "none", 
              plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
              axis.ticks = element_blank(), axis.text = element_blank(), 
              axis.line.x = element_blank()) + scale_x_continuous(expand = c(0,  0)) + 
        scale_y_continuous(expand = c(0, 0))
    
    if (length(geneSetID) == 1) {
        # v <- seq(1, sum(gsdata$position), length.out = 9)
        # inv <- findInterval(rev(cumsum(gsdata$position)), v)
        # if (min(inv) == 0) 
        #   inv <- inv + 1
        # col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
        # ymin <- min(p2$data$ymin)
        # yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
        # xmin <- which(!duplicated(inv))
        # xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        # d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
        #                 xmax = xmax, col = col[unique(inv)])
        #p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
        #                          ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
        #                     alpha = 0.9, inherit.aes = FALSE)
    }
    
    
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + 
        geom_segment(data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + 
        ylab("Ranked List Metric") + 
        xlab("Rank in Ordered Dataset") + 
        theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
    
    if (!is.null(title) && !is.na(title) && title != "") 
        p.res <- p.res + ggtitle(title)
    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }
        else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
    
    if (pvalue_table) {
        pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
        rownames(pd) <- pd$Description
        pd <- pd[, -1]
        pd <- round(pd, 4)
        tp <- tableGrob2(pd, p.res)
        p.res <- p.res + 
            theme(legend.position = "none") + 
            annotation_custom(tp,
                              xmin = quantile(p.res$data$x, 0.5), 
                              xmax = quantile(p.res$data$x,  0.95), 
                              ymin = quantile(p.res$data$runningScore,0.75), 
                              ymax = quantile(p.res$data$runningScore,  0.9))
    }
    
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                           axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1) 
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                          r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    if (length(rel_heights) > length(subplots)) 
        rel_heights <- rel_heights[subplots]
    plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}


scRNA_seq <- readRDS('DATA/Figure1H_GSEA_scRNA_seq.Rds')
bulk <- readRDS(file = 'DATA/Figure1H_GSEA_bulk.Rds')

library(cowplot)
p1 <- mygseaplot2(scRNA_seq, geneSetID = 1, color="#E5614CFF", ES_geom='line', pvalue_table=F, title=NULL, rel_heights = c(1.5, 0.2), base_size=8, subplots=c(1,2)) +
    geom_text(size = 3,aes(x=0.8, y=0.84, label=paste0('p-value: ', signif(scRNA_seq@result['TRAIL', 'pvalue'], 2))), color='black' ) +
    geom_text(size = 3,aes(x=0.8, y=0.9, label=paste0('ES: ', round(scRNA_seq@result['TRAIL', 'enrichmentScore'], 2))), color='black') +
    geom_text(size = 3,aes(x=0.9, y=-0.015, label='R'), color='black')+
    geom_text(size = 3,aes(x=0.2, y=-0.015, label='NR'), color='black')+
    ggtitle('scRNA-seq')+
    theme(plot.title = element_text(vjust=0.5,size = 8),
          plot.margin = margin(t = 0.1, r = 0.2, b = 0.5, l = 0.2, unit = "cm"))

p2 <- mygseaplot2(bulk, geneSetID = 1, color="#E5614CFF",
                  ES_geom='line',
                  pvalue_table=F, title=NULL,
                  rel_heights = c(1.5, 0.2),
                  base_size=8,subplots=c(1,2))+
    geom_text(size = 3,aes(x=0.8, y=0.84, label=paste0('p-value: ', signif(bulk@result['TRAIL', 'pvalue'], 2))), color='black')+
    geom_text(size = 3,aes(x=0.8, y=0.9, label=paste0('ES: ', round(bulk@result['TRAIL', 'enrichmentScore'], 2))), color='black')+
    geom_text(size = 3,aes(x=0.9, y=-0.015, label='R'), color='black')+
    geom_text(size = 3,aes(x=0.2, y=-0.015, label='NR'), color='black')+
    ggtitle('Bulk')+
    theme(plot.title = element_text(vjust=0.5,size = 8),
          plot.margin = margin(t = 0.1, r = 0.2, b = 0.5, l = 0.2, unit = "cm"))

pdf(file = paste0(output_path,'/Figure1I.pdf'), width = 7,height =4)
p1 + p2
dev.off()

