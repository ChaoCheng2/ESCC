#' Figure 1

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
workdir <- "/Volumes/Camellia/Project_ESCC" 
output_path = 'Figure_output/Figure1_original'
dir.create(file.path(workdir,output_path))

# --------------------------------------------------------------------------------

All_cell_post_harmony <- readRDS("./Final_RDS/All_cells_after_harmony.rds")
Tcell_srt <- readRDS("./Final_RDS/T_cell.rds")
T_cell_no_NK <- readRDS("./Final_RDS/T_cell_no_NK.rds")

# -------------------------------------------------------------------------------

gene_test = c("IGFBP5","ACKR1","FABP5","SLCO2A1","KCNIP4","CD74","HLA-DRA","PLAT","TFF3","HLA-DRB1","FABP4","SELE","MT2A","LYZ","FLRT2","ANO2","ADIRF","ZNF385D","NR2F2","IFI30","S100A10","S100A6","GPM6A","NRG3","VCAM1","HLA-DQB1","NNMT","CEBPD")
png('gene.test.png',width = 1200,height = 1500, res = 100)
FeaturePlot(All_cell_post_harmony, features = gene_test,raster=FALSE)
dev.off()
All_cell_post_harmony$gene_test_mean = colMeans(as.matrix(All_cell_post_harmony@assays[['RNA']]@data[gene_test,]))

png('gene.test2.png',width = 500,height = 400, res = 100)
FeaturePlot(All_cell_post_harmony, features = 'gene_test_mean',raster=FALSE)
dev.off()


dir.create('Test')
gene_test = read.csv('/Volumes/Project_ESCC/DATA/endo_c5_sigs_1.csv')
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

# --------------------------------------------------------------------------------
# Figure 1B

P1 = DimPlot(Tcell_srt, reduction = "umap",label = F, pt.size = 0.4, group.by = "cluster_name", raster=FALSE,cols = zylcolor40) + NoLegend()
P2 = DimPlot(Tcell_srt, reduction = "umap",label = F, pt.size = 0.4, group.by = "response",raster=FALSE) +  scale_color_manual(values = color_palette_Response) +  NoLegend()

pdf(paste0(output_path,'/Figure1B.pdf'), width = 18, height=7)
P1+P2
dev.off()

# --------------------------------------------------------------------------------
# Figure 1C 
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
pdf('Figure1/Fig1C.pdf',width = 18,height = 2)
p1%>%insert_top(p2,height = 0.5)
dev.off()

# ----------------------------------------------------------------------------- 
####Figure1D#######
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
# Figure1E
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

###BULK####
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
TRAILmark <- read.gmt('./Final_RDS/TRAIL.gmt')
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


pdf(paste0(output_path,'/Figure1E.pdf'),width = 3.5,height = 4)
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
# Figure 1F
Tcell_srt = readRDS('./Final_RDS/Tcell_withTCR.rds')
metadata = Tcell_srt@meta.data 
metadata = metadata[!is.na(metadata$raw_clonotype_id), ]
metadata$Clono_cell_num = cut(metadata$clono_cell_num2, c(0,1,2,3,4,5,6,Inf),right = FALSE, labels = c('0','1','2', '3', '4','5','>=6'))
clonal_cell = metadata %>% group_by(response, Clono_cell_num) %>% dplyr::summarise(clono_num=n())%>% as.data.frame()

clonal_cell[which(clonal_cell$response == 'Responder'),'response'] <- 'R'
clonal_cell[which(clonal_cell$response == 'Non_responder'),'response'] <- 'NR'

clonal_cell$response = factor(clonal_cell$response, levels=c('NR', 'R'))

pdf(paste0(output_path, '/Figure1F.pdf'),width = 3,height = 4)
ggplot(data = clonal_cell, aes(x=response, y=clono_num, fill=Clono_cell_num)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  scale_fill_manual(values= color_palette_clone )+
  labs(y='Cell proportion', x='Group')+
  guides(fill=guide_legend(title='Clone size'))+
  theme_bw()+
  theme(
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"),
    panel.border = element_rect(fill=NA,color="black",linetype="solid"))
dev.off()

#-----------------------------------------------------------
# Figure 1G
cyto_res <- t(read.table('./Final_RDS/Cytoseq_output/CloneTCR/output_expand.Zscore', sep='\t'))
cyto_res <- as.data.frame(cyto_res[, c('TGFB1', 'PGE2', 'TRAIL')])
cyto_res$Clono_cell_num = sapply(rownames(cyto_res), function(x)strsplit(x, '[\\.]+')[[1]][2])
line_data <- cyto_res %>% group_by(Clono_cell_num) %>% dplyr::summarise(TRAIL = median(TRAIL)) %>% as.data.frame()
Figure2B <- cyto_res %>%mutate(Clono_cell_num = factor(sapply(Clono_cell_num, function(x) ifelse(x==6, '>=6', as.character(x))), levels = c('0','1','2', '3', '4','5','>=6')))  

pdf(paste0(output_path, '/Figure1G.pdf'),width = 2.8,height = 4)
ggplot(data = Figure2B,aes(x=Clono_cell_num, y=TRAIL))+
  stat_boxplot(Fgeom = 'errorbar',width = 0.4) +
  geom_boxplot(outlier.fill="white",outlier.color="white",width = 0.8) +
  geom_dotplot(binaxis='y', aes(fill = Clono_cell_num), stackdir='center', stackratio=1, dotsize=0.7, binwidth = 0.6)+
  scale_fill_manual(values=color_palette_clone) +
  scale_color_manual(values = color_palette_clone)+
  geom_smooth(data=line_data, mapping=aes(x=as.numeric(Clono_cell_num), y=TRAIL,linetype = 'dashed'),
              method = "lm", se=F,  formula = y~x, size=1.2, inherit.aes = F, color = "#FF0000")+
  stat_cor(data=line_data, mapping = aes(x = as.numeric(Clono_cell_num), y= TRAIL), method='pearson', size = 4, color='black', inherit.aes = F)+
  xlab('Clone Size') + 
  ylab( 'TRAIL Signaling') +
  theme_few( ) +
  theme(
    legend.position = 'no',
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"))
dev.off() 

#-----------------------------------------------------------
# Figure 1H
All_cell<-readRDS("./Final_RDS/All_cells_after_harmony.rds")
All_cell$cell.type[is.na(All_cell$cell.type)]='Tumor cell'
All_cell@meta.data[All_cell@meta.data[, 'cell.type'] == 'NA',  'cell.type'] = 'Tumor cell'
All_cell@meta.data[All_cell@meta.data[, 'cell.type'] %in% c("Tumor or epithelial cells"),  'cell.type'] = "Tumor cell"

IL7R_exp = FetchData(All_cell, vars="TNFSF10")
IL7R_exp$cell_type = All_cell@meta.data[rownames(IL7R_exp), "cell.type"]
IL7R_exp$sample = All_cell@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = All_cell@meta.data[rownames(IL7R_exp), "response"]

IL7R_exp_mean = IL7R_exp %>%
  group_by(cell_type,sample,response) %>%
  summarise_all(list(mean)) %>% as.data.frame()
colnames(IL7R_exp_mean)<-c("cell_type","sample","response","TRAIL_expression")

# ----------------------------------------------------------------------------
IL7R_exp_mean$Group <- 'b'
Group1 <- c('Mast cell', 'Endothelial cell',"Tumor cell","Macrophage/Monocyte")
IL7R_exp_mean[IL7R_exp_mean$cell_type%in%Group1,'Group'] <- 'a'


my_comparisons <- list(c('a','b'))
library(rstatix)

P_sig <- IL7R_exp_mean %>%  
  wilcox_test(TRAIL_expression ~ Group)%>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "cell_type", dodge = 0.8) 


pdf(paste0(output_path, '/Figure1H.pdf'),width = 4,height = 4)
ggplot(IL7R_exp_mean, aes(x=reorder(cell_type,-TRAIL_expression),y=TRAIL_expression)) +
  stat_boxplot(geom = 'errorbar',width = 0.5) +
  geom_boxplot(outlier.fill="white",outlier.color="white") +
  geom_jitter( width = 0.1,aes(color = cell_type)) +
  annotate('text',x=6, y = 2.5, label=paste0('P = ',P_sig$p.adj))+
  # stat_compare_means(data =IL7R_exp_mean ,mapping = aes(x = Group,TRAIL_expression),comparisons = my_comparisons,
  #                    method = "wilcox.test",size=4,label = "p.format",angle=0,vjust = 0, hjust=0)+
  scale_color_manual(values = color_cell_type) +
  theme_cowplot(font_size = 12) +
  NoLegend() +
  theme(
    axis.text.x = element_text(colour="black",angle=30,hjust = 1,vjust = 1), 
    axis.line = element_line(colour = "black"), 
    legend.text = element_text(  colour="black"),
    legend.title = element_text( colour="black"))
dev.off() 

#-----------------------------------------------------------
# Figure 1I
Tumor <- subset(All_cell,cell.type=="Tumor cell")
TNFSF10_exp = FetchData(Tumor, vars="TNFSF10")
TNFSF10_exp$orig.ident = Tumor@meta.data[rownames(TNFSF10_exp), "orig.ident"]
TNFSF10_exp$response = Tumor@meta.data[rownames(TNFSF10_exp), "response"]
TNFSF10_exp_Tumor = TNFSF10_exp %>%
  group_by(orig.ident,response) %>%
  dplyr::summarise_all(list(mean)) %>% as.data.frame()
TNFSF10_exp_Tumor

output_dir = './Final_RDS//Cytoseq_output/new_updated/' # output path

cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample_new<-cyto_res[,"TRAIL"]
TRAIL_sample_new<-as.data.frame(TRAIL_sample_new)
plot_counts_Tumor_sample<-cbind(TNFSF10_exp_Tumor,TRAIL_sample_new)
plot_counts_Tumor_sample_no_16_pre<-subset(plot_counts_Tumor_sample,orig.ident!="EC_16_pre_T")
plot_counts_Tumor_sample_no_16_pre<-subset(plot_counts_Tumor_sample_no_16_pre,orig.ident!="EC_07_post_T")

Figure1I <- plot_counts_Tumor_sample_no_16_pre
Figure1I[which(Figure1I$response=='Responder'),'response'] <- 'R'
Figure1I[which(Figure1I$response=='Non_responder'),'response'] <- 'NR'


pdf(paste0(output_path, '/Figure1I.pdf'),width = 3,height = 4)
ggplot(Figure1I,mapping =aes_string(x='TNFSF10', y='TRAIL_sample_new'),fill="response") +
  geom_point(aes_string(color="response"),size=2) + 
  geom_smooth(method = "lm", se=F, formula = y~x, size=1,color = 'red')+
  stat_cor(method='pearson', size=5) +
  scale_color_manual(values = color_palette_Response) +
  xlab('TRAIL expression of tumor cells') +
  ylab('TRAIL activity of T cells') +
  theme_bw() +
  theme(
    legend.position = "no",
    panel.border = element_rect(color = "black",size = 1, fill=NA),
    axis.text.x = element_markdown(colour = 'black'),
    axis.text.y = element_markdown(colour = 'black'),
    legend.text = element_text(size = 8),
    panel.grid = element_blank())
dev.off() 

#' ---------------------------------------------------------
#' Figure1J
Tumor <- readRDS("./Final_RDS/All_tumor_3.22.rds")
Tumor_pre <- subset(Tumor,time=="Pre_T")
Tumor_post <- subset(Tumor,time=="Post_T")

IL7R_exp = FetchData(Tumor, vars="TNFSF10")
IL7R_exp$sample = Tumor@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = Tumor@meta.data[rownames(IL7R_exp), "response"]
IL7R_exp$time = Tumor@meta.data[rownames(IL7R_exp), "time"]
IL7R_exp_mean = IL7R_exp %>%
  group_by(sample,response,time) %>%
  dplyr::summarise_all(list(mean)) %>% as.data.frame()

scRNA_seq <- IL7R_exp_mean
scRNA_seq$Title <- 'scRNA-seq'

# bulk
sample_info = read.table('ESCC_bulk_final/sample_info.txt', header = F, sep='\t')
colnames(sample_info) = c('sample', 'time', 'response')
sample_info$sample = gsub('-', '_', sample_info$sample)
rownames(sample_info) = sample_info$sample
sample_info[sample_info[,"response"] %in% c("CR","CR "),  'response'] = 'pCR'

sample_info[sample_info[,"response"] %in% c("pCR","MPR"),  'response2'] = 'Responder'
sample_info[sample_info[,"response"] %in% c("PR","SD","PD","NA"),  'response2'] = 'Non_responder'
sample_info[sample_info[,"response"] %in% c('pCR'),  'response3'] = 'pCR'
sample_info[sample_info[,"response"] %in% c("MPR","PR","SD","PD","NA"),  'response3'] = 'Non_pCR'
sample_info[sample_info[,"response"] %in% c('pCR',"MPR","PR","SD"),  'response4'] = 'Non_PD'
sample_info[sample_info[,"response"] %in% c("PD","NA"),  'response4'] = 'PD'
sample_info[sample_info[,"response"] %in% c('pCR',"MPR"),  'response5'] = 'pCR+MPR'
sample_info[sample_info[,"response"] %in% c("PR","SD","PD","NA"),  'response5'] = 'PR+PD'

table(sample_info$response4)
dim(sample_info)
sample_info2<-subset(sample_info,response2%in%c('Responder','Non_responder'))
#sample_info_Pre<-subset(sample_info2,time%in%c('Pre-T'))
sample_info_Pre_Post<-subset(sample_info,time%in%c('Pre-T',"Post-T"))
sample_info_Pre_Post2<-subset(sample_info_Pre_Post,sample!='ESCA_16_pre_T')
sample_info_Pre_Post2<-subset(sample_info_Pre_Post2,sample!='ESCA_16_post_T')
sample_info_Pre<-subset(sample_info_Pre_Post2,time=='Pre-T')
sample_info_Post<-subset(sample_info_Pre_Post2,time=='Post-T')
#bulk_tpm = read.table('/Users/lab5/Desktop/All samples_finally/RNA-seq/ESCC_bulk_final/EC_bulk_tmp.txt', header = T, sep='\t')
bulk_tpm = read.table('ESCC_bulk_final/EC_bulk_count.txt', header = T, sep='\t')
#sample_info_Pre[sample_info_Pre[,"sample"] %in% c("EC_45_pre_T","ESCA_25_pre_T","ESCA_01_pre_T"),  'response2'] = 'Non_responder'
#sample_info_Pre$response2=factor(sample_info_Pre$response2,levels = c("Responder","Non_responder"))
my_comparisons <- list(c("Responder","Non_responder"))
###all samples
m6A_genes=c("TNFSF10")
data <- bulk_tpm[m6A_genes, ] %>%
  log1p()%>%
  mutate(gene=m6A_genes) %>%
  reshape2::melt() %>%
  merge(sample_info2, by.x='variable', by.y='sample')

bulk <- data[,c('response2','value')]
bulk$Title <- "Bulk RNA-seq"
colnames(bulk) <- c('Response','TRAIL_expression','Title')

scRNA_seq <- scRNA_seq[,c( 'response','TNFSF10','Title')]
colnames(scRNA_seq) <- c('Response','TRAIL_expression','Title')

Figure1J <- rbind(scRNA_seq,bulk)
Figure1J$Title <- factor(Figure1J$Title,levels = c('scRNA-seq','Bulk RNA-seq'))

head(Figure1J)
Figure1J[which(Figure1J$Response=='Responder'),'Response'] <- 'R'
Figure1J[which(Figure1J$Response=='Non_responder'),'Response'] <- 'NR'

pdf(paste0(output_path, '/Figure1J.pdf'),width = 4,height = 4)
ggplot(Figure1J, aes(x = Response, y = TRAIL_expression,color=Response)) +
  stat_boxplot(geom = 'errorbar',width = 0.4 ) +
  geom_boxplot(width=0.5, size=0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1,size = 1.5) +
  stat_compare_means(comparisons = list(c('R', 'NR')), label = 'p.format',face="bold") +
  facet_wrap(~Title,scales = 'free') + 
  scale_color_manual(values =color_palette_Response) +
  theme_bw()+
  xlab('') +
  ylab('TRAIL expression')+
  theme(
    legend.position = "top",
    strip.text.y = element_text( size = 14),
    strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
    panel.border = element_rect(color = "black",size = 1, fill=NA),
    axis.text.x = element_markdown(colour = 'black'),
    axis.text.y = element_markdown(colour = 'black'),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()) + 
  NoLegend() 
dev.off()


#' ---------------------------------------------------------
#' Figure1K
#' correlation of DR-5 expression in T cells and TRAIL signaling
TRAIL_receptor<-c("TNFRSF10B")
IL7R_exp = FetchData(T_cell_no_NK, vars=TRAIL_receptor)
IL7R_exp$orig.ident = T_cell_no_NK@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = T_cell_no_NK@meta.data[rownames(IL7R_exp), "response"]
IL7R_exp$time = T_cell_no_NK@meta.data[rownames(IL7R_exp), "time"]
IL7R_exp_mean = IL7R_exp %>%
  group_by(orig.ident,response,time) %>%
  summarise_all(list(mean)) %>% as.data.frame()
output_dir = './Final_RDS/Cytoseq_output/new_updated/' # output path

# bulk
cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample_new<-cyto_res[,"TRAIL"]
TRAIL_sample_new<-as.data.frame(TRAIL_sample_new)
data<-cbind(IL7R_exp_mean,TRAIL_sample_new)
data<-melt(data, id=c("orig.ident","response","time","TRAIL_sample_new"))
colnames(data)[colnames(data)=="variable"]<- "TRAIL_receptor"
colnames(data)[colnames(data)=="value"]<- "TRAIL_receptor_exp"
data_no_16_pre <- subset(data,orig.ident!="EC_16_pre_T")
data_no_16_pre <- subset(data_no_16_pre,orig.ident!="EC_07_post_T")

Figure1K <- data_no_16_pre
Figure1K[which(Figure1K$response=='Responder'),'response'] <- 'R'
Figure1K[which(Figure1K$response=='Non_responder'),'response'] <- 'NR'
head(Figure1K)

pdf(paste0(output_path, '/Figure1K.pdf'),width = 3,height = 4)
ggplot(Figure1K,mapping =aes_string(x='TRAIL_receptor_exp', y='TRAIL_sample_new'),fill="response") +
  geom_point(aes_string(color="response"),size=2) + 
  geom_smooth(method = "lm", se=F, formula = y~x, size=1,color = 'red') +
  stat_cor(method='spearman', size=4)+
  scale_color_manual(values = color_palette_Response) +
  xlab('Average expression in T cells') +
  ylab('TRAIL signaling activity of T cells') +
  ggtitle('DR-5')+
  theme_bw() +
  theme(
    legend.position = "no",
    panel.border = element_rect(color = "black",size = 1, fill=NA),
    axis.text.x = element_markdown(colour = 'black'),
    axis.text.y = element_markdown(colour = 'black'),
    legend.text = element_text(size = 8),
    panel.grid = element_blank())
dev.off() 
