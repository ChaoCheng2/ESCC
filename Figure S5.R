
#-------------------------------------------------------------------------------
#Figure S5
# Setwd path
output_path = 'Figure_output/FigureS5_original'
dir.create(file.path(workdir,output_path))

#--------------------------------------------------------
color_palette_Response <- c('NR' = "#0099CCFF",'R' = "#FFC20AFF")
color_palette_clone <- c('#F7CCBF','#F19386', '#F06976','#C76578','#975C75',"#8E3F65")
color_cell_type = c("B cell" = "#E64B35FF", "Endothelial cell"= "#4DBBD5FF", "Fibroblast" = "#00A087FF", "Macrophage/Monocyte" = "#3C5488FF", 
                    "Mast cell" = "#F39B7FFF", "Neutrophils"= "#8491B4FF", "Plasma" = "#91D1C2FF",  "T cell" = "#fdbf6f", "Tumor or epithelial cells" = "#7E6148FF")
color_palette_Treatment = c('Post_T'='#91d1c2ff', 'Pre_T'= '#f768a1')

# ----------------------------------------------------------
# Source
library(tidyr)
library(ggpubr)
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggsci)
library(cowplot)
library(egg)
library(reshape2)
library(RColorBrewer)
library(SingleR)
library(SingleCellExperiment)
library(pheatmap)
library(harmony) 

#' Figure S5A
#' Comparing TRAIL signaling in different sub cell type
#' CD4 T
T_cell_no_NK <- readRDS("./Final_RDS/T_cell_no_NK.rds")
CD4_T<-subset(T_cell_no_NK,T_subtype=="CD4_T")
data<- FetchData(CD4_T,vars = c("cluster_name","TRAIL_activity_new"))
my_comparisons <- list(c("CD4_C1_CCR7","CD4_C3_IL17A"))
data2 <- data %>% group_by(cluster_name) %>% summarise_all(list(mean)) %>% as.data.frame()
data3 <- merge(data2,data,by="cluster_name")
colnames(data3) <- c("cluster_name","mean","TRAIL_Signaling")
Figure2C_CD4 <- data3
Figure2C_CD4$Title <- 'CD4'

#CD8 T
CD8_T<-subset(T_cell_no_NK,T_subtype=="CD8_T")
data<- FetchData(CD8_T,vars = c("cluster_name","TRAIL_activity_new"))
my_comparisons=list(c("CD8_C3_ANXA1","CD8_C5_CX3CR1"))
data2 <- data %>% group_by(cluster_name) %>% summarise_all(list(mean)) %>% as.data.frame()
data3 <- merge(data2,data,by="cluster_name")
colnames(data3) <- c("cluster_name","mean","TRAIL_Signaling")
Figure2C_CD8 <- data3
Figure2C_CD8$Title <- 'CD8'
Figure2C_data <- rbind(Figure2C_CD4,Figure2C_CD8)
my_comparisons = list(c("CD4_C1_CCR7","CD4_C3_IL17A"),c('CD8_C3_ANXA1','CD8_C5_CX3CR1'))

library(ggpubr)
library(rstatix)

sig <- Figure2C_CD4[Figure2C_CD4$cluster_name%in%c("CD4_C1_CCR7","CD4_C3_IL17A"),]

x = round(wilcox.test(TRAIL_Signaling~cluster_name,data = sig)$p.value,5)

sig_x <- data.frame(cluster_name = 'CD4_C1_CCR7', y = 9,Title = 'CD4',pvalue = paste0('P = ',x))
sig_x

pdf(paste0(output_path, '/Figure2C.pdf'),width = 4,height = 4)
ggplot(Figure2C_data,aes(x=reorder(cluster_name,-TRAIL_Signaling), y = TRAIL_Signaling, fill = mean)) +
  geom_violin(width = 0.7,size = 0.1) +
  geom_boxplot(width = 0.1,color = '#525252',alpha = 0.3,outlier.size = 0.4,size = 0.1)+
  #stat_compare_means( comparisons = list(c("CD4_C1_CCR7","CD4_C3_IL17A")), method = "wilcox.test",label = "p.format")+
  geom_text(aes(x = cluster_name,y = y,label = pvalue),data = sig_x,inherit.aes = F) +
  stat_compare_means( comparisons = list(c('CD8_C3_ANXA1','CD8_C5_CX3CR1')), method = "wilcox.test",label = "p.format")+
  facet_grid(~Title,scales = 'free_x',space = 'free') +
  scale_fill_gradientn( colours = c( "#F8B195","#F67280","#C06C84","#6C5B7B")) +
  theme_few() +
  xlab("") +
  scale_y_continuous(position = 'left') +
  theme(
    legend.position = 'no',
    panel.spacing.x = unit(0, "cm"),
    strip.background = element_rect(color="black", fill="white",linetype="solid" ),
    panel.border = element_rect(color = "black",  fill=NA),
    axis.text.x = element_markdown(colour = 'black', angle = 30,hjust = 1,vjust = 1),
    axis.text.y = element_markdown(colour = 'black'),
    legend.text = element_text(),
    panel.grid = element_blank())
dev.off()


# -------------------------------------------------------------------------------
#' Figure S5B 
#' TRAIL activity CD4-CCR7
CD4_CCR7_Pre <- subset(T_cell_no_NK,cluster_name=="CD4_C1_CCR7"&time=="Pre_T")
data<- FetchData(CD4_CCR7_Pre,vars = c("response","TRAIL_activity_new"))
my_comparisons=list(c("Responder","Non_responder"))
data2 = data %>%
  group_by(response) %>%
  summarise_all(list(mean)) %>% as.data.frame()
data3 <- merge(data2,data,by="response")
colnames(data3) <- c("response","mean","TRAIL_Signaling")

Figure2D_CD4 <- data3
Figure2D_CD4$Title <- 'CD4_C1_CCR7'

CD8_ANXA1_Pre<-subset(T_cell_no_NK,cluster_name=="CD8_C3_ANXA1"&time=="Pre_T")
data<- FetchData(CD8_ANXA1_Pre,vars = c("response","TRAIL_activity_new"))
my_comparisons=list(c("Responder","Non_responder"))
data2 = data %>%
  group_by(response) %>%
  summarise_all(list(mean)) %>% as.data.frame()
data3<-merge(data2,data,by="response")
colnames(data3)<-c("response","mean","TRAIL_Signaling")

Figure2D_CD8 <- data3
Figure2D_CD8$Title <- 'CD8_C3_ANXA1'


FigureS5B_data <- rbind(Figure2D_CD4,Figure2D_CD8)
FigureS5B_data$cluster_name = FigureS5B_data$response


pdf(paste0(output_path, '/FigureS5B.pdf'),width = 2.8 ,height = 4)
ggplot(FigureS5B_data,aes(x= response, y = TRAIL_Signaling, fill = mean)) +
  geom_violin(trim = FALSE,width = 0.7, position = position_dodge(0.8),size = 0.1) +
  geom_boxplot(width = 0.1,color = '#525252',outlier.size = NULL,alpha = 0.8,size = 0.1, position = position_dodge(0.8) ) +
  stat_compare_means(ref.group = NULL,method = "wilcox.test",size=3,label = "p.format",paired = FALSE) +
  scale_fill_gradientn( colours = c( "#F8B195","#F67280","#C06C84","#6C5B7B")) +
  facet_wrap(~Title) +
  theme_bw() +
  xlab("") +
  scale_y_continuous(position = 'left') +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "cm"),
    strip.text.y = element_text( size = 10),
    strip.background = element_rect(color="black", fill="white", size=1, linetype="solid" ),
    panel.border = element_rect(color = "black",size = 1, fill=NA),
    axis.text.x = element_markdown(colour = 'black'),
    axis.text.y = element_markdown(colour = 'black'),
    legend.text = element_text(size = 8),panel.grid = element_blank())
dev.off()
FigureS5B_data_c <- FigureS5B_data[,c("cluster_name","mean","TRAIL_Signaling","Title")]
#' FigureS5C
All_cell<-readRDS("./Final_RDS/All_cells_after_harmony_10.13.rds")
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


pdf(paste0(output_path, '/FigureS5C.pdf'),width = 4,height = 4)
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

# Figure S5D  correlation between  TRAIL Expr in other cell and TRAIL activity of T cell 

cellType_Data <- subset(All_cell,cell.type=="Mast cell" | cell.type=="Endothelial cell"| cell.type=="Neutrophils"); head(cellType_Data)

TNFSF10_exp = FetchData(cellType_Data, vars="TNFSF10")
TNFSF10_exp$sample = cellType_Data@meta.data[rownames(TNFSF10_exp), "orig.ident"]
TNFSF10_exp$response = cellType_Data@meta.data[rownames(TNFSF10_exp), "response"]
TNFSF10_exp$Celltype = cellType_Data@meta.data[rownames(TNFSF10_exp), "cell.type"]

TNFSF10_exp = TNFSF10_exp %>%
  group_by(sample,response,Celltype) %>%
  dplyr::summarise_all(list(mean)) %>% as.data.frame()

output_dir = './Final_RDS//Cytoseq_output/new_updated/' # output path

cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample_new <-  data.frame(sample = rownames(cyto_res),TRAIL_signaling = cyto_res[,"TRAIL"])
plot_counts_Tumor_sample <- merge(TNFSF10_exp,TRAIL_sample_new,by ="sample" )
Plot_df <- subset(plot_counts_Tumor_sample, sample!="EC_16_pre_T")
Plot_df <- subset(Plot_df, sample!="EC_07_post_T")
Plot_df$Celltype <- factor(Plot_df$Celltype,levels = c("Mast cell","Endothelial cell","Neutrophils"))


pdf(file = 'Supple_Figure/Supp_Figure5/FigS5D.pdf',width = 6.5,height = 3)
ggplot(Plot_df, mapping =aes_string(x = 'TNFSF10', y='TRAIL_signaling'),fill="response") +
  geom_point(aes_string(color="response"),size=2) + 
  geom_smooth(method = "lm", se=F, formula = y~x, size=1,color = 'red')+
  stat_cor(method='pearson')+
  scale_color_manual(values = color_palette_Response)+
  facet_wrap(~Celltype,scales = 'free_x') +
  xlab("TRAIL expression")+
  ylab("TRAIL signaling of T cells") +
  theme_classic() +
  theme(axis.text = element_text(colour="black"))
dev.off()      


# ------------------------------------------------------------------------------
# Figure S5E

T_cell_no_NK <- readRDS("./Final_RDS/T_cell_no_NK_updated_8.8.rds")
TRAIL_receptor <- c("TNFRSF10A","TNFRSF10B","TNFRSF10C","TNFRSF11B")

IL7R_exp = FetchData(T_cell_no_NK, vars=TRAIL_receptor)
IL7R_exp$orig.ident = T_cell_no_NK@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = T_cell_no_NK@meta.data[rownames(IL7R_exp), "response"]
IL7R_exp$time = T_cell_no_NK@meta.data[rownames(IL7R_exp), "time"]
IL7R_exp_mean = IL7R_exp %>%
  group_by(orig.ident,response,time) %>%
  summarise_all(list(mean)) %>% as.data.frame()

output_dir = './Final_RDS/Cytoseq_output/new_updated/' # output path

cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample_new <- cyto_res[,"TRAIL"]
TRAIL_sample_new <- as.data.frame(TRAIL_sample_new)
data <- cbind(IL7R_exp_mean,TRAIL_sample_new)
data <- melt(data, id=c("orig.ident","response","time","TRAIL_sample_new"))
colnames(data)[colnames(data)=="variable"]<- "TRAIL_receptor"
colnames(data)[colnames(data)=="value"]<- "TRAIL_receptor_exp"
data_no_16_pre<-subset(data,orig.ident!="EC_16_pre_T")
data_no_16_pre<-subset(data_no_16_pre,orig.ident!="EC_07_post_T")

pdf(file = 'Supple_Figure/Supp_Figure3/FigS3C.pdf',width = 7,height = 6)
ggplot(data_no_16_pre,mapping =aes_string(x='TRAIL_receptor_exp', y='TRAIL_sample_new'),fill="response") +
  geom_point(aes_string(color="response"),size=2) + 
  scale_color_manual(values = color_palette_Response)+
  stat_cor(method='spearman')+
  geom_smooth(method = "lm", se=F, formula = y~x)+
  facet_wrap(~TRAIL_receptor, ncol = 2, scale='free')+
  xlab("TNFRSF10B expression of T cells")+
  ylab("TRAIL signaling of T cells") +
  theme_classic() +
  theme(axis.text = element_text(colour="black")) 
dev.off()

# ---------------------------------------------------------------------------------
# Figure S5G

Tcell_srt = readRDS('./Final_RDS/Tcell_withTCR.rds')

ccgene_zzm = c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
               'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
               'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1')

Tcell_srt = AddModuleScore(Tcell_srt, features = list(ccgene_zzm), name='zzm')
#
CD4_srt = subset(Tcell_srt, cluster_name%in%grep('CD4', unique(Tcell_srt$cluster_name), value=T))
cc_exp2 = CD4_srt@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarise(zzm=median(zzm1)) %>% as.data.frame()
cc_exp = cc_exp2$zzm
names(cc_exp) = cc_exp2$orig.ident
target_gene = c('CXCL2', 'CCL20', 'CXCL8', 'TNFAIP3', 'BIRC3', 'IRF1', 'CXCL1', 'NFKBIA', 'FOS', 'JUNB', 'PPP1R15A', 'EGR1', 'PLK3', 'RHOB')
trail_target_exp_CD4 = AverageExpression(CD4_srt, features = target_gene, group.by = 'orig.ident')$RNA %>% log1p()

cor_data_CD4 = c()
cor_data_CD4_pvale = c()

for(i in 1:nrow(trail_target_exp_CD4)){
  ct = cor.test(trail_target_exp_CD4[i,], cc_exp[colnames(trail_target_exp_CD4)])
  cor_data_CD4 = rbind(cor_data_CD4, ct$estimate)
  cor_data_CD4_pvale = rbind(cor_data_CD4_pvale, ct$p.value)
}
#
CD8_srt = subset(Tcell_srt, cluster_name%in%grep('CD8', unique(Tcell_srt$cluster_name), value=T))
cc_exp2 = CD8_srt@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarise(zzm=median(zzm1)) %>% as.data.frame()
cc_exp = cc_exp2$zzm
names(cc_exp) = cc_exp2$orig.ident
target_gene = c('CXCL2', 'CCL20', 'CXCL8', 'TNFAIP3', 'BIRC3', 'IRF1', 'CXCL1', 'NFKBIA', 'FOS', 'JUNB', 'PPP1R15A', 'EGR1', 'PLK3', 'RHOB')
trail_target_exp_CD8 = AverageExpression(CD8_srt, features = target_gene, group.by = 'orig.ident')$RNA %>% log1p()
cor_data_CD8 = cor(t(trail_target_exp_CD8), cc_exp[colnames(trail_target_exp_CD8)])

cor_data_CD8 = c()
cor_data_CD8_pvale = c()
for(i in 1:nrow(trail_target_exp_CD8)){
  ct = cor.test(trail_target_exp_CD8[i,], cc_exp[colnames(trail_target_exp_CD8)])
  cor_data_CD8 = rbind(cor_data_CD8, ct$estimate)
  cor_data_CD8_pvale = rbind(cor_data_CD8_pvale, ct$p.value)
}


cor_data = as.data.frame(cbind(cor_data_CD4, cor_data_CD8))
cor_data_pvalue = as.data.frame(cbind(cor_data_CD4_pvale, cor_data_CD8_pvale))
colnames(cor_data) = c('CD4', 'CD8')
rownames(cor_data) = target_gene

pos1 = cor_data_pvalue<0.01
pos2 = cor_data_pvalue<0.05
pos3 = cor_data_pvalue>=0.05
cor_data_pvalue[pos1] = '**'
cor_data_pvalue[pos2] = '*'
cor_data_pvalue[pos3] = ''

width_in = 50 / 25.4
height_in = 100 / 25.4

pdf(file = 'Supple_Figure/Supp_Figure5/FigS5G.pdf',width = 8.5,height = 3)
bk <- c(seq(-0.5,-0.1,by=0.01),seq(0,0.5,by=0.01))
pheatmap::pheatmap(cor_data,
                   cluster_rows = T, cluster_cols = F,
                   treeheight_row = 5,
                   show_column_names = T,
                   display_numbers = cor_data_pvalue,
                   fontsize_number = 15,
                   number_color = 'white',
                   border_color = 'white',
                   cellwidth  = 16,
                   cellheight = 16,
                   color = c(colorRampPalette(c("#053061","#2166ac","#4393c3","white"))(length(bk)/2),colorRampPalette(c("white","#fddbc7","#d6604d","#b2182b","#67001f"))(length(bk)/2)),
                   breaks=bk
)
dev.off()

# ---------------------------------------------------------------------------------
# Figure S5H
CDT_cell_seurat <- subset(T_cell_no_NK,T_subtype=="CD4_T" | T_subtype=="CD8_T")

data <- FetchData(CDT_cell_seurat,vars =c("NFKBIA","EGR1","TRAIL_activity_new"))

data$sample = CDT_cell_seurat@meta.data[rownames(data), "orig.ident"]
data$response = CDT_cell_seurat@meta.data[rownames(data), "response"]
data$T_subtype = CDT_cell_seurat@meta.data[rownames(data), "T_subtype"]

data2 = data %>%
  group_by(sample,response,T_subtype) %>%
  summarise_all(list(mean)) %>% as.data.frame()

head(data2)
library(tidyr)
data_long <- gather(data2, key = 'gene',value = 'Expr', -sample,-response,-T_subtype,-TRAIL_activity_new )
head(data_long)

data_long$gene <- factor(data_long$gene,levels = c("NFKBIA","EGR1"))

pdf(file = 'Supple_Figure/Supp_Figure5/FigS5H.pdf',width = 10,height = 3)
ggplot(data = data_long, mapping =aes_string(x = 'Expr', y = 'TRAIL_activity_new'),fill="response") + 
  geom_point(aes_string(color="response"),size=2) + 
  geom_smooth(method = "lm", se=F, formula = y~x, size=1,color = 'red')+
  stat_cor(method='pearson', size=5)+
  scale_color_manual(values = color_palette_Response) +
  facet_wrap(~T_subtype+gene,  scale='free',ncol = 4)+
  theme_classic()

dev.off()

# ---------------------------------------------------------------------------------
# Figure S5I
NFKBIA <- FetchData(T_cell_no_NK,vars =c("NFKBIA","TRAIL_activity_new"))
NFKBIA$sample = T_cell_no_NK@meta.data[rownames(NFKBIA), "orig.ident"]
NFKBIA$response = T_cell_no_NK@meta.data[rownames(NFKBIA), "response"]
NFKBIA = NFKBIA %>%
  group_by(sample,response) %>%
  summarise_all(list(mean)) %>% as.data.frame()

colnames(NFKBIA) <-c('Sample','Response','Expression','TRAIL_Signaling')
NFKBIA$Gene <- 'NFKBIA'

EGR1 <- FetchData(T_cell_no_NK,vars =c("EGR1","TRAIL_activity_new"))
EGR1$sample <- T_cell_no_NK@meta.data[rownames(EGR1), "orig.ident"]
EGR1$response <- T_cell_no_NK@meta.data[rownames(EGR1), "response"]
EGR1 <- EGR1%>%
  group_by(sample,response) %>%
  summarise_all(list(mean)) %>% as.data.frame()

colnames(EGR1) <-c('Sample','Response','Expression','TRAIL_Signaling')
EGR1$Gene <- 'EGR1'
Figure2J_data <- rbind(NFKBIA,EGR1)

Figure2J_data[which(Figure2J_data$Response=='Responder'),'Response'] <- 'R'
Figure2J_data[which(Figure2J_data$Response=='Non_responder'),'Response'] <- 'NR'

Figure2J_data$Gene <- factor(Figure2J_data$Gene,levels = c('NFKBIA','EGR1'))

pdf(paste0(output_path, '/FigureS5I.pdf'),width = 5.5,height = 4)
ggplot(Figure2J_data , mapping = aes_string(x = 'Expression', y = 'TRAIL_Signaling'),fill = "Response") +
  geom_point(aes_string(color="Response"),size=2) + 
  scale_color_manual(values = color_palette_Response)+
  geom_smooth(method = "lm", se=F, formula = y~x, size=1,col = 'red')+
  stat_cor(method='pearson')+
  facet_wrap(~Gene,scales = 'free')+
  theme_classic()+
  xlab("Average expression in T cells")+
  ylab("TRAIL signaling activity of T cells") +
  theme(
    legend.position = "no",
    strip.text.y = element_text( size = 14),
    strip.background = element_rect(fill="white"),
    panel.grid = element_blank())
dev.off()   

# ---------------------------------------------------------------------------------
# Figure S5J
IL7R_exp = FetchData(T_cell_no_NK, vars=c("EGR1","NFKBIA"))
IL7R_exp$sample = T_cell_no_NK@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$response = T_cell_no_NK@meta.data[rownames(IL7R_exp), "response"]

result = c()
row_name = c()
for(s in unique(IL7R_exp$sample)){
  dp = 0
  dn = 0
  sp_il17r = 0
  sp_cd48 = 0
  
  tmp_data = subset(IL7R_exp, sample==s)
  for(i in 1:nrow(tmp_data)){
    if(tmp_data[i, 'EGR1'] > 0 && tmp_data[i, 'NFKBIA'] > 0){
      dp=dp+1
    }else if(tmp_data[i, 'EGR1'] > 0 && tmp_data[i, 'NFKBIA'] == 0){
      sp_il17r=sp_il17r + 1
    }else if(tmp_data[i, 'EGR1'] == 0 && tmp_data[i, 'NFKBIA'] == 0){
      dn=dn+1
    }else if(tmp_data[i, 'EGR1'] == 0 && tmp_data[i, 'NFKBIA'] > 0){
      sp_cd48=sp_cd48 + 1
    }
  }
  row_name = c(row_name, s)
  result = rbind(result, c(dp, dn, sp_il17r, sp_cd48))
}

result = as.data.frame(result)
rownames(result) = row_name
result$sum <- rowSums(result[,c(1:4)])
colnames(result)=c("dp","EGR1","dn","NFKBIA","sum")

result$prop <-result$dp / result$sum
data<-IL7R_exp %>% group_by(sample, response) %>% dplyr::summarise(n=n()) %>% as.data.frame()
rownames(data)=data[,1]
result$response<-data[rownames(result),"response"]
colnames(result)=c("dp","EGR1","dn","NFKBIA","sum","prop","response")
compare=c("R","NR")

result[which(result$response=='Responder'),'response'] <- 'R'
result[which(result$response=='Non_responder'),'response'] <- 'NR'


pdf(paste0(output_path, '/FigureS5J.pdf'),width = 2,height = 3)
ggplot(data = result,aes(x = response, y = prop,color = response)) +
  stat_boxplot(geom = 'errorbar',width = 0.5 ,color = '#4d4d4d') +
  geom_boxplot(width=0.7, size=0.6,outlier.shape = NA,color = '#4d4d4d')+
  geom_jitter(width = 0.3,size = 2) +
  stat_compare_means(comparisons = list(c('R', 'NR')), label = 'p.format',face="bold")+
  scale_color_manual(values =color_palette_Response) +
  theme_classic() +
  xlab("") +
  ylab("EGR1+NFKBIA+ T cell proporation") +
  theme(legend.position = "no",
        axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        axis.title.x=element_text(colour="black"),
        axis.title.y=element_text(colour="black"))
dev.off() 

