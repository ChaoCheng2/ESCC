# Figure S3------------------------------------------------------------------------------
rm(list = ls())
gc()

load('Final_RDS/supplement_Figure_S1A_All_cell.Rdata')
dir.create(path =  'Supple_Figure/Supp_Figure3')

color_cell_type = c("B cell" = "#E64B35FF",
                    "Endothelial cell"= "#4DBBD5FF",
                    "Fibroblast" = "#00A087FF",
                    "Macrophage/Monocyte" = "#3C5488FF", 
                    "Mast cell" = "#F39B7FFF",
                    "Neutrophils"= "#8491B4FF",
                    "Plasma" = "#91D1C2FF", 
                    "T cell" = "#fdbf6f",
                    "Tumor or epithelial cells" = "#7E6148FF")
color_cell_type2 = c("B cell" = "#E64B35FF",
                    "Endothelial cell"= "#4DBBD5FF",
                    "Fibroblast" = "#00A087FF",
                    "Macrophage" = "#3C5488FF", 
                    "Mast cell" = "#F39B7FFF",
                    "Neutrophils"= "#8491B4FF",
                    "Plasma" = "#91D1C2FF", 
                    "T cell" = "#fdbf6f",
                    "Tumor cell" = "#7E6148FF")

 
color_palette_clone <- c('#F7CCBF','#F19386', '#F06976','#C76578','#975C75')
table(data3$cell.type)
library(ggthemes)
table(All_cell_post_harmony$time)


# ------------------------------------------------------------------------------
# Figure S3A-1
All_cell_post_harmony2 <- subset(All_cell_post_harmony,cell.type!="Tumor or epithelial cells")
All_cell_post_harmony2@meta.data$response=factor(All_cell_post_harmony2@meta.data$response,levels = c("Responder","Non_responder"))
cell.prop <- as.data.frame(prop.table(table(All_cell_post_harmony2$cell.type, All_cell_post_harmony2$orig.ident)))
colnames(cell.prop)<-c("cell.type","orig.ident","proportion")

P1 = ggplot(cell.prop,aes(orig.ident,proportion,fill=cell.type)) +  
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values = color_cell_type) +
    theme_few() +
    ggtitle(" ") + 
    guides(fill=guide_legend(title=NULL)) +
    theme(
        legend.position = 'no',
        axis.text.x=element_text(colour="black", angle = 90,hjust =1), 
        axis.text.y=element_text(colour="black"))

# Figure S3A-2
cell.prop <- as.data.frame(prop.table(table(All_cell_post_harmony2$cell.type, All_cell_post_harmony2$group)))
colnames(cell.prop)<-c("cell.type","group","proportion")
cell.prop$group <- factor(cell.prop$group,levels = c("Non_responder_pre","Non_responder_post" ,"Responder_pre","Responder_post" ))

P2 = ggplot(cell.prop,aes(group,proportion,fill=cell.type)) + 
    geom_bar(stat="identity",position="fill") +
    scale_fill_manual(values = color_cell_type) +
    theme_few() +
    ggtitle(" ") + 
    scale_x_discrete(labels = c("NR(Pre-T)","NR(Post-T)","R(Pre-T)","R(Post-T)")) +
    guides(fill=guide_legend(title=NULL)) +
    theme(
        axis.text.x=element_text(colour="black", angle = 90,hjust =1), 
        axis.text.y=element_text(colour="black"))
 
pdf(file = 'Supple_Figure/Supp_Figure3/FigS3A.pdf',width = 10,height = 3.5 )
ggarrange(P1,P2,widths = c(5,1))
dev.off()

 
# ------------------------------------------------------------------------------
# Figure S3B 
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

# Figure S3C
data <- FetchData(All_cell_post_harmony,vars = c("cell.type","PDCD1"))
my_comparisons <- list(c('T cell', 'B cell'),c('T cell', 'Macrophage/Monocyte'))
data2 = data %>%
    group_by(cell.type) %>%
    summarise_all(list(mean)) %>% as.data.frame()
data3<-merge(data2,data,by="cell.type")
colnames(data3)<-c("cell.type","mean","PDCD1")

pdf(file = 'Supple_Figure/Supp_Figure3/FigS3C.pdf',width = 4,height = 3)
ggplot(data3,aes(x=reorder(cell.type,-PDCD1),y=PDCD1,fill=mean))+
    geom_violin(color="black",width = 0.9,scale = 'width',size = 0.01)+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.label") +
    scale_fill_gradientn( colours = color_palette_clone) +
    xlab("")+
    ylab("PDCD1 expression")+
    theme_classic() +
    theme(
        axis.text.x=element_text(colour="black", angle = 90,hjust =1,vjust = 0.5))
dev.off()

# ------------------------------------------------------------------------------
# Figure S3D

IL7R_exp = FetchData(All_cell_post_harmony, vars="PDCD1")
IL7R_exp$sample = All_cell_post_harmony@meta.data[rownames(IL7R_exp), "orig.ident"]
IL7R_exp$cell.type = All_cell_post_harmony@meta.data[rownames(IL7R_exp), "cell.type"]
#IL7R_exp$time = All_cell@meta.data[rownames(IL7R_exp), "time"]
IL7R_exp_mean = IL7R_exp %>%
    group_by(sample,cell.type) %>%
    summarise_all(list(mean)) %>% as.data.frame()
IL7R_exp_mean[IL7R_exp_mean[, 'cell.type'] == 'Tumor or epithelial cells',  'cell.type'] = 'Tumor cell'
IL7R_exp_mean[IL7R_exp_mean[, 'cell.type'] == 'Macrophage/Monocyte',  'cell.type'] = 'Macrophage'
my_comparisons <- list(c('T cell', 'B cell'),c('T cell', 'Macrophage'))
IL7R_exp_mean <- IL7R_exp_mean[-37,]

pdf('Supple_Figure/Supp_Figure3/FigS3D.pdf',width = 5,height = 4)
ggplot(IL7R_exp_mean, aes(x=reorder(cell.type,-PDCD1),y=PDCD1)) +
    geom_boxplot(position = "dodge2")+
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
    geom_jitter(shape= 21,width = 0.2,aes(fill = cell.type),size = 2)+
    scale_fill_manual(values = color_cell_type2) +
    stat_compare_means(comparisons = my_comparisons,method = "t.test",size=4,label = "p.format",angle=0,vjust = 0, hjust=0)+
    theme_classic()+
    theme(legend.position = 'no',
          axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = 0.5), 
          axis.text.y=element_text(colour="black"))
dev.off() 


# ------------------------------------------------------------------------------
# Figure S3E 
load('Final_RDS/supplement_Figure_S1A_Tcell.Rdata')

pdf('Supple_Figure/Supp_Figure2/FigS2D.pdf',width = 12,height = 4)
DotPlot(T_cell, features=c("CD3D","CD3E","CD4","CD8A","CCR7","TCF7","IL7R","FOXP3", "TNFRSF4","IL2RA","IL17A","IL22","CXCL13","TOX2","PDCD1",
                             "IFNG","ZNF683","CCL4L2","TIGIT","LAG3","ANXA1","LMNA","MT1X","IFIT2","IFI6","IFIT3","IFIT1","CX3CR1","FGFBP2","FCGR3A","MKI67","TYMS","STMN1",
                             "TYROBP","FCER1G","TRDC"), group.by = "cluster_name")+ 
    scale_color_gradientn(colours = c("#053061","#2166ac","#4393c3","#d1e5f0","#fddbc7","#d6604d","#b2182b","#67001f"))+
    theme(axis.text.x=element_text(colour="black",angle = 90,hjust = 1,vjust = 0.5))
dev.off() 


colours = c("#67001f","#b2182b","#d6604d","#fddbc7","#d1e5f0","#4393c3","#2166ac","#053061")
 

# ------------------------------------------------------------------------------ 
# Figure S3F  
library(tidyr)
cell_distribution <- function(obj, group_by1, group_by2, show.sign=FALSE, angle_col=45,order=NULL){
    
    metadata = obj@meta.data
    # Roe for each cluster
    observe_data1 = metadata %>% 
        group_by_(group_by1, group_by2) %>% 
        summarise(cell_num=n())
    colnames(observe_data1) = c('group_by1', 'group_by2', 'cell_num')
    
    observe_data1 = observe_data1 %>% tidyr::spread(key =group_by1, value = cell_num, fill=0) %>% as.data.frame()
    
    rownames(observe_data1) = observe_data1[,1]
    observe_data1 = observe_data1[, 2:ncol(observe_data1)]
    observe_data2 = rowSums(observe_data1) - observe_data1
    
    expected_data = observe_data1
    pvalue = c()
    for(i in 1:ncol(observe_data1)){
        chisq_table = matrix(c(observe_data1[,i], observe_data2[,i]), nrow = 2, byrow = TRUE)
        chisq_res = chisq.test(chisq_table)
        expected_data[, i] = chisq_res$expected[1,]
        pvalue = c(pvalue, chisq_res$p.value)
    }
    Reo = observe_data1 / expected_data
    ####
    # plot
    
    if(show.sign=='signif'){
        annote.heatmap = t(Reo)
        annote.heatmap[annote.heatmap>1] = '+++'
        annote.heatmap[annote.heatmap<=1&annote.heatmap>0.8] = '++'
        annote.heatmap[annote.heatmap<=0.8&annote.heatmap>=0.2] = '+'
        annote.heatmap[annote.heatmap<0.2&annote.heatmap>0] = '+/-'
        annote.heatmap[annote.heatmap==0] = '-'
    }else if(show.sign=='value'){
        annote.heatmap = round(t(Reo),2)
    }else{
        annote.heatmap = t(Reo)
        annote.heatmap[annote.heatmap<Inf] = '' 
    }
    bk <- c(seq(0,0.99,by=0.01),seq(1.01,2,by=0.01))
    if(is.null(order)){
        order = colnames(Reo)
    }
    annote.heatmap = annote.heatmap[order,]
    pheatmap::pheatmap(t(Reo[,order]), 
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       treeheight_row = 0,
                       color = c(colorRampPalette(c("#053061","#2166ac","#4393c3","white"))(length(bk)/2),colorRampPalette(c("white","#fddbc7","#d6604d","#b2182b","#67001f"))(length(bk)/2)),
                       breaks = bk,
                       display_numbers = annote.heatmap,
                       cellwidth = 20, cellheight = 20,
                       fontsize = 10,
                       border_color = '#ffffff',
                       angle_col=angle_col,
                       main = 'Roe',
    )
    res = list('Reo'= t(Reo), 'pvalue'=pvalue, 'annote'=annote.heatmap)
    return(res)
}
srt = readRDS('./Final_RDS/T_cell.rds')
srt$response[srt$response=='Non_responder'] = 'Non-responder'
srt$response = factor(srt$response, levels = c('Responder', 'Non-responder'))
# 1
pdf('./Supple_Figure/Supp_Figure3/FigS3F_a.pdf', width= 3, height= 4)
cell_distribution(srt, show.sign = 'signif',
                  'cluster_name', 'response',  angle_col=90,
                  order=c('CD4_C3_IL17A', 'CD8_C4_IFI6', 'NK_C2_FGFBP2', 'CD4_C1_CCR7','High proliferation T','CD8_C2_CXCL13',
                          'CD8_C5_CX3CR1', 'CD8_C3_ANXA1', 'CD4_C2_FOXP3', 'CD4_C4_CXCL13', 'CD8_C1_IFNG','NK_C1_TYROBP'))
dev.off()

# 2
srt$new_group = paste0(srt$response, '_', sapply(srt$time, function(x)strsplit(x, '_')[[1]][1]))
srt$new_group = factor(srt$new_group, levels = c('Responder_Pre', 'Non-responder_Pre','Responder_Post', 'Non-responder_Post'))
pdf('./Supple_Figure/Supp_Figure3/FigS3F_b.pdf', width= 3, height= 4)
cell_distribution(srt,
                  'cluster_name', 'new_group',  show.sign = 'signif', angle_col=90,
                  order=c('CD8_C4_IFI6', 'NK_C2_FGFBP2',  'High proliferation T', 'CD4_C1_CCR7','CD8_C2_CXCL13',
                          'CD8_C3_ANXA1', 'CD4_C2_FOXP3', 
                          'CD4_C3_IL17A', 'CD8_C1_IFNG', 'NK_C1_TYROBP',
                          'CD4_C4_CXCL13', 'CD8_C5_CX3CR1'))
dev.off()
