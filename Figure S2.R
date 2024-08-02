# Figure S2
# CNV
library(infercnv)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(ggsci)
library(ggplot2)

#########
# cnv funtion
col_map_func <- function(dat){
  res = list()
  for(n in colnames(dat)){
    tmp_name = sort(unique(dat[, n]))
    tmp_col = pal_igv()(length(tmp_name))
    names(tmp_col) = tmp_name
    res[[n]] = tmp_col
  }
  return(res)
}

seurat = readRDS('Final_RDS/updata_tumors.rds')
DimPlot(seurat, group.by = 'orig.ident') + ggsci::scale_color_igv()

cnv_data_list = readRDS('Final_RDS/all.cnv.rds')

srt_meta = seurat@meta.data[, c('response', 'time','patient')]
srt_meta = srt_meta %>% arrange(response, time,patient)

cnv_data = t(cnv_data_list$data)
bk = cnv_data_list$bk
col = cnv_data_list$col

######
genefile_dir = 'Final_RDS/hg38_filtered_order_gene.txt'
geneFile <- read.table(genefile_dir)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(colnames(cnv_data),geneFile$V1),]
expr=cnv_data[,intersect(colnames(cnv_data),geneFile$V1)]
meta_data = seurat@meta.data
meta_data$response[meta_data$response=='Non_responder'] = 'Non-responder'
meta_data$response = factor(meta_data$response, levels=c('Responder', 'Non-responder'))
meta_data$CB = rownames(meta_data)

#
expr = expr[rownames(meta_data),]

meta_data = meta_data[order(meta_data$response, meta_data$time, meta_data$patient, expr[, 'TNFSF10']),]
cell_order = rownames(meta_data)
expr = expr[cell_order,]

#注释
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), 
                                               labels = 1:22,
                                               labels_gp = gpar(cex = 1.5)))
patient_color = pal_igv()(length(unique(meta_data$patient)))
names(patient_color) = unique(meta_data$patient)
color_palette_Response <- c('Non-responder' = "#0099CCFF",'Responder' = "#FFC20AFF")

left_anno <- rowAnnotation(df = meta_data[,c('response','time', 'patient')],
                           col=list('response'=color_palette_Response,
                                    'time'=c('Pre_T'='#f768a1', 'Post_T'='#91d1c2ff'),
                                    'patient'=patient_color))



pdf('Supple_Figure/Supp_FigureS2/CNV_legend.pdf',width = 2,height = 20)
draw(left_anno)
dev.off()

scale_color_manual(values = c('#91d1c2ff','#f768a1'))  

png("FigS2_scCNV_with_reference.png",width = 15,height = 8)
normal_cnv = wx_subset_CNV(cnv_data_list, class = 'normal')
normal_expr = normal_cnv$data
normal_expr = t(normal_expr[colnames(expr), ])

ht_normal = Heatmap(normal_expr, 
                    col=circlize::colorRamp2(bk, col),
                    cluster_rows = F,cluster_columns = F, use_raster = F,
                    show_column_names = F,show_row_names = F,
                    column_split = factor(normal_cnv$gene_ann$V2, paste("chr",1:22,sep = "")),
                    column_gap = unit(2, "mm"),
                    show_heatmap_legend=FALSE,
                    row_split = rep('Reference', nrow(normal_expr)),
                    column_title = NULL)

ht = Heatmap(expr[cell_order,], use_raster = F,
             col=circlize::colorRamp2(bk, col),
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(
               title = "CNV"),
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
# draw(ht, heatmap_legend_side = "right")
draw(ht_normal%v%ht, padding = unit(c(20, 10, 10, 10), "mm"), heatmap_legend_side = "right")
dev.off()
