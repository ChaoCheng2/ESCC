# Figure S4---------------------------------------------------------------------
rm(list = ls())
gc()


# Figure S4A 
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

pdf(paste0(output_path,'/FigureS4A.pdf'),width = 6,height = 4.5)
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

##### FigureS4B
ccgene = c('CDC16', 'CDC7', 'CDC45', 'GADD45B', 'DBF4', 'ANAPC1', 'CREBBP', 
           'MDM2', 'ABL1', 'SMC1B', 'GADD45G', 'ATM', 'ATR', 'ANAPC7', 'RBL2',
           'ANAPC5', 'RBL1', 'MYC', 'CDC14B', 'SMC1A', 'CDC14A', 
           'SKP1', 'TGFB2', 'TGFB1', 'GADD45A', 'STAG1', 'PLK1', 'RBX1', 
           'STAG2', 'TGFB3', 'MCM4', 'ORC6', 'CCND1', 'MAD1L1', 'MCM3', 'MCM6', 
           'MCM5', 'YWHAB', 'CCNA2', 'MCM7', 'BUB1', 'CHEK1', 'WEE2', 'MCM2', 
           'PTTG1', 'CDC27', 'CDC25B', 'CDC25C', 'CDC25A', 'CDC6', 'CDC20', 
           'BUB3', 'YWHAZ', 'CCND2', 'YWHAH', 'CCNB1', 'YWHAG', 'YWHAE',
           'CCNE1', 'ORC3', 'CCND3', 'PTTG2', 'SFN', 'E2F1', 'TFDP1', 'ZBTB17',
           'CDK1', 'ESPL1', 'ANAPC10', 'RAD21', 'BUB1B', 'ANAPC11', 'RB1', 
           'SKP2', 'CUL1', 'SMAD3', 'SMAD4', 'ANAPC2', 'TFDP2', 'PRKDC', 
           'MAD2L1', 'ANAPC4', 'SMAD2', 'YWHAQ', 'CHEK2', 'CDC23', 'EP300', 
           'GSK3B', 'CDKN2A', 'CDKN1C', 'CDKN1B', 'CDKN1A', 'CDKN2D', 'CCNA1',
           'CDKN2B', 'CDKN2C', 'FZR1', 'SMC3', 'ANAPC13', 'PCNA', 'TTK', 'PKMYT1',
           'CDK2', 'CDC26', 'E2F5', 'CDK4', 'WEE1', 'E2F4', 'E2F3', 'TP53', 
           'E2F2', 'ORC1', 'ORC2', 'CCNE2', 'CDK6', 'ORC4', 'CCNB2', 'CDK7', 
           'MAD2L2', 'ORC5', 'HDAC1', 'HDAC2', 'CCNH', 'CCNB3', 'DNA2', 'POLE4', 
           'POLE3', 'PRIM1', 'PRIM2', 'POLD4', 'RFC4', 'RFC5', 'RPA1', 'POLA1', 
           'RPA3', 'POLD3', 'RNASEH2B', 'RPA2', 'PCNA', 'RPA4'
           # 'ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
           # 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
           # 'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1'#,
           #'MKI67','TYMS','STMN1','TUBB','TUBA1B'
)
#ccgene = unique(ccgene)
#ccgene = c(cc.genes$s.genes, cc.genes$g2m.genes)
ccgene_zzm = c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
               'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
               'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1')

Tcell_srt = readRDS('./Final_RDS/T_cell.rds')
tmeta = Tcell_srt@meta.data[,c('orig.ident', 'response', 'time')]
tmeta = tmeta[!duplicated(tmeta),]
rownames(tmeta) = tmeta$orig.ident

cyto_res = t(read.table('./Final_RDS/Cytoseq_output/output_only_cell_expand.Zscore', sep='\t'))
cyto_res = cyto_res[, c('TGFB1', 'PGE2', 'TRAIL')]
rownames(cyto_res) = gsub('\\.', '-',rownames(cyto_res))

# prolifreation
Tcell_srt$group = colnames(Tcell_srt)
cc_exp = AverageExpression(Tcell_srt, group.by = 'group', features = ccgene_zzm)$RNA
cc_exp = colMeans(cc_exp) %>% log1p()
cyto_res = as.data.frame(cyto_res)
cyto_res$proliferation = cc_exp[rownames(cyto_res)]
plodata = cyto_res
plodata$response = Tcell_srt@meta.data[rownames(plodata), 'response']
plodata$sample = Tcell_srt@meta.data[rownames(plodata), 'orig.ident']

sample_order = plodata %>%
  filter(sample!='EC_07_post_T') %>%
  filter(sample!='EC_16_pre_T') %>%
  group_by(sample) %>%
  dplyr::summarise(cor_res = cor(TRAIL, proliferation)) %>%
  arrange(cor_res)

plodata = plodata %>%
  filter(sample!='EC_07_post_T') %>%
  filter(sample!='EC_16_pre_T') %>%
  mutate(sample = factor(sample, levels = sample_order$sample))

table(plodata$sample)

pdf('./Supple_Figure/Supp_Figure4/FigS4B.pdf', width= 13, height= 6)
ggplot(data = plodata, aes(x = TRAIL, y = proliferation)) +
  geom_point(size=0.2, color='gray', alpha=0.5)+
  geom_smooth(data=plodata,mapping = aes(x=TRAIL, y=proliferation),
              method = "lm", se=F,  formula = y~x, size=0.5, inherit.aes = F, color='darkred')+
  #geom_bin2d(aes(fill = after_stat(count)))+
  facet_wrap(~sample, ncol = 8, scales = 'free') +
  stat_cor(method='pearson', size=3, label.sep='\n')+
  theme_classic()+
  labs(x='Suppression(TRAIL)', y='Proliferation')+
  theme(
    strip.background = element_rect(color="white", fill="white"),
    axis.text.x = element_text(color='black'),
    axis.text.y = element_text(color='black'))
dev.off()

# ------------------------------------------------------------------------------
# Figure S4C 

cyto_res = t(read.table('Final_RDS/Cytoseq_output/new_updated/output_cell_expand_noNK.Zscore', sep='\t'))
cyto_res_sample = t(read.table('Final_RDS/Cytoseq_output/new_updated/output_sample_expand_noNK.Zscore', sep='\t'))
cyto_res_TRAIL = cyto_res_sample[, c('TRAIL')]
cyto_res = cyto_res[, c('TGFB1', 'PGE2', 'TRAIL')]
rownames(cyto_res) = gsub('\\.', '-',rownames(cyto_res))
#
T_cell_no_NK<-readRDS("Final_RDS/T_cell_no_NK_updated_8.8.rds")

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
  summarise(TGFB1=cor(TGFB1, proliferation), 
            PDCD1=cor(PDCD1, proliferation), 
            PGE2=cor(PGE2, proliferation),
            TRAIL=cor(TRAIL, proliferation))
save[is.na(save)] = 0

colnames(save) <- c("sample","TGFB1","PDCD1-Tres","PGE2","TRAIL-Tres")

library(tidyr)
saveData_long <- gather(save[,c('sample','PDCD1-Tres','TRAIL-Tres')],key = 'gene',value = 'sample')

pdf('./Supple_Figure/Supp_Figure4/FigS4C.pdf', width= 4, height= 3)
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


# base::plot(density(save$PDCD1,width = 0.3),ylim=c(0,6),xlim=c(-0.7,0.3),col="blue",main = "",xlab="",lwd=4,lty=1)
# par(new=T)
# base::plot(density(save$TRAIL,width = 0.5),ylim=c(0,6),xlim=c(-0.7,0.3),col="darkred",main = "",xlab="",lwd=4,lty=1)
# legend("topleft",legend = c("PDCD1-Tres","TRAIL-Tres"),
#        lty=c(1,1),col=c("blue","darkred"),bty="n")

# ------------------------------------------------------------------------------
# Figure S4D 

output_dir = './Final_RDS/Cytoseq_output/new_updated/' # output path
# draw
cyto_res = t(read.table(paste0(output_dir, 'output_sample_expand_noNK.Zscore'), sep='\t'))
TRAIL_sample<-cyto_res[,"TRAIL"]
TRAIL_sample<-as.data.frame(TRAIL_sample)
TRAIL_exp = FetchData(T_cell_no_NK, vars="TNFSF10")
TRAIL_exp$sample = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "orig.ident"]
TRAIL_exp$time = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "time"]
TRAIL_exp$response = T_cell_no_NK@meta.data[rownames(TRAIL_exp), "response"]
TRAIL_exp_mean = TRAIL_exp %>%
  group_by(sample,time,response) %>%
  summarise_all(list(mean)) %>% as.data.frame()

TRAIL_sample_all<-cbind(TRAIL_exp_mean,TRAIL_sample)
TRAIL_sample_all<-TRAIL_sample_all[,c(1,2,3,5)]
colnames(TRAIL_sample_all)<-c("orig.ident","time","response","TRAIL_activity")
TRAIL_sample_all_no_16pre<-subset(TRAIL_sample_all,orig.ident!="EC_16_pre_T")
TRAIL_sample_all_no_16pre<-subset(TRAIL_sample_all_no_16pre,orig.ident!="EC_07_post_T")

color_palette_Response <- c('Non_responder' = "#0099CCFF",'Responder' = "#FFC20AFF")


pdf('./Supple_Figure/Supp_Figure4/FigS4D.pdf', width=2.5, height= 4)
ggplot(data = TRAIL_sample_all_no_16pre, aes(x=response, y=TRAIL_activity))+
  geom_boxplot(width = 0.7,color = '#4d4d4d',  outlier.size = 0.0,outlier.color = 'white',size = 0.6)+
  geom_jitter(aes(color = response),width = 0.15,size = 3)+
  scale_color_manual(values = color_palette_Response)+
  scale_x_discrete( labels = c("NR","R"))+
  stat_compare_means(comparisons = list(c('Responder', 'Non_responder')), label = 'p.format',face="bold") +
  xlab("") +
  ylab("TRAIL activity") +
  theme_classic() +
  theme(
    legend.position = 'no',
    axis.text.x = element_markdown(colour = 'black'),
    axis.text.y = element_markdown(colour = 'black'))
dev.off()

# ------------------------------------------------------------------------------
# Figure S4G-H
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

pdf(file = paste0(output_path,'/FigureS4G-H.pdf'), width = 7,height =4)
p1 + p2
dev.off()
