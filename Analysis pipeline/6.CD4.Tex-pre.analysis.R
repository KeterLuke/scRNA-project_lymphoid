####load packages####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
library(colourpicker)
library(cowplot)
library(ggsci)
library(vegan)
library(CellChat)
library(SCENIC)
library(Nebulosa)
library(tidyverse)
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source("function/clusterProfiler.enricher.R")
dir.create('./6.CD4.Tex-pre.analysis')
####load CD4.Tex-pre####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
table(Idents(CD4))
texpre <- subset(CD4,idents = 'C9 CD4Tex_pre-MIR155HG')
table(texpre$minor_celltype)
DimPlot(texpre,group.by = 'minor_celltype',reduction = 'umap')

dir.create('./6.CD4.Tex-pre.analysis/enrichment')
####enrichment####
plan('multisession',workers = 11)
celltype.DE <- FindAllMarkers(CD4,only.pos = F,group.by = 'minor_celltype',min.diff.pct = 0.1,logfc.threshold = 0.1
                             )
plan('multisession',workers = 1)
Idents(CD4) <- CD4$minor_celltype
idents <- levels(CD4)
saveFormat <- lapply(idents, function(x){
  index <- which(celltype.DE$cluster == x)
  DEGs <- celltype.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "6.CD4.Tex-pre.analysis/enrichment/celltype.all.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(celltype.DE, file = "6.CD4.Tex-pre.analysis/enrichment/celltype.all.DE.rds")

##-- Functional enrichment analysis
DEGs <- lapply(saveFormat, function(x){
  x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.1) %>% arrange(desc(avg_log2FC))
  return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "6.CD4.Tex-pre.analysis/enrichment/celltype.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "6.CD4.Tex-pre.analysis/enrichment/celltype.DEGs.rds")
####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
DEGs <- readRDS('6.CD4.Tex-pre.analysis/enrichment/celltype.DEGs.rds')
DEGs <- DEGs[['C9 CD4Tex_pre-MIR155HG']]
deg <- DEGs$gene
source("function/clusterProfiler.enricher.R")
####enrich REACTOME####
dir.create('6.CD4.Tex-pre.analysis/enrichment/REACTOME')
pdf('6.CD4.Tex-pre.analysis/enrichment/REACTOME/celltype.clusterProfiler.REACTOME.pdf')
{
  res <- cluterProfiler.enricher(gene = deg, geneType = "SYMBOL", db.type = "MsigDB",MsigDB.category = list(C2 = 'REACTOME'),saveDir = paste0(getwd(),'/6.CD4.Tex-pre.analysis/enrichment/REACTOME/'),
                                 filename = 'REACTOME.xlsx',title = 'CD4 Tex-pre', qvalueCutoff = 0.2, pvalueCutoff = 0.05,minGSSize = 10)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio

}
dev.off()
write.xlsx(res, file = "6.CD4.Tex-pre.analysis/enrichment/REACTOME/clusterProfiler.reactome.result.xlsx", sheetName = 'Tex-pre', rowNames = F)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/REACTOME/clusterProfiler.reactome.result.rds')
####show results
# res <- readRDS('./6.CD4.Tex-pre.analysis/enrichment/REACTOME/clusterProfiler.reactome.result.rds')
# res <- res %>% filter(p.adjust<0.001)
# pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
# gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
# res$gene.ratio <- gene.ratio
# write.csv(res,file = './6.CD4.Tex-pre.analysis/enrichment/REACTOME/CD4.Tex-pre.analysis_filtered.reactome.csv',quote = F,col.names = T)
res <- openxlsx::read.xlsx('./6.CD4.Tex-pre.analysis/enrichment/REACTOME/clusterProfiler.reactome.result.selected.xlsx')
enrich.res <- res
enrich.res$Description <- gsub(enrich.res$Description,pattern = '_',replacement = " ")
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

library(ggthemes)

enrich.res$wrap <- wrapText(enrich.res$Description, 30)

#barplot
enrich.res <- enrich.res[order(enrich.res$Count,decreasing = T),]
enrich.res$wrap <- factor(enrich.res$wrap,ordered = T,levels = enrich.res$wrap)
pdf('6.CD4.Tex-pre.analysis/enrichment/REACTOME/REACTOME.pdf')
p1 <- ggplot(enrich.res, 
             aes(x = Count, 
                 y = wrap, 
                 color = p.adjust,
                 fill = p.adjust))
p2 <- p1 + guides(color='none') + geom_bar(stat = 'identity') + theme_classic() + scale_color_gradient(low = "#FF9375", high = "red") +
  scale_fill_gradient(high = "yellow", low = "red")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 8)) + xlab("Count") + ylab("") + labs(fill='p.adjust',title = 'CD4 Tex-pre REACTOME enrichment') + scale_x_continuous(breaks = seq(0,24,2))
print(p2)

dev.off()
####enrich HALLMARK####
dir.create('6.CD4.Tex-pre.analysis/enrichment/HALLMARK')
pdf('6.CD4.Tex-pre.analysis/enrichment/HALLMARK/celltype.clusterProfiler.HALLMARK.pdf')
{
  res <- cluterProfiler.enricher(gene = deg, geneType = "SYMBOL", db.type = "MsigDB",MsigDB.category = list(H = c('All')),saveDir = paste0(getwd(),'/6.CD4.Tex-pre.analysis/enrichment/HALLMARK/'),
                                 filename = 'HALLMARK.xlsx',title = 'CD4 Tex-pre', qvalueCutoff = 0.2, pvalueCutoff = 0.1,minGSSize = 5)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio
  
}
dev.off()
write.xlsx(res, file = "6.CD4.Tex-pre.analysis/enrichment/HALLMARK/clusterProfiler.hallmark.result.xlsx", sheetName = 'Tex-pre', rowNames = F)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/HALLMARK/clusterProfiler.hallmark.result.rds')
####show results
res <- readRDS('./6.CD4.Tex-pre.analysis/enrichment/HALLMARK/clusterProfiler.hallmark.result.rds')
enrich.res <- res
enrich.res$Description <- gsub(enrich.res$Description,pattern = '_',replacement = " ")
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

library(ggthemes)

enrich.res$wrap <- wrapText(enrich.res$Description, 30)

#barplot
enrich.res <- enrich.res[order(enrich.res$Count,decreasing = T),]
enrich.res$wrap <- factor(enrich.res$wrap,ordered = T,levels = enrich.res$wrap)
pdf('6.CD4.Tex-pre.analysis/enrichment/HALLMARK/HALLMARK.pdf')
p1 <- ggplot(enrich.res, 
             aes(x = Count, 
                 y = wrap, 
                 color = p.adjust,
                 fill = p.adjust))
p2 <- p1 + guides(color='none') + geom_bar(stat = 'identity') + theme_classic() + scale_color_gradient(low = "#FF9375", high = "red") +
  scale_fill_gradient(high = "yellow", low = "red")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 8)) + xlab("Count") + ylab("") + labs(fill='p.adjust',title = 'CD4 Tex-pre HALLMARK enrichment') + scale_x_continuous(breaks = seq(0,50,10))
print(p2)

dev.off()
####enrich KEGG####
dir.create('6.CD4.Tex-pre.analysis/enrichment/KEGG')
pdf('6.CD4.Tex-pre.analysis/enrichment/KEGG/celltype.clusterProfiler.KEGG.pdf')
{
  res <- cluterProfiler.enricher(gene = deg, geneType = "SYMBOL", db.type = "MsigDB",MsigDB.category = list(C2= c('KEGG')),saveDir = paste0(getwd(),'/6.CD4.Tex-pre.analysis/enrichment/KEGG/'),
                                 filename = 'KEGG.xlsx',title = 'CD4 Tex-pre', qvalueCutoff = 0.2, pvalueCutoff = 0.2,minGSSize = 5)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio
  
}
dev.off()
write.xlsx(res, file = "6.CD4.Tex-pre.analysis/enrichment/KEGG/clusterProfiler.KEGG.result.xlsx", sheetName = 'Tex-pre', rowNames = F)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/KEGG/clusterProfiler.KEGG.result.rds')
####show results
res <- openxlsx::read.xlsx('./6.CD4.Tex-pre.analysis/enrichment/KEGG/clusterProfiler.KEGG.result.xlsx')
enrich.res <- res
enrich.res$Description <- gsub(enrich.res$Description,pattern = '_',replacement = " ")
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

library(ggthemes)

enrich.res$wrap <- wrapText(enrich.res$Description, 30)

#barplot
enrich.res <- enrich.res[order(enrich.res$Count,decreasing = T),]
enrich.res$wrap <- factor(enrich.res$wrap,ordered = T,levels = enrich.res$wrap)
pdf('6.CD4.Tex-pre.analysis/enrichment/KEGG/KEGG.pdf')
p1 <- ggplot(enrich.res, 
             aes(x = Count, 
                 y = wrap, 
                 color = p.adjust,
                 fill = p.adjust))
p2 <- p1 + guides(color='none') + geom_bar(stat = 'identity') + theme_classic() + scale_color_gradient(low = "#FF9375", high = "red") +
  scale_fill_gradient(high = "yellow", low = "red")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 8)) + xlab("Count") + ylab("") + labs(fill='p.adjust',title = 'CD4 Tex-pre KEGG enrichment') + scale_x_continuous(breaks = seq(0,22,10))
print(p2)

dev.off()
####enrichgo####
dir.create('6.CD4.Tex-pre.analysis/enrichment/GO')
pdf('6.CD4.Tex-pre.analysis/enrichment/GO/celltype.clusterProfiler.GO.pdf')
{
  res <- cluterProfiler.enricher(gene = deg, geneType = "SYMBOL", db.type = "GO",MsigDB.category = list(C5= c('BP')),saveDir = paste0(getwd(),'/6.CD4.Tex-pre.analysis/enrichment/GO/'),
                                 filename = 'GO.xlsx',title = 'CD4 Tex-pre', simplify = F,qvalueCutoff = 0.2, pvalueCutoff = 0.05,minGSSize = 5)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio
  
}
dev.off()
write.xlsx(res, file = "6.CD4.Tex-pre.analysis/enrichment/GO/clusterProfiler.GO.result.xlsx", sheetName = 'Tex-pre', rowNames = F)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/GO/clusterProfiler.GO.result.rds')

# deg['MIR155HG'] <- 'Homo sapiens (human) MIR155 host gene (MIR155HG)'
# GODB <- data.table::fread('./6.CD4.Tex-pre.analysis/enrichment/GOBP.DB.personal.txt',data.table = F,sep = '\t')
# 
# 
# res <- enricher(deg, TERM2GENE = GODB.ALL[,c('GO TERM','SYMBOL')],TERM2NAME = GODB.ALL[,c('GO TERM','GO NAME')], qvalueCutoff = 0.2, pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500)
# result <- res@result
# write.csv(result,file = './6.CD4.Tex-pre.analysis/enrichment/GO/unfiltered.go.result.csv',quote = F,col.names = T)
# res <- res@result %>% filter(p.adjust<0.1)
# pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
# gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
# res$gene.ratio <- gene.ratio
# write.csv(res,file = './6.CD4.Tex-pre.analysis/enrichment/GO/CD4.Tex-pre.analysis_enrichgoResult.csv',quote = F,col.names = T)
res <- openxlsx::read.xlsx('./6.CD4.Tex-pre.analysis/enrichment/GO/clusterProfiler.GO.result.selected.xlsx')
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

library(ggthemes)
enrich.res <- res

enrich.res$wrap <- wrapText(enrich.res$Description, 45)

#barplot
enrich.res <- enrich.res[order(enrich.res$Count,decreasing = T),]
enrich.res$wrap <- factor(enrich.res$wrap,ordered = T,levels = enrich.res$wrap)
pdf('6.CD4.Tex-pre.analysis/enrichment/GO/GO.pdf')
p1 <- ggplot(enrich.res, 
             aes(x = Count, 
                 y = wrap, 
                 fill = p.adjust))
p2 <- p1 + guides(color='none') + geom_bar(stat = 'identity') + theme_classic() +
  scale_fill_gradient(high = "yellow", low = "red",limits = c(0,0.05))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 8)) + xlab("Count") + ylab("") + labs(fill='p.adjust',title = 'CD4 Tex-pre GO enrichment') + scale_x_continuous(breaks = seq(0,34,2))
print(p2)

dev.off()
####GSEA####
####GSEAGO####
dir.create('./6.CD4.Tex-pre.analysis/enrichment/GSEA')
dir.create('./6.CD4.Tex-pre.analysis/enrichment/GSEA/GO')
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
gmtfile = '6.CD4.Tex-pre.analysis/enrichment/c5.go.bp.v2022.1.Hs.symbols.gmt'
geneset <- read.gmt(gmtfile)
deg <- readRDS('6.CD4.Tex-pre.analysis/enrichment/celltype.all.DE.rds')
deg <- deg[deg$cluster=='C9 CD4Tex_pre-MIR155HG',]
{
  genelist <- deg$avg_log2FC
  names(genelist) <- toupper(deg$gene)
  genelist <- sort(genelist,decreasing = T)
  res <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.9,verbose = F)
  cat(paste0('GSEA of cluster',unique(deg$cluster),'_finished','\n'))
}
head(res)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/GSEA/GO/GSEA_result.rds')
res_df <- res@result
res_df.sig <- res_df %>% dplyr::filter(abs(NES) > 1,pvalue < 0.05, p.adjust < 0.4)
write.xlsx(res_df.sig,file = '6.CD4.Tex-pre.analysis/enrichment/GSEA/GO/GSEA_sig.result.xlsx',sheetName = 'Tex-pre')

gene <- c(c('IRF4','CD7'),c('LAG3','CTLA4','CCR7','CHI3L2','HIF1A','IL7R','TIGIT'))
library(GseaVis)
id <- c('GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY','GOBP_PROGRAMMED_CELL_DEATH_INVOLVED_IN_CELL_DEVELOPMENT','GOBP_NEGATIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION',
'GOBP_IMMUNE_SYSTEM_DEVELOPMENT','GOBP_LEUKOCYTE_DIFFERENTIATION','GOBP_MONONUCLEAR_CELL_DIFFERENTIATION')
pdf('./6.CD4.Tex-pre.analysis/enrichment/GSEA/GO/GSEA.pdf',width = 12,height = 8)
gseaNb(object = res,
       geneSetID = id,
       curveCol = Palettes[['mycols_6']],
       subPlot = 3,
       addPval = T,
       rmHt = T,
       rmSegment = T,
       pvalX = 1,
       pvalY = 1,
       rank.gene = gene,
       rank.gene.nudgey = .5,termWidth = 60,arrowAngle = 20)
dev.off()
####GSEAREACTOME####
deg <- readRDS('6.CD4.Tex-pre.analysis/enrichment/celltype.all.DE.rds')
deg <- deg[deg$cluster=='C9 CD4Tex_pre-MIR155HG',]
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
dir.create('./6.CD4.Tex-pre.analysis/enrichment/GSEA/REACTOME')
gmtfile = '6.CD4.Tex-pre.analysis/enrichment/c2.cp.reactome.v2022.1.Hs.symbols.gmt'
geneset <- read.gmt(gmtfile)
{
  genelist <- deg$avg_log2FC
  names(genelist) <- toupper(deg$gene)
  genelist <- sort(genelist,decreasing = T)
  res <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.9,verbose = F)
  cat(paste0('GSEA of cluster',unique(deg$cluster),'_finished','\n'))
}
head(res)
saveRDS(res,file = '6.CD4.Tex-pre.analysis/enrichment/GSEA/REACTOME/GSEA_result.rds')
res_df <- res@result
res_df.sig <- res_df %>% dplyr::filter(abs(NES) > 1,pvalue < 0.05, p.adjust < 0.4)
write.xlsx(res_df.sig,file = '6.CD4.Tex-pre.analysis/enrichment/GSEA/REACTOME/GSEA_sig.result.xlsx',sheetName = 'Tex-pre')
gene <- c(c('IRF4','CD7'),c('LAG3','CTLA4','CCR7','CHI3L2','HIF1A','IL7R','TIGIT'))
library(GseaVis)
id <- c('REACTOME_GLUCOSE_METABOLISM','REACTOME_METABOLISM_OF_CARBOHYDRATES','REACTOME_SIGNALING_BY_INTERLEUKINS',
        'REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY','REACTOME_MAPK_FAMILY_SIGNALING_CASCADES')
pdf('./6.CD4.Tex-pre.analysis/enrichment/GSEA/REACTOME/GSEA.pdf',width = 12,height = 8)
gseaNb(object = res,
       geneSetID = id,
       curveCol = Palettes[['zissou']],
       subPlot = 3,
       addPval = T,
       rmHt = T,
       rmSegment = T,
       pvalX = 1,
       pvalY = 1,
       rank.gene = gene,
       rank.gene.nudgey = .5,termWidth = 60,arrowAngle = 20)
dev.off()
gseaplot(res,id[5])

####cellchat####
dir.create('./6.CD4.Tex-pre.analysis/Cellchat')
plan('multisession',workers = 1)
####redefine celltype(tailored for cellchat analysis)####
data.merge <- readRDS('./5.merge.minor.major/data.merge.pro.rds')
data.merge <- data.merge[,data.merge$large_annotation=="Immune"]
table(data.merge$large_annotation)
table(data.merge$minor_celltype,data.merge$group)
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
metadata <- data.merge@meta.data
anno_CD4 <- CD4@meta.data
anno_CD4$minor_celltype <- sapply(str_split(anno_CD4$minor_celltype,pattern = ' '),"[",2)
table(anno_CD4$minor_celltype,anno_CD4$group)
anno_CD4 <- anno_CD4 %>% dplyr::mutate(celltype = case_when(
  grepl('Tn',minor_celltype)~'CD4Tn',
  grepl('Tm',minor_celltype)~'CD4Tm',
  grepl('Treg',minor_celltype)~'CD4Treg',
  grepl('PDCD1',minor_celltype)~'CD4Tex',
  grepl('Tr1',minor_celltype)~'CD4Tr1',
  grepl('pre',minor_celltype)~'CD4Texpre',
  grepl('Tstd',minor_celltype)~'CD4Tstd',
  grepl('Texpre',minor_celltype)~'CD4Texpre',
  T ~ minor_celltype)) %>% dplyr::select(group,minor_celltype,celltype)
table(anno_CD4$celltype,anno_CD4$group,useNA = 'always')

anno_CD8 <- CD8@meta.data
anno_CD8$minor_celltype <- sapply(str_split(anno_CD8$minor_celltype,pattern = ' '),"[",2)
table(anno_CD8$minor_celltype,anno_CD8$group)
anno_CD8 <- anno_CD8 %>% dplyr::mutate(celltype = case_when(
   grepl('Tn',minor_celltype)~'CD8Tn',
   grepl('Tem',minor_celltype)~'CD8Tem',
   grepl('Teff',minor_celltype)~'CD8Teff',
   T ~ minor_celltype
  ))%>% dplyr::select(group,minor_celltype,celltype)
table(anno_CD8$celltype,anno_CD8$group,useNA = 'always')

anno_B <- B@meta.data
table(anno_B$minor_celltype,anno_B$group)
anno_B <- anno_B %>% dplyr::mutate(celltype = ifelse(minor_celltype=='B&Tdoublet-CD3D',NA,'B'))%>% dplyr::select(group,minor_celltype,celltype)

#   dplyr::mutate(celltype = 
#   case_when(
#   grepl('Bgc',minor_celltype)~'Bgc',
#   grepl('Bmem',minor_celltype)~'Bmem',
#   grepl('Bn',minor_celltype)~'Bn',
#   grepl('Bprol',minor_celltype)~'Bpro',
#   grepl('Bstre',minor_celltype)~'Bstd',
#   grepl('Plasma',minor_celltype)~'Plasma',
#   T ~ NA_character_
# ))
table(anno_B$celltype,anno_B$group,useNA = 'always')

df <- rbind(anno_CD4,anno_CD8,anno_B)
metadata$celltype <- as.character(metadata$minor_celltype)
table(metadata$celltype,metadata$group)
index <- match(rownames(df),rownames(metadata))
metadata$celltype[index] <- df$celltype
table(metadata$celltype,metadata$group,useNA='always')
metadata$minor_celltype[index] <- df$minor_celltype
table(metadata$minor_celltype,metadata$group,useNA='always')
metadata <- metadata[!is.na(metadata$celltype),]
data.merge <- data.merge[,colnames(data.merge)%in%rownames(metadata)]
data.merge@meta.data <- metadata
table(data.merge$celltype,data.merge$group,useNA='always')


data.merge <- NormalizeData(data.merge)
saveRDS(data.merge,file = './6.CD4.Tex-pre.analysis/Cellchat/data.merge.redefined.rds')
rm(list = ls())
gc()
#########################cell_chat analysis###############################
#load data
cell_chat_obj <- readRDS('./6.CD4.Tex-pre.analysis/Cellchat/data.merge.redefined.rds')
library('CellChat')
library(patchwork)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
showDatabaseCategory(CellChatDB)

meta=cell_chat_obj@meta.data
table(meta$celltype,meta$group)
normal_obj=subset(cell_chat_obj,subset=group=="normal")
tumor_obj=subset(cell_chat_obj,subset=group=="tumor")

table(normal_obj$celltype,normal_obj$group)
table(tumor_obj$celltype,tumor_obj$group)

cellchat1 <- createCellChat(object= normal_obj, group.by = "celltype")
cellchat1@DB=CellChatDB.use
cellchat2 <- createCellChat(object= tumor_obj, group.by = "celltype")
cellchat2@DB=CellChatDB.use

do_process_cellchat=function(cellchat){
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 1) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 6)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}


cellchat1=do_process_cellchat(cellchat1)
cellchat2=do_process_cellchat(cellchat2)

table(cellchat1@idents,cellchat1@meta$group)
table(cellchat2@idents,cellchat2@meta$group)

save(cellchat1,cellchat2,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.list.RData')

pdf("./6.CD4.Tex-pre.analysis/Cellchat/overcheck.pdf",15,15)
showDatabaseCategory(CellChatDB)
groupSize <- as.numeric(table(cellchat1@idents))
{par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Normal:Number of interactions")
netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Normal:Interaction weights/strength")}

groupSize <- as.numeric(table(cellchat2@idents))
{par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Tumor:Number of interactions")
netVisual_circle(cellchat2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Tumor:Interaction weights/strength")}

object.list <- list(normal = cellchat1, tumor = cellchat2)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}

dev.off()

pdf('./6.CD4.Tex-pre.analysis/Cellchat/overcheck of each celltype.pdf',width = 15,height=32)
mat <- cellchat1@net$count
groupSize <- as.numeric(table(cellchat1@idents))
par(mfrow=c(6,3),xpd=T)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = rownames(mat)[i])
}
title(main = "Normal:Number of interactions",outer = T,line = -2)

mat <- cellchat1@net$weight
groupSize <- as.numeric(table(cellchat1@idents))
par(mfrow=c(6,3),xpd=T)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = rownames(mat)[i])
}
title(main = "Normal:Interaction weights/strength",outer = T,line = -2)

mat <- cellchat2@net$count
groupSize <- as.numeric(table(cellchat2@idents))
par(mfrow=c(6,3),xpd=T)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = rownames(mat)[i])
}
title(main = "Tumor:Number of interactions",outer = T,line = -2)

mat <- cellchat2@net$weight
groupSize <- as.numeric(table(cellchat2@idents))
par(mfrow=c(6,3),xpd=T)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = rownames(mat)[i])
}
title(main = "Tumor:Interaction weights/strength",outer = T,line = -2)
dev.off()

pdf('./6.CD4.Tex-pre.analysis/Cellchat/Texpre.income.signal.pdf')
groupSize <- as.numeric(table(cellchat1@idents))
netVisual_circle(cellchat1@net$weight, weight.scale = T,targets.use = 'CD4Texpre' ,label.edge= F, edge.width.max = 12, title.name = paste0("Interaction weights/strength - ",names(object.list[1])))
groupSize <- as.numeric(table(cellchat2@idents))
netVisual_circle(cellchat2@net$weight, weight.scale = T,targets.use = 'CD4Texpre' ,label.edge= F, edge.width.max = 12, title.name = paste0("Interaction weights/strength - ",names(object.list[2])))
dev.off()



object.list <- list(normal = cellchat1, tumor = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
df <- subsetCommunication(cellchat)
idents <- names(df)
openxlsx::write.xlsx(df,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.all.xlsx',sheetName = idents)

df1 <- subsetCommunication(cellchat,sources.use = 'CD4Texpre')
openxlsx::write.xlsx(df1,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.CD4Texpre.out.xlsx',sheetName = idents)

df2 <- subsetCommunication(cellchat,targets.use = 'CD4Texpre')
openxlsx::write.xlsx(df2,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.CD4Texpre.in.xlsx',sheetName = idents)

pdf("./6.CD4.Tex-pre.analysis/Cellchat/splited_overcheck.pdf",10,10)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

{par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")}

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

dev.off()

df1 <- subsetCommunication(cellchat1)
df2 <- subsetCommunication(cellchat2)
idents1 <- levels(cellchat1@idents)
idents2 <- levels(cellchat2@idents)
group.celltypes <- factor(c('B','CD4T','CD4T','CD4T','CD4T','CD4T','CD4T','CD4T','CD4T','CD8T','CD8T','CD8T','Myeloid','Myeloid','Myeloid','Myeloid','NK','Myeloid'))
names(group.celltypes) <- levels(cellchat1@idents)

pdf('./6.CD4.Tex-pre.analysis/Cellchat/all.signal.chord.pdf',8,8)
netVisual_chord_cell(cellchat1,signaling = unique(df1$pathway_name),title.name = 'normal',group = group.celltypes,big.gap = 10)
netVisual_chord_cell(cellchat2,signaling = unique(df2$pathway_name),title.name = 'tumor',group = group.celltypes,big.gap = 10)
dev.off()
#####################
library(wesanderson)
col=wes_palette("Royal1")[1:2]
pdf("./6.CD4.Tex-pre.analysis/Cellchat/cell_specific_sig_overcheck.pdf",10,10)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4Texpre")
gg1
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4Tex")
gg1
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC2")
gg1
dev.off()
#Compare the major sources and targets in 2D space
#Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.
#netAnalysis_computeCentrality(cellchat2)
pdf("./6.CD4.Tex-pre.analysis/Cellchat/Compare the major sources and targets.pdf",10,10)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  object.list[[i]]=netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()
##Part II: Identify the conserved and context-specific signaling pathways
#function 
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
saveRDS(cellchat,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.celltype.rds')
##the following steps are processed via jupyternote book####
#> Compute signaling network similarity for datasets 1 2
#cellchat <- netEmbedding(cellchat, type = "functional") ####process via jupyter notebook
#> Manifold learning of the signaling networks for datasets 1 2
#cellchat <- netClustering(cellchat, type = "functional")####process via jupyter notebook
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,title = "functional_cluster")
#> 2D visualization of signaling networks from datasets 1 2
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
#cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
#cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,title = "structural_cluster")
######
cellchat <- readRDS('./6.CD4.Tex-pre.analysis/Cellchat/cellchat.celltype.rds')
#> 2D visualization of signaling networks from datasets 1 2
pdf("./6.CD4.Tex-pre.analysis/Cellchat/signaling groups based on their functional_structure similarity.pdf",10,10)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,title = "functional_cluster")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,title = "structural_cluster")
dev.off()



pdf("./6.CD4.Tex-pre.analysis/Cellchat/signaling groups based on their functional_structure similarity_distance.pdf",10,10)
library(viridis)
col=viridis(n=20)
rankSimilarity(cellchat, type = "functional",color.use =col )
rankSimilarity(cellchat, type = "structural")
dev.off()

#Compare the overall information flow of each signaling pathway
pdf("./6.CD4.Tex-pre.analysis/Cellchat/Compare the overall information flow of each signaling pathway.pdf",10,10)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pdf("./6.CD4.Tex-pre.analysis/Cellchat/combining all the identified signaling pathways.pdf",15,30)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.8, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1],width = 15, height = 18, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.8, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 18, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.8, "cm"))
dev.off()
##Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
cellchat
pdf("./6.CD4.Tex-pre.analysis/Cellchat/Identify the upgulated and down-regulated signaling ligand-receptor pairs bubble.pdf",10,10)
netVisual_bubble(cellchat, sources.use = 'CD4Texpre', comparison = c(1, 2), 
                 angle.x = 45, remove.isolate = F,thresh = 0.05)
gg1 <- netVisual_bubble(cellchat, sources.use = 'CD4Texpre',comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in tumor", angle.x = 45, remove.isolate = F,thresh = 0.05)
gg1
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 'CD4Texpre',  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in tumor", angle.x = 45, remove.isolate = F,thresh = 0.05)
gg2
#> Comparing communications on a merged object


netVisual_bubble(cellchat, targets.use = 'CD4Texpre', comparison = c(1, 2), 
                 angle.x = 45, remove.isolate = F,thresh = 0.05)
gg1 <- netVisual_bubble(cellchat, targets.use = 'CD4Texpre',comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in tumor", angle.x = 45, remove.isolate = F,thresh = 0.05)
gg1
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use = 'CD4Texpre',  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in tumor", angle.x = 45, remove.isolate = F,thresh = 0.05)
gg2
#> Comparing communications on a merged object

dev.off()

saveRDS(cellchat,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.celltype.rds')
saveRDS(object.list,file = './6.CD4.Tex-pre.analysis/Cellchat/cellchat.list.rds')

##########
cellchat <- readRDS('./6.CD4.Tex-pre.analysis/Cellchat/cellchat.celltype.rds')
object.list <- readRDS('./6.CD4.Tex-pre.analysis/Cellchat/cellchat.list.rds')
load('./6.CD4.Tex-pre.analysis/Cellchat/cellchat.list.RData')

pdf("./6.CD4.Tex-pre.analysis/Cellchat/Identify_specific_path.pdf",10,10)
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("normal","tumor")) # set factor level
plotGeneExpression(cellchat, signaling = "MIF", split.by = "datasets", colors.ggplot = T)+scale_color_aaas()
plotGeneExpression(cellchat, signaling = "MHC-II", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "GALECTIN", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "ICOS", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "ICAM", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "SELPLG", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "IL16", split.by = "datasets", colors.ggplot = T)

dev.off()

pdf('./6.CD4.Tex-pre.analysis/Cellchat/compare selected pairways between groups in CD4Texpre.pdf')
pairLR <- CellChat::extractEnrichedLR(cellchat,signaling = c('GALECTIN','CD86'))
pairLR.use1 <- data.frame(interaction_name = pairLR[c(1,3),])
p1 <- netVisual_bubble(cellchat1,pairLR.use = pairLR.use1,sources.use = c('cDC1','cDC2','LAMP3+cDC','Macro'),targets.use = c('CD4Texpre'),title.name = 'normal')
p2 <- netVisual_bubble(cellchat2,pairLR.use = pairLR.use1,sources.use = c('cDC1','cDC2','LAMP3+cDC','Macro'),targets.use = c('CD4Texpre'),title.name = 'tumor')
p1 + p2
pairLR.use2 <- data.frame(interaction_name = pairLR[c(4,5),])
p3 <- netVisual_bubble(cellchat1,pairLR.use = pairLR.use2,sources.use = c('cDC1','cDC2','LAMP3+cDC','Macro'),targets.use = c('CD4Texpre'),title.name = 'normal')
p4 <- netVisual_bubble(cellchat2,pairLR.use = pairLR.use2,sources.use = c('cDC1','cDC2','LAMP3+cDC','Macro'),targets.use = c('CD4Texpre'),title.name = 'tumor')
p3 + p4
dev.off()

save(object.list, cellchat1, cellchat2, cellchat, file = "./6.CD4.Tex-pre.analysis/Cellchat/cellchat.rData")


