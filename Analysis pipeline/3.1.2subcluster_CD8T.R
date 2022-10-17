#' @description: annotate CD8T subclusters

####load packages####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(openxlsx)
library(Nebulosa)
library(cowplot)
set.seed(101)
library(future)
plan("multisession", workers = 11) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
source('./function/variableFeatureSelection.R')
source('./function/Integrate_multisamples.R')
####primarily annotated results####
CD8 <- readRDS("./3.1Subcluster_T/Annotate/CD8/CD8.rds")
table(Idents(CD8))
## Remove previous clustering results
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(CD8@meta.data))
CD8@meta.data <- CD8@meta.data[,-index]
View(CD8@meta.data)
# assay=RNA
DefaultAssay(CD8) <- "RNA"
CD8 <- NormalizeData(CD8, verbose = FALSE)
CD8 <- FindVariableFeatures(CD8,nfeatures = 3000,verbose = F)
CD8 <- ScaleData(CD8, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(CD8))
####correct batch effect####
####RNA.Harmony.PC15####
pdf("3.1Subcluster_T/Annotate/CD8/RNA.Harmony.PC15.Integration.pdf")
CD8.harmony <- Harmony.integration.reduceDimension(seurat.object = CD8, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 15, npcs = 100)
dev.off()
saveRDS(CD8.harmony,file = './3.1Subcluster_T/Annotate/CD8/RNA.Harmony.Integration.PC15.rds')
####RNA.Harmony.PC10####
pdf("3.1Subcluster_T/Annotate/CD8/RNA.Harmony.PC10.Integration.pdf")
CD8.harmony <- Harmony.integration.reduceDimension(seurat.object = CD8, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 10, npcs = 100)
dev.off()
saveRDS(CD8.harmony,file = './3.1Subcluster_T/Annotate/CD8/RNA.Harmony.Integration.PC10.rds')
####RNA.Harmony.PC20####
pdf("3.1Subcluster_T/Annotate/CD8/RNA.Harmony.PC20.Integration.pdf")
CD8.harmony <- Harmony.integration.reduceDimension(seurat.object = CD8, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(CD8.harmony,file = './3.1Subcluster_T/Annotate/CD8/RNA.Harmony.Integration.PC20.rds')
####select RNA Harmony PC10,res:0.5####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/RNA.Harmony.Integration.PC10.rds')
CD8$seurat_clusters <- CD8$RNA_snn_res.0.5
length(unique(CD8$seurat_clusters))
table(CD8$seurat_clusters)
dir.create('./3.1Subcluster_T/Annotate/CD8/Annotate')
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster.pdf')
DimPlot(CD8,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_6']])
DimPlot(CD8,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_6']])
dev.off()
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.rds')
#### Differential expression####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
Idents(CD8) <- CD8$seurat_clusters
table(Idents(CD8))
cluster.pos.markers <- FindAllMarkers(CD8, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,test.use = 'MAST',
                                      latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD8/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD8/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/CD8/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
CD8 <- ScaleData(CD8,features = c(VariableFeatures(CD8),unique(top5.genes$gene)))
pdf("3.1Subcluster_T/Annotate/CD8/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(CD8, features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['darjeeling']],group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

top10.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
CD8 <- ScaleData(CD8,features = c(VariableFeatures(CD8),unique(top10.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('3.1Subcluster_T/Annotate/CD8/Annotate/cluster.top10genes.mean.pdf',height = 12)
DoHeatmap_my(CD8,features = unique(top10.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.rds')
####add module score####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
Idents(CD8) <- CD8$seurat_clusters
signature <- list(Naive = c('CCR7','LEF1','CD55','TCF7'),
                  Effector = c('IFNG','GNLY','GZMA','GZMB','GZMH','NKG7','CCL4','CCL5','CST7'),
                  Exhaustion = c('PDCD1','TIGIT','LAG3','HAVCR2','TOX2','MAF','ITM2A','LGALS3'))
sce <- AddModuleScore(CD8,features = signature,name = names(signature))
pdf('3.1Subcluster_T/Annotate/CD8/Annotate/Signature_score.pdf',width = 8,height = 12)
VlnPlot(sce,features = c('Naive1','Effector2','Exhaustion3'),ncol = 2,pt.size = 0,cols = Palettes[['darjeeling']])
FeaturePlot(sce,features = c('Naive1','Effector2','Exhaustion3'),reduction = 'umap',ncol = 2,order = T,min.cutoff = 0,cols = Palettes[['greyMagma']])
dev.off()
CD8 <- sce
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.rds')
####celltypist predict####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
library(reticulate)
py_config()
numpy <- import('numpy')
pandas <- import('pandas')
scanpy <- import('scanpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(CD8[['RNA']]@counts)))),
                       obs = pandas$DataFrame(CD8@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(CD8[['RNA']]@counts),
                                                         row.names = rownames(CD8[['RNA']]@counts)))
)
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- CD8@meta.data$seurat_clusters
CD8 <- AddMetaData(CD8,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(CD8$celltypist_minor,CD8$seurat_clusters)

write.xlsx(predicted_minor,file = '3.1Subcluster_T/Annotate/CD8/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(CD8,file = '3.1Subcluster_T/Annotate/CD8/CD8.rds')
####enrichment####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
cluster.DE <- FindAllMarkers(CD8,only.pos = F,group.by = 'seurat_clusters',min.diff.pct = 0.1,test.use = 'MAST',
                             latent.vars = 'orig.ident')
idents <- levels(CD8)
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "3.1Subcluster_T/Annotate/CD8/Annotate/cluster.all.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "3.1Subcluster_T/Annotate/CD8/Annotate/cluster.all.DE.rds")

##-- Functional enrichment analysis
DEGs <- lapply(saveFormat, function(x){
  x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.25) %>% arrange(desc(avg_log2FC))
  return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "3.1Subcluster_T/Annotate/CD8/Annotate/clusters.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "3.1Subcluster_T/Annotate/CD8/Annotate/clusters.DEGs.rds")
####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
DEGs <- readRDS('3.1Subcluster_T/Annotate/CD8/Annotate/clusters.DEGs.rds')
idents <- levels(CD8)
DEGs <- lapply(DEGs, function(x){
  return(x$gene)
})
source("function/clusterProfiler.enricher.R")
####enricher
dir.create('3.1Subcluster_T/Annotate/CD8/Annotate/enrichment')
dir.create('3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/msigDB')
pdf('3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/msigDB/cluster.clusterProfiler.enricher.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB",saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/msigDB'),
                                 filename = 'enricherResult.xlsx',title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio
  return(res)
})
dev.off()
names(res) <- names(DEGs)
write.xlsx(res, file = "3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.rds')
####enrichgo
dir.create('3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/GO')
pdf('3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/GO/cluster.clusterProfiler.enrichgo.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "GO",GO.ont = 'BP',filename = 'enrichgoResult.csv',
                                 saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/GO/'),
                                 title = x, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
  # gene ratio
  res <- res$em.res.genesymbol@result %>% filter(p.adjust<0.05) #fdr adjust
  pathway.geneNum <- unlist(strsplit(res$BgRatio, "/"))[seq(1, nrow(res),by=2)]
  gene.ratio <- as.numeric(res$Count)/as.numeric(pathway.geneNum)
  res$gene.ratio <- gene.ratio
  return(res)
})
dev.off()
names(res) <- names(DEGs)
write.xlsx(res, file = "3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/CD8/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.rds')
##cluster similarity#####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
expMatrix <- GetAssayData(CD8, slot = "scale.data")
highVariableGenes <- VariableFeatures(CD8)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, CD8$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.1Subcluster_T/Annotate/CD8/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
##expression of functional markers####
cell.type.markers <- read.table(file = "3.1Subcluster_T/Annotate/CD8/Annotate/CD8_markers.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(CD8, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, CD8$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) ##perform 0-1 standardization for all clusters per gene
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$Type, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)##perform 0-1 standardization for all clusters per celltype marker
annotation.colors <- Palettes$stallion2[1:length(unique(cell.type.markers$Type))]
names(annotation.colors) <- unique(cell.type.markers$Type)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Type, 
                                               levels = unique(cell.type.markers$Type)),
                                 col = list(Type = annotation.colors),show_annotation_name = F)
pdf("3.1Subcluster_T/Annotate/CD8/Annotate/cluster.functional.signature.pdf",height = 10)
col_fun1 <- colorRamp2(c(0, 1), c("grey", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$Type, levels = unique(cell.type.markers$Type))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, left_annotation = row.annotations,
        width = unit(10, "cm"), height = unit(25, "cm"), cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1,
        width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , 
        cluster_rows = T, show_column_names = T, show_row_names = T, 
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6),rect_gp = gpar(col = 'white'))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), 
        height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 6), rect_gp = gpar(col = 'white'),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
dev.off()
####annotate celltype####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.rds')
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster0.pdf')
plot_density(CD8,features = c('GZMK'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster1.pdf')
plot_density(CD8,features = c('CCR7'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster2.pdf')
plot_density(CD8,features = c('ANXA1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster3.pdf')
plot_density(CD8,features = c('HSPA1B'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cluster4.pdf')
plot_density(CD8,features = c('HAVCR2'),reduction = 'umap')
dev.off()
####annotate
####annotate
minor_celltype <- CD8@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "C0 CD8Teff-GZMK", minor_celltype)
minor_celltype <- gsub("^1$", "C1 CD8Tn-CCR7", minor_celltype)
minor_celltype <- gsub("^2$", "C2 CD8Tm-ANXA1", minor_celltype)
minor_celltype <- gsub("^3$", "C3 CD8Tefs-HSPA1B", minor_celltype)
minor_celltype <- gsub("^4$", "C4 CD8Tex-HAVCR2", minor_celltype)

table(minor_celltype)
CD8 <- AddMetaData(CD8,minor_celltype,col.name = 'minor_celltype')
table(CD8$minor_celltype)
CD8$minor_celltype <- factor(CD8$minor_celltype,levels = c('C1 CD8Tn-CCR7','C2 CD8Tm-ANXA1','C0 CD8Teff-GZMK','C3 CD8Tefs-HSPA1B','C4 CD8Tex-HAVCR2'))
Idents(CD8) <- CD8$minor_celltype
DefaultAssay(CD8) <- "RNA"

source('./function/do_dimplot.R') 
pdf("3.1Subcluster_T/Annotate/CD8/Annotate/cellType.pdf",height = 10,width = 12)
DoDimplot(CD8,groupby = 'minor_celltype',colors = Palettes[['mycols_6']],pt.size = 1)
dev.off()

pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/cellType.clean.pdf',width = 10,height = 8)
DimPlot(CD8,group.by = 'minor_celltype',reduction = 'umap',cols = Palettes[['mycols_6']],pt.size = 1)
DimPlot(CD8,group.by = 'minor_celltype',reduction = 'tsne',cols = Palettes[['mycols_6']],pt.size = 1)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = CD8, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = CD8, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_6']])
ratio.plot(seurat.object = CD8, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(CD8,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_6']])
dev.off()

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(CD8,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_6']],ncol = 2)
dev.off()

saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
####expression of celltype markers####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
genelist <- c('CD55','CCR7','LEF1','SELL','ANXA1','IL7R','S100A4','GZMK','GZMH','GZMA','NKG7','CCL4','CCL5',
              'HSPA1B','HSPA1A','HAVCR2','LGALS3','LYST','PMAIP1')
Idents(CD8) <- CD8$minor_celltype
exp.matrix <- GetAssayData(CD8, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, CD8$minor_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform standardization for all clusters per gene
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/celltype_expression.pdf',width = 10,height = 18)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(30, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Expression", at = c(0,1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/celltype_expression.dotplot.pdf',width = 12,height = 8)
Idents(CD8) <- CD8$minor_celltype
DotPlot(CD8,features = genelist,cols = c('grey','red'),scale.by = 'size',col.min = 0,scale.min = 0) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(CD8) <- CD8$minor_celltype
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
########test of celltype percentage between groups####
CD8 <- readRDS('3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
df <- as.data.frame(table(CD8$minor_celltype,CD8$sample))
df$percent <- round(df$Freq/ncol(CD8),4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype.percentage.group.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    stat_compare_means(comparisons = comparison,method = 'wilcox.test',paired = T) +
    theme_classic() +
    labs(title = i,x = 'Group',y = '% of CD8+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype.percentage.group.unpaired.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggplot(tmp,aes(x = group,y = percent,fill = group)) + geom_boxplot(width = .4) + geom_point() + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of CD8+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD8/Annotate/celltype.percentage.group.barplot.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggbarplot(tmp,x = 'group',y = 'percent',fill = 'group',add = 'mean_se') + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of CD8+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

saveRDS(CD8,file = '3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
####observe immune-exhaustion score in Teff and Tex clusters####
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
Idents(CD8) <- CD8$minor_celltype
df <- FetchData(CD8,vars = c('Exhaustion3','minor_celltype'),cells = WhichCells(CD8,idents = c('C0 CD8Teff-GZMK','C3 CD8Tefs-HSPA1B','C4 CD8Tex-HAVCR2')))
table(df$minor_celltype)
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/Immune.exhaustion.score.compare.pdf')
ggviolin(df,x = 'minor_celltype',y = 'Exhaustion3',fill = 'minor_celltype',palette = 'npg',add = 'boxplot') + 
  stat_compare_means(comparisons = list(c('C0 CD8Teff-GZMK','C3 CD8Tefs-HSPA1B'),
                                        c('C3 CD8Tefs-HSPA1B','C4 CD8Tex-HAVCR2'),
                                        c('C0 CD8Teff-GZMK','C4 CD8Tex-HAVCR2')),
                     method = 'wilcox.test',label = 'p.signif') + ylab('Immune exhaustion score') + 
  stat_summary(fun.data = function(x) data.frame(y=-.8, label = paste("Mean=",round(mean(x),3))), geom="text")
dev.off()

df2 <- FetchData(CD8,vars = c('Exhaustion3','minor_celltype','group'),cells = WhichCells(CD8,idents = c('C0 CD8Teff-GZMK','C3 CD8Tefs-HSPA1B','C4 CD8Tex-HAVCR2')))
pdf('./3.1Subcluster_T/Annotate/CD8/Annotate/Immune.exhaustion.score.group.compare.pdf',width = 10)
ggviolin(df2,x = 'group',y = 'Exhaustion3',fill = 'minor_celltype',palette = 'npg',add = 'boxplot',facet.by = 'minor_celltype') + 
  stat_compare_means(comparisons = list(c('normal','tumor')),
                     method = 'wilcox.test',label = 'p.signif') + ylab('Immune exhaustion score') + 
  stat_summary(fun.data = function(x) data.frame(y=-.8, label = paste("Mean=",round(mean(x),3))), geom="text")
dev.off()
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
