####add normal CD4T to this work 
library(Seurat)
library(devtools)
library(clustree)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(ggplot2)
library(ggExtra)
library(DoubletFinder)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(miscTools)
library(ggthemes)
library(corrplot)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(viridis)
library(circlize)
library(plotly)
library(reshape)
library(pheatmap)
require(gbm)
library(future)
library(scater)
library(ggsci)
library(biomaRt)
library(ComplexHeatmap)
library(clusterProfiler)
set.seed(123)
options(future.globals.maxSize = 128000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
dir.create('./4.3add.normal.cd4.T')
CD4.external <- readRDS('./4.2process.external.T.sub/external.CD4.sub.rds')
CD4.thiswork <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
table(CD4.thiswork$minor_celltype)
table(CD4.external$CD4.CD8)
CD4.thiswork$source <- 'thiswork'
CD4.merge <- merge(CD4.external,CD4.thiswork)
table(CD4.merge$source)
CD4.merge <- CD4.merge %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
####cca
split.object <- SplitObject(CD4.merge,split.by = 'source')
pca_cca_pipline=function(listfiles){
  features <- SelectIntegrationFeatures(object.list = listfiles)
  
  listfiles <- lapply(X = listfiles, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  anchors <- FindIntegrationAnchors(object.list = listfiles, reduction =  "rpca")
  listfiles.integrated <- IntegrateData(anchorset = anchors)
  return(listfiles.integrated)
}
integrated.data <- pca_cca_pipline(split.object)
integrated.data <- ScaleData(integrated.data)
integrated.data <- RunPCA(integrated.data,npcs = 100)
DimPlot(CD4.merge,reduction = 'pca',split.by = 'source',group.by = 'source')
DimPlot(integrated.data,reduction = 'pca',split.by = 'source',group.by = 'source')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(integrated.data,harmony = F)
DefaultAssay(integrated.data) <- 'integrated'
pc <- 45
integrated.data <- FindNeighbors(integrated.data,dims = 1:pc,reduction = 'pca',verbose = F)
set.resolutions <- seq(0.1,1,0.1)
integrated.data <- FindClusters(integrated.data,resolution = set.resolutions,verbose = F)
integrated.data <- RunUMAP(integrated.data,reduction = 'pca',dims = 1:pc,verbose = F)
pdf('./4.3add.normal.cd4.T/cca.integrated.CD4.pdf')
clustree(integrated.data)
DimPlot(integrated.data,reduction = 'umap',split.by = 'source',ncol = 3,group.by = 'source')
DimPlot(integrated.data,group.by = 'source',reduction = 'umap')
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = integrated.data, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x))
  print(p)
})
dev.off()
# ####harmonny
# library(harmony)
# harmony_seurat <- RunHarmony(CD4.merge,group.by.vars = 'source',verbose = F)
# DimPlot(harmony_seurat,reduction = 'pca',split.by = 'source',group.by = 'source')
# DimPlot(harmony_seurat,reduction = 'harmony',split.by = 'source',group.by = 'source')
# do_Elbow_quantitative(harmony_seurat,harmony = T)
# pc <- 15
# harmony_seurat <- FindNeighbors(harmony_seurat,dims = 1:pc,reduction = 'harmony',verbose = F)
# set.resolutions <- seq(0.1,1,0.1)
# harmony_seurat <- FindClusters(harmony_seurat,resolution = set.resolutions,verbose = F)
# harmony_seurat <- RunUMAP(harmony_seurat,reduction = 'harmony',dims = 1:pc,verbose = F)
# pdf(file = "4.3add.normal.cd4.T/harmony.integrated.CD4.pdf",width = 8)
# clustree(harmony_seurat)
# DimPlot(harmony_seurat,reduction = 'umap',split.by = 'source',ncol = 3,group.by = 'source')
# sce.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = harmony_seurat, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
#   print(p)
# })
# dev.off()


DefaultAssay(integrated.data) <- 'RNA'
integrated.data <- integrated.data %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
integrated.data$seurat_clusters <- integrated.data$integrated_snn_res.0.5
Idents(integrated.data) <- integrated.data$seurat_clusters
table(Idents(integrated.data))
integrated.data$source <- factor(integrated.data$source,levels = c('GSE131907_nLN','obj_gse182436_nln','r1_nln','r2_nln','r3_nln','thiswork'),ordered = T)
pdf('./4.3add.normal.cd4.T/CD4.dimplot.pdf')
DimPlot(integrated.data,group.by = 'seurat_clusters',cols = sample(Palettes[['mycols_22']],size = 21))
DimPlot(integrated.data,cells = WhichCells(integrated.data,idents = '12'),pt.size = .5,group.by = 'source',order = c('GSE131907_nLN','obj_gse182436_nln','r1_nln','r2_nln','r3_nln','thiswork'))
dev.off()
plan("multisession",workers = 11)
cluster.pos.markers <- FindAllMarkers(integrated.data, only.pos = TRUE, group.by = "seurat_clusters",test.use = 'MAST',latent.vars = 'source',
                                      min.diff.pct = 0.1,logfc.threshold = 0.1)
plan('multisession',workers = 1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.3add.normal.cd4.T/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.3add.normal.cd4.T/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.3add.normal.cd4.T/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

library(Nebulosa)
pdf('./4.3add.normal.cd4.T/marker.texpre.pdf')
FeaturePlot(integrated.data,features = c('MIR155HG','LAG3','TSHZ2','C1orf228'),order = T)
plot_density(integrated.data,features = c('MIR155HG'),reduction = 'umap')
genelist <- c('CCR7','PDCD1','TIGIT','HAVCR2','CTLA4','MIR155HG','LAG3','TSHZ2','C1orf228')
exp.matrix <- GetAssayData(integrated.data, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, integrated.data$seurat_clusters, mean)
  return(a)
})
library(vegan)
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
Heatmap(t(cluster.score.normailzed), 
        column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

pdf('./4.3add.normal.cd4.T/marker.texpre.dotplot.pdf',width = 10,height = 5)
DotPlot(integrated.data,features = c('MIR155HG','LAG3','TSHZ2','C1orf228'),group.by = 'seurat_clusters',scale.by = 'size') + scale_color_gradientn(colors = Palettes[['blueRed']]) +coord_flip()
dev.off()

selected.cluster <- subset(integrated.data,seurat_clusters=='12')
selected.cluster$cluster <- '12'
source('./function/ratio_plot.R')
pdf('./4.3add.normal.cd4.T/texpre.ratio.pdf')
ratio.plot(seurat.object = selected.cluster, id.vars1 = "source", id.vars2 = "cluster", angle = 60)
dev.off()
table(selected.cluster$source)
sum(!is.na(match(colnames(CD4.thiswork[,CD4.thiswork@meta.data$minor_celltype=='C9 CD4Tex_pre-MIR155HG']),colnames(selected.cluster))))
Idents(integrated.data) <- integrated.data$seurat_clusters
saveRDS(integrated.data,file = './4.3add.normal.cd4.T/CD4.merge.rds')

integrated.data <- readRDS('./4.3add.normal.cd4.T/CD4.merge.rds')
integrated.data$paper <- as.character(integrated.data$source)
table(integrated.data$paper)
integrated.data$paper[grep(integrated.data$paper,pattern = '^r[0-9]')] <- 'roider_nln'
metadata <- integrated.data@meta.data

df1 <- data.frame(paper = rep('a',4),cluster12 = 0,CD4 = 0,percentage = 0)
df1$paper <- names(table(metadata$paper))
df1$CD4 <- as.numeric(table(metadata$paper))
df1$cluster12 <- as.numeric(table(metadata$paper[metadata$seurat_clusters=='12']))
df1$percentage <- round(df1$cluster12/df1$CD4,3)*100

index <- which(!is.na(metadata$group))
metadata$group[index] <- paste('thiswork_',metadata$group[index])
metadata$group[-index] <- metadata$paper[-index]
table(metadata$group)
df2 <- data.frame(paper = rep('a',2),cluster12 = 0,CD4 = 0,percentage = 0)
metadata_2 <- metadata[metadata$paper=='thiswork',]
df2$paper <- names(table(metadata_2$group))
df2$CD4 <- as.numeric(table(metadata_2$group))
df2$cluster12 <- as.numeric(table(metadata_2$group[metadata_2$seurat_clusters=='12']))
df2$percentage <- round(df2$cluster12/df2$CD4,3)*100

df <- rbind(df1,df2)
pdf('./4.3add.normal.cd4.T/paper.texpre.percentage.pdf')
ggbarplot(df,x = 'paper',y = 'percentage',fill = 'paper',palette = 'npg',label = T,title = 'MIR155HG+ Tex-pre') + 
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(y = '% of CD4 T cells')
dev.off()

saveRDS(integrated.data,file = './4.3add.normal.cd4.T/CD4.merge.rds')
