#####process external T sub, devided into CD4&CD8 T subclusters
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
T.sub <- readRDS('./4.1process.external.data/external.T.sub.rds')
table(T.sub$source)
table(T.sub$major_celltype)
T.sub <- T.sub %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
####cca
split.object <- SplitObject(T.sub,split.by = 'source')
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
DimPlot(T.sub,reduction = 'pca',split.by = 'source')
DimPlot(integrated.data,reduction = 'pca',split.by = 'source')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(integrated.data,harmony = F)
pc <- 15
integrated.data <- FindNeighbors(integrated.data,dims = 1:pc,reduction = 'pca',verbose = F)
set.resolutions <- seq(0.1,1,0.1)
integrated.data <- FindClusters(integrated.data,resolution = set.resolutions,verbose = F)
integrated.data <- RunUMAP(integrated.data,reduction = 'pca',dims = 1:pc,verbose = F)
pdf('./4.2process.external.T.sub/cca.integrated.external.T.pdf')
clustree(integrated.data)
DimPlot(integrated.data,reduction = 'umap',split.by = 'source',ncol = 3)
DimPlot(integrated.data,group.by = 'source',reduction = 'umap')
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = integrated.data, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x))
  print(p)
})
dev.off()

DefaultAssay(integrated.data) <- 'RNA'
integrated.data <- integrated.data %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
FeaturePlot(integrated.data,features = c('CD4','CD8A','CD8B'),order = T,min.cutoff = 0)
saveRDS(integrated.data,file = './4.2process.external.T.sub/cca.integrated.T.sub.rds')
####annotate
integrated.data <- readRDS('./4.2process.external.T.sub/cca.integrated.T.sub.rds')
####select: res 0.6
integrated.data$seurat_clusters <- integrated.data$integrated_snn_res.0.6
table(integrated.data$seurat_clusters)
pdf('./4.2process.external.T.sub/T.dimplot.pdf')
DimPlot(integrated.data,group.by = 'source')
DimPlot(integrated.data,split.by = 'source',group.by = 'source',ncol = 3)
DimPlot(integrated.data,group.by = 'seurat_clusters',label = T)
dev.off()

pdf('./4.2process.external.T.sub/marker.obseve.pdf')
FeaturePlot(integrated.data,features = c('CD4','CD8A','CD8B'),order = T,min.cutoff = 0)
FeaturePlot(integrated.data,features = c('MS4A1','CD19','CD79A','CD79B'),order = T,min.cutoff = 1)
VlnPlot(integrated.data,features = c('CD4','CD8A','CD8B'),group.by = 'seurat_clusters',ncol = 2)
dev.off()

genelist <- c('CD4','CD8A','CD8B')
pdf('4.2process.external.T.sub/cluster.marker.heatmap.pdf',width = 15)
GEB::DoHeatmap_my(integrated.data,features = genelist,group.by = 'seurat_clusters')
dev.off()

pdf('4.2process.external.T.sub/cluster.marker.heatmap2.pdf',width = 15)
genelist <- c('CD4','CD8A','CD8B')
Idents(integrated.data) <- integrated.data$seurat_clusters
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
        column_names_gp = gpar(fontsize =15),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

Idents(integrated.data) <- integrated.data$seurat_clusters
table(Idents(integrated.data))
plan("multisession",workers = 12)
cluster.pos.markers <- FindAllMarkers(integrated.data, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.1,test.use = 'MAST',latent.vars = 'source')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.2process.external.T.sub/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.2process.external.T.sub/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.2process.external.T.sub/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

df <- FetchData(integrated.data,c('seurat_clusters',c('CD4','CD8A','CD8B','GZMK','NKG7')),slot = 'data')
pdf('./4.2process.external.T.sub/marker.boxplot.pdf',width = 15)
plotlist <- lapply(genelist,function(i){
  plot <- ggboxplot(df,x = 'seurat_clusters',y = i,fill = 'seurat_clusters') + 
    theme_classic() + 
    labs(title = i,x = 'clusters',y = 'expression',fill = 'clusters')
  return(plot)
})
cowplot::plot_grid(plotlist = plotlist,ncol = 2)
dev.off()

integrated.data$pt.CD4 <- PercentageFeatureSet(integrated.data,features = 'CD4')
VlnPlot(integrated.data,group.by = 'seurat_clusters',features = 'pt.CD4')
####annotate cd4&cd8
major_celltype <- integrated.data@meta.data$seurat_clusters
major_celltype <- gsub("^3$", "CD8", major_celltype)
major_celltype <- gsub("^10$", "CD8", major_celltype)
major_celltype <- gsub("^11$", "CD8", major_celltype)
major_celltype <- gsub("^12$", "CD8", major_celltype)
major_celltype <- gsub("^15$", "CD8", major_celltype)
major_celltype[major_celltype != "CD8"] <- 'CD4'
table(major_celltype)
table(integrated.data$seurat_clusters)
integrated.data <- AddMetaData(integrated.data,major_celltype,col.name = 'CD4.CD8')
table(integrated.data$CD4.CD8)

pdf('./4.2process.external.T.sub/celltype.pdf')
DimPlot(integrated.data,group.by = 'CD4.CD8',label = T)
dev.off()

pdf('./4.2process.external.T.sub/celltype.marker.observe.pdf')
VlnPlot(integrated.data,group.by = 'CD4.CD8',features = c('CD4','CD8A','CD8B'),pt.size = .1,flip = T)
GEB::DoHeatmap_my(integrated.data,features = c('CD4','CD8A','CD8B','NKG7','GZMK'),group.by = 'CD4.CD8')
DoHeatmap(integrated.data,group.by = 'CD4.CD8',features = c('CD4','CD8A','CD8B','NKG7','GZMK'))
dev.off()

saveRDS(integrated.data,file = './4.2process.external.T.sub/annotated.T.sub.rds')

CD4 <- subset(integrated.data,subset = CD4.CD8 == "CD4")
CD4 <- DietSeurat(CD4,assays = 'RNA')
saveRDS(CD4,file = './4.2process.external.T.sub/external.CD4.sub.rds')

CD8 <- subset(integrated.data,subset = CD4.CD8 == "CD8")
CD8 <- DietSeurat(CD8,assays = 'RNA')

saveRDS(CD8,file = './4.2process.external.T.sub/external.CD8.sub.rds')


integrated.data <- readRDS('./4.2process.external.T.sub/annotated.T.sub.rds')
library(Nebulosa)
pdf('./4.2process.external.T.sub/texpre.marker.pdf')
plot_density(integrated.data,features = c('MIR155HG','LAG3','TSHZ2','C1orf228'))
dev.off()
