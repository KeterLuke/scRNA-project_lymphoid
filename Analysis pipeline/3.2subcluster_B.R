#' @description: annotate B subclusters

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
data.merge <- readRDS("./2.Cluster/data.merge.pro.rds")
sub.scRNA <- subset(data.merge, subset = major_celltype == 'B')
rm(data.merge)
DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")####remove SCT assay, scale.data and dimension reduction 
# observe the batch effect
# sub.scRNA <- SCTransform(sub.scRNA, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
# sub.scRNA <- RunPCA(sub.scRNA, npcs = 30, verbose = FALSE) %>%
# RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindClusters(resolution = seq(0.2, 1.2, by = 0.1), verbose = FALSE)

## Remove previous clustering results
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(sub.scRNA@meta.data))
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-index]
View(sub.scRNA@meta.data)

DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- NormalizeData(sub.scRNA, verbose = FALSE)
sub.scRNA <- FindVariableFeatures(sub.scRNA,nfeatures = 3000,verbose = F)
sub.scRNA <- ScaleData(sub.scRNA, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(sub.scRNA))
####correct batch effect####
dir.create('./3.2Subcluster_B')
sub.scRNA <- RunPCA(sub.scRNA,features = VariableFeatures(sub.scRNA),npcs = 100,verbose = F)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(sub.scRNA,harmony = F)
####RNA.Harmony.PC15####
pdf("3.2Subcluster_B/RNA.Harmony.PC15.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 15, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC15.rds')
####RNA.Harmony.PC20####
pdf("3.2Subcluster_B/RNA.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC20.rds')
####RNA.Harmony.PC25####
pdf("3.2Subcluster_B/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC25.rds')
####RNA.Harmony.PC30####
pdf("3.2Subcluster_B/RNA.Harmony.PC30.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC30.rds')
####RNA.Harmony.PC10####
pdf("3.2Subcluster_B/RNA.Harmony.PC10.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 10, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC10.rds')
####RNA.Harmony.PC13####
pdf("3.2Subcluster_B/RNA.Harmony.PC13.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 13, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC13.rds')
####RNA.Harmony.PC35####
pdf("3.2Subcluster_B/RNA.Harmony.PC35.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 35, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.2Subcluster_B/RNA.Harmony.Integration.PC35.rds')
####select RNA Harmony pc25 resolution:0.1####
sub.B <- readRDS('./3.2Subcluster_B/RNA.Harmony.Integration.PC25.rds')
sub.B$seurat_clusters <- sub.B$RNA_snn_res.0.1
table(sub.B$seurat_clusters)
sub.B <- RunUMAP(sub.B,dims = 1:25,reduction = 'harmony',seed.use = 40,verbose = F)
sub.B <- RunTSNE(sub.B,dims = 1:25,reduction = 'harmony',seed.use = 40,verbose = F)
pdf('3.2Subcluster_B/RNA.Harmony.Integration.PC25.seed40.pdf')
set.resolutions <- seq(0.1,1,0.1)
p <- clustree(sub.B)
print(p)
p <- DimPlot(object = sub.B, reduction = 'umap',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.B, reduction = 'umap',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
p <- DimPlot(object = sub.B, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.B, reduction = 'tsne',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')

Idents(sub.B) <- sub.B$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.B, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
table(sub.B$seurat_clusters)
length(unique(sub.B$seurat_clusters))
dir.create('./3.2Subcluster_B/Annotate')
pdf('./3.2Subcluster_B/Annotate/cluster.pdf')
DimPlot(sub.B,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']],shuffle = T,pt.size = .2)
DimPlot(sub.B,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']],shuffle = T,pt.size = .2)
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')

#### Differential expression####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
Idents(sub.B) <- sub.B$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.B, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sub.B <- ScaleData(sub.B,features = c(VariableFeatures(sub.B),unique(top5.genes$gene)))
pdf("3.2Subcluster_B/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(sub.B, features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['mycols_12']],group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

top8.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
sub.B <- ScaleData(sub.B,features = c(VariableFeatures(sub.B),unique(top8.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('./3.2Subcluster_B/Annotate/cluster.top8genes.mean.pdf',height = 12)
DoHeatmap_my(sub.B,features = unique(top8.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
#### classical marker#### 
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
pdf('./3.2Subcluster_B/Annotate/classical.marker.B.pdf')
FeaturePlot(sub.B,features  = c('CD79A','CD79B','CD19','MS4A1'),reduction = 'umap',order = T,ncol = 2)
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.Naive.pdf')
FeaturePlot(sub.B,features  = c('MS4A1','IGHD','TCL1A','IL4R'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']],min.cutoff = 0)##Naive
FeaturePlot(sub.B,features  = c('MS4A1','IGHD','TCL1A','IL4R'),reduction = 'tsne',order = T,cols = Palettes[['greyMagma']],min.cutoff = 0)##Naive
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.Memory.pdf')
FeaturePlot(sub.B,features  = c('MS4A1','CD27','AIM2','TNFRSF13B'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']])##memory
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.GC.pdf')
FeaturePlot(sub.B,features  = c('AICDA','LRMP','RGS13','SUGCT'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']],min.cutoff = 0)##GC
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
##cluster similarity#####
expMatrix <- GetAssayData(sub.B, slot = "scale.data")
highVariableGenes <- VariableFeatures(sub.B)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, sub.B$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.2Subcluster_B/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.2Subcluster_B/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
####add module score####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
Idents(sub.B) <- sub.B$seurat_clusters
signature <- list(Naive = c('YBX3','IGHD','TCL1A','FCER2'),
                  Memory = c('CD27','TNFRSF13B','AIM2'),
                  GC = c('AICDA','MKI67','RGS13','LRMP','SUGCT'))
sce <- AddModuleScore(sub.B,features = signature,name = names(signature))
pdf('./3.2Subcluster_B/Annotate/Signature_score.pdf')
VlnPlot(sce,features = c('Naive1','Memory2','GC3'),ncol = 2,pt.size = 0,cols = Palettes[['mycols_12']])
FeaturePlot(sce,features =c('Naive1','Memory2','GC3'),reduction = 'umap',ncol = 2,order = T,min.cutoff = 0,cols = Palettes[['greyMagma']],pt.size = 0.01)
dev.off()
sub.B <- sce
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
###############celltypist predict########
library(reticulate)
py_config()
scanpy <- import('scanpy')
pandas <- import('pandas')
numpy <- import('numpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(sub.B[['RNA']]@counts)))),
                       obs = pandas$DataFrame(sub.B@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(sub.B[['RNA']]@counts),
                                                         row.names = rownames(sub.B[['RNA']]@counts)))
)
adata
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- sub.B@meta.data$seurat_clusters
sub.B <- AddMetaData(sub.B,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(sub.B$celltypist_minor,sub.B$seurat_clusters)

write.xlsx(predicted_minor,file = '3.2Subcluster_B/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(sub.B,file = '3.2Subcluster_B/sub.B.rds')

####expression of functional markers####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
cell.type.markers <- read.table(file = "./3.2Subcluster_B/Annotate/B_markers.txt", header = T, stringsAsFactors = F, sep = "\t")
exp.matrix <- GetAssayData(sub.B, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.B$seurat_clusters, mean)
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
pdf("3.2Subcluster_B/Annotate/cluster.functional.signature.pdf",height = 10)
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
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
pdf('./3.2Subcluster_B/Annotate/cluster0.pdf')
plot_density(sub.B,features = c('TNFRSF13B'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster1.pdf')
plot_density(sub.B,features = c('IGHD'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster2.pdf')
plot_density(sub.B,features = c('HSPA1A'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster3.pdf')
plot_density(sub.B,features = c('MKI67'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster4.pdf')
plot_density(sub.B,features = c('CD3D'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster5.pdf')
plot_density(sub.B,features = c('LRMP'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster6.pdf')
plot_density(sub.B,features = c('XBP1'),reduction = 'umap')
dev.off()
####annotate
minor_celltype <- sub.B@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "Bmem-TNFRSF13B", minor_celltype)
minor_celltype <- gsub("^1$", "Bn-IGHD", minor_celltype)
minor_celltype <- gsub("^2$", "Bstressed-HSPA1A", minor_celltype)
minor_celltype <- gsub("^3$", "Bproliferating-MKI67", minor_celltype)
minor_celltype <- gsub("^4$", "B&Tdoublet-CD3D", minor_celltype)
minor_celltype <- gsub("^5$", "Bgc-LRMP", minor_celltype)
minor_celltype <- gsub("^6$", "Plasma-XBP1", minor_celltype)


table(minor_celltype)
sub.B <- AddMetaData(sub.B,minor_celltype,col.name = 'minor_celltype')
table(sub.B$minor_celltype)
sub.B$minor_celltype <- factor(sub.B$minor_celltype,levels = names(table(sub.B$minor_celltype)))
Idents(sub.B) <- sub.B$minor_celltype
DefaultAssay(sub.B) <- "RNA"

source('./function/do_dimplot.R')
pdf("3.2Subcluster_B/Annotate/cellType.pdf",height = 10,width = 12)
DoDimplot(sub.B,groupby = 'minor_celltype',colors = Palettes[['summerNight']])
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.B, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = sub.B, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['summerNight']])
ratio.plot(seurat.object = sub.B, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.B,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['summerNight']])
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.B,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['summerNight']],ncol = 2)
dev.off()

saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.pro.rds')

####use plot1cell to visualize####
library(plot1cell)
sub.B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
table(Idents(sub.B))
circ_data <- prepare_circlize_data(sub.B, scale = 0.8 )
cluster_colors< Palettes[['summerNight']]
group_colors<-rand_color(length(names(table(sub.B$group))))
rep_colors<-rand_color(length(names(table(sub.B$orig.ident))))
pdf('3.2Subcluster_B/Annotate/circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()
####expression of celltype markers####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
genelist <- c('CD3D','CD3E','CD7','IL7R','LRMP','RGS13','AICDA','SUGCT','CD27','TNFRSF13B','AIM2','VIM','IGHD','TCL1A','FCER2','YBX3','MKI67','TOP2A','STMN1','TUBB','HMGB2','HSPA1A',
              'HSPA1B','DNAJB1','HSPE1','HSPH1','HSPD1','XBP1','MZB1','SDC1','IGHG1','IGHG2','IGKC'
              )
sub.B@meta.data$celltype <- sub.B@meta.data$minor_celltype
Idents(sub.B) <- sub.B$minor_celltype
exp.matrix <- GetAssayData(sub.B, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.B$minor_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('./3.2Subcluster_B/Annotate/celltype_expression.pdf',width = 12,height = 22)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(24, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

Idents(sub.B) <- sub.B$minor_celltype
pdf('./3.2Subcluster_B/Annotate/celltype_expression.dotplot.pdf',width = 10)
DotPlot(sub.B,features = genelist,scale.by = 'size',col.min = 0) + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(sub.B) <- sub.B$minor_celltype
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.pro.rds')
####test of celltype percentage between groups####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
df <- as.data.frame(table(sub.B$minor_celltype,sub.B$sample))
df$total <- apply(df,1,function(x){sum(df$Freq[df$Var2==x[2]])})
df$percent <- round(df$Freq/df$total,4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('./3.2Subcluster_B/Annotate/celltype.percentage.group.pdf',height = 12,width = 15)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    theme_classic() + scale_y_continuous(breaks = scales::pretty_breaks(n=10)) + 
    labs(title = i,x = 'Group',y = '% of all T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 4)
dev.off()
