#' @description: annotate the cell type

####load packages####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
library(colourpicker)
library(cowplot)
plan("multisession", workers = 2) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')

#### Harmony corrected result####
data.merge <- readRDS('./2.Cluster/data.merge.rds')
DefaultAssay(data.merge) <- "RNA"

####res = 0.2####
pdf("2.Cluster/res_0.2/cluster.pdf")
data.merge@meta.data$seurat_clusters <- data.merge@meta.data$RNA_snn_res.0.2
length(unique(data.merge$seurat_clusters))
mycols_15 <- c("#00CDCD", "#A4D3EE", "#FFA500", "#FFB6C1", "#607B8B", "#8F8F8F", "#4F94CD",
               "#228B22", "#CD5C5C", "#EEDFCC", "#36648B", "#CD950C", "#009ACD", "#AB82FF", "#8B6969")
mycols_19 <- c("#66CDAA", "#458B74", "#7FFFD4", "#8EE5EE", "#7AC5CD", "#104E8B", "#009ACD", "#EEA2AD", 
               "#CD8C95", "#EEDD82", "#CDBE70", "#FF7F00", "#CD6600", "#CD2626", "#8B7D6B", "#9F79EE", "#ADADAD", "#FFD700", "#9BCD9B")
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters",cols = mycols_15)+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "seurat_clusters",cols = mycols_15)+NoLegend()
DimPlot(object = data.merge,reduction = 'umap',label = TRUE,group.by = "seurat_clusters",cols = mycols_15,split.by = 'orig.ident',ncol = 3) + NoLegend()
dev.off()

## Plot the ratio of each cell type and sample situation##
pdf("2.Cluster/res_0.2/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()
pdf("2.Cluster/res_0.2/sample.cell.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "seurat_clusters", id.vars2 = "orig.ident",color.len = mycols_19)
dev.off()


expMatrix <- GetAssayData(data.merge, slot = "scale.data")
highVariableGenes <- VariableFeatures(data.merge)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, data.merge$seurat_clusters, mean)
  return(mean.value)
})

corrMatrix <- (1- cor(t(meanExpCluster), method="pearson"))/2
library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("2.Cluster/res_0.2/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')
#### Differential expression##
data.merge <- readRDS('./2.Cluster/data.merge.rds')
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25,latent.vars = "orig.ident",test.use = 'MAST')
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "2.Cluster/Annotate/res_0.2/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "2.Cluster/Annotate/res_0.2/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "2.Cluster/Annotate/res_0.2/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
data.merge <- ScaleData(data.merge,features = c(VariableFeatures(data.merge),unique(top5.genes$gene)))
pdf("2.Cluster/Annotate/res_0.2/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(data.merge, features = unique(top5.genes$gene), size = 2,group.colors = mycols_15,group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

####The expression of the classic marker##
cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) ##perform 0-1 standardization for all clusters per gene
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$Celltype, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)##perform 0-1 standardization for all clusters per celltype marker
annotation.colors <- Palettes$stallion2[1:length(unique(cell.type.markers$Celltype))]
scales::show_col(annotation.colors)
names(annotation.colors) <- unique(cell.type.markers$Celltype)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Celltype, 
                                               levels = unique(cell.type.markers$Celltype)),
                                 col = list(Type = annotation.colors),show_annotation_name = F)
pdf("2.Cluster/Annotate/res_0.2/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 1.5), c("grey", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$Celltype, levels = unique(cell.type.markers$Celltype))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, left_annotation = row.annotations,
        width = unit(10, "cm"), height = unit(16, "cm"), cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))

a <- data.merge
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(T = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='T'],
                  B = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='B'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Plasma = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Plasma'],
                  Macro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Macro'],
                  DC = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='DC'],
                  Epi = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Epi'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo'],
                  Fibro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Fibro'])
library(GEB)
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "RNA_snn_res.0.2",dot.scale = 4)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')
####res = 0.4####
data.merge <- readRDS('./2.Cluster/data.merge.rds')
DefaultAssay(data.merge) <- "RNA"
pdf("2.Cluster/res_0.4/cluster.pdf")
data.merge@meta.data$seurat_clusters <- data.merge@meta.data$RNA_snn_res.0.4
length(unique(data.merge$seurat_clusters))
mycols_15 <- c("#00CDCD", "#A4D3EE", "#FFA500", "#FFB6C1", "#607B8B", "#8F8F8F", "#4F94CD",
               "#228B22", "#CD5C5C", "#EEDFCC", "#36648B", "#CD950C", "#009ACD", "#AB82FF", "#8B6969")
mycols_19 <- c("#66CDAA", "#458B74", "#7FFFD4", "#8EE5EE", "#7AC5CD", "#104E8B", "#009ACD", "#EEA2AD", 
                 "#CD8C95", "#EEDD82", "#CDBE70", "#FF7F00", "#CD6600", "#CD2626", "#8B7D6B", "#9F79EE", "#ADADAD", "#FFD700", "#9BCD9B")
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters",cols = mycols_19)+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "seurat_clusters",cols = mycols_19)+NoLegend()
DimPlot(object = data.merge,reduction = 'umap',label = TRUE,group.by = "seurat_clusters",cols = mycols_19,split.by = 'orig.ident',ncol = 3) + NoLegend()
dev.off()

#### Plot the ratio of each cell type and sample situation##
pdf("2.Cluster/res_0.4/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()
pdf("2.Cluster/res_0.4/sample.cell.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "seurat_clusters", id.vars2 = "orig.ident",color.len = mycols_19)
dev.off()


expMatrix <- GetAssayData(data.merge, slot = "scale.data")
highVariableGenes <- VariableFeatures(data.merge)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, data.merge$seurat_clusters, mean)
  return(mean.value)
})

pdf('2.Cluster/Annotate/res_0.4/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
cor <- (cor - min(cor))/(max(cor)-min(cor))
pheatmap::pheatmap(cor,angle_col = 45,cluster_cols = T,cluster_rows = T)
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2
library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("2.Cluster/res_0.4/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')
#### Differential expression####
data.merge <- readRDS('./2.Cluster/data.merge.rds')
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25,latent.vars = "orig.ident",test.use = 'MAST')
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "2.Cluster/Annotate/res_0.4/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "2.Cluster/Annotate/res_0.4/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "2.Cluster/Annotate/res_0.4/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
data.merge <- ScaleData(data.merge,features = c(VariableFeatures(data.merge),unique(top5.genes$gene)))
pdf("2.Cluster/Annotate/res_0.4/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(data.merge, features = unique(top5.genes$gene), size = 2,group.colors = mycols_15,group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

####The expression of the classic marker##
cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) ##perform 0-1 standardization for all clusters per gene
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$Celltype, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)##perform 0-1 standardization for all clusters per celltype marker
annotation.colors <- Palettes$stallion2[1:length(unique(cell.type.markers$Celltype))]
scales::show_col(annotation.colors)
names(annotation.colors) <- unique(cell.type.markers$Celltype)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Celltype, 
                                               levels = unique(cell.type.markers$Celltype)),
                                 col = list(Type = annotation.colors),show_annotation_name = F)
pdf("2.Cluster/Annotate/res_0.4/cluster.signature.expression.pdf",height = 15)
col_fun1 <- colorRamp2(c(0, 3.5), c("grey", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$Celltype, levels = unique(cell.type.markers$Celltype))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, left_annotation = row.annotations,
        width = unit(10, "cm"), height = unit(25, "cm"), cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
        title = "Expression", at = c(0, 1), 
        labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))

a <- data.merge
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(T = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='T'],
                  B = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='B'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Plasma = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Plasma'],
                  Macro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Macro'],
                  DC = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='DC'],
                  Epi = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Epi'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo'],
                  Fibro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Fibro'],
                  Smooth = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Smooth'])
library(GEB)
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "RNA_snn_res.0.4",dot.scale = 3)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')


#####annotated the cell type####
#major celltype####
data.merge <- readRDS('./2.Cluster/data.merge.rds')
table(Idents(data.merge))
table(data.merge$seurat_clusters)
major_celltype <- data.merge@meta.data$seurat_clusters
major_celltype <- gsub("^0$", "T", major_celltype)
major_celltype <- gsub("^1$", "T", major_celltype)
major_celltype <- gsub("^2$", "B", major_celltype)
major_celltype <- gsub("^3$", "B", major_celltype)
major_celltype <- gsub("^4$", "T", major_celltype)
major_celltype <- gsub("^5$", "T", major_celltype)
major_celltype <- gsub("^6$", "T", major_celltype)
major_celltype <- gsub("^7$", "B", major_celltype) 
major_celltype <- gsub("^8$", "Myeloid", major_celltype)
major_celltype <- gsub("^9$", "Myeloid", major_celltype)
major_celltype <- gsub("^10$", "Cycling", major_celltype)
major_celltype <- gsub("^11$", "NK", major_celltype)
major_celltype <- gsub("^12$", "T", major_celltype)
major_celltype <- gsub("^13$", "B", major_celltype)
major_celltype <- gsub("^14$", "Endo", major_celltype)
major_celltype <- gsub("^15$", "Myeloid", major_celltype)
major_celltype <- gsub("^16$", "Fibro", major_celltype)
major_celltype <- gsub("^17$", "Myeloid", major_celltype)
major_celltype <- gsub("^18$", "Smooth", major_celltype)
table(major_celltype)
data.merge <- AddMetaData(data.merge, major_celltype, col.name = "major_celltype")
table(data.merge$major_celltype)
data.merge@meta.data$large_annotation <- ifelse(data.merge@meta.data$major_celltype %in% c('T','B','NK','Myeloid','Cycling'),'Immune',
                                                ifelse(data.merge@meta.data$major_celltype%in%c('Endo','Fibro','Smooth'),'Stromal',""))
table(data.merge$large_annotation)
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')

####verify the annotated results by Celltypist####
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
library(reticulate)
py_config()
scanpy <- import('scanpy')
pandas <- import('pandas')
numpy <- import('numpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(data.merge[['RNA']]@counts)))),
                       obs = pandas$DataFrame(data.merge@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(data.merge[['RNA']]@counts),
                                                         row.names = rownames(data.merge[['RNA']]@counts)))
)
adata
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Immune_All_High.pkl', majority_voting = T)
predicted_major <- as.data.frame(predictions$predicted_labels)
predicted_major$seurat_cluster <- data.merge@meta.data$seurat_clusters
data.merge <- AddMetaData(data.merge,predicted_major$majority_voting,col.name = 'celltypist_major')
table(data.merge$celltypist_major,data.merge$seurat_clusters)

predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- data.merge@meta.data$seurat_clusters
data.merge <- AddMetaData(data.merge,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(data.merge$celltypist_minor,data.merge$seurat_clusters)

write.xlsx(predicted_major,file = './2.Cluster/Annotate/res_0.4/celltypist_predicted.major.xlsx',rowNames = T)
write.xlsx(predicted_minor,file = './2.Cluster/Annotate/res_0.4/celltypist_predicted.minor.xlsx',rowNames = T)

View(data.merge@meta.data)
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')

####observe the annotated results####
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
source('./function/colorPalettes.R')
pdf("2.Cluster/Annotate/res_0.4/cellType.pro.pdf")
length(unique(data.merge$major_celltype))
DimPlot(object = data.merge, reduction = 'tsne',label = FALSE, group.by = "large_annotation")
DimPlot(object = data.merge, reduction = 'umap',label = FALSE, group.by = "large_annotation")
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "major_celltype",cols = Palettes[['mycols_8']])+NoLegend()
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "major_celltype",cols = Palettes[['mycols_8']])
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "major_celltype",cols = Palettes[['mycols_8']])+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "major_celltype",cols = Palettes[['mycols_8']])

dev.off()

##Plot--- celltype marker plot
pdf("2.Cluster/Annotate/res_0.4/cellType.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "major_celltype", angle = 60)
dev.off()

pdf("2.Cluster/Annotate/res_0.4/cellType.sample.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "major_celltype", id.vars2 = "orig.ident", angle = 60,color.len = Palettes[['mycols_8']])
dev.off()

cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres_v2.txt", header = T, stringsAsFactors = F, sep = "\t")
exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(
                  B = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='B'],
                  T = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='T'],
                  Cycing = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Proliferating'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo'],
                  Fibro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Fibro'],
                  Myeloid = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Myeloid'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Smooth = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Smooth']
                  )
data.merge@meta.data$major_celltype <- factor(data.merge@meta.data$major_celltype,levels = c(
  'B','T','Cycling','Endo','Fibro','Myeloid','NK','Smooth'
),ordered = T)
library(GEB)
pdf("2.Cluster/Annotate/res_0.4/cellType.marker.dotplot.pdf", height = 10,width = 8)
DotPlot_ByColumnList(data.merge,group.by = 'major_celltype',cols = c('grey','red'),features = gene_list,scale.by = 'size',col.min = 1.5)
dev.off()

group <- ifelse(grepl(pattern = 'T',data.merge@meta.data$sample),'tumor','normal')
data.merge <- AddMetaData(data.merge,metadata = group,col.name = 'group')
table(data.merge$group)
pdf("2.Cluster/Annotate/res_0.4/cellType.group.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "group", id.vars2 = "major_celltype", angle = 60,color.len = Palettes[['mycols_8']])
dev.off()

pdf("2.Cluster/Annotate/res_0.4/group.cellratio.pdf", height = 10, width = 4)
ratio.plot(seurat.object = data.merge, id.vars1 = "major_celltype", id.vars2 = "group", angle = 60,color.len = Palettes[['mycols_8']])
dev.off()
####
library(Nebulosa)
pdf('./2.Cluster/Annotate/res_0.4/T.pdf',width = 16,height = 16)
plot_density(data.merge,features = c('CD3D','CD3E','CD3G','CD7'),joint = F,reduction = 'umap')
plot_density(data.merge,features = c('CD3D','CD3E','CD3G','CD7'),joint = F,reduction = 'tsne')

dev.off()
####
pdf('./2.Cluster/Annotate/res_0.4/B.pdf',width = 16,height = 16)
plot_density(data.merge,features = c('MS4A1','CD19','CD79A','CD79B'),joint = F,reduction = 'umap',begin = 0.2)
plot_density(data.merge,features = c('MS4A1','CD19','CD79A','CD79B'),joint = F,reduction = 'tsne',begin = 0.2)
VlnPlot(data.merge,features = c('CD79A','CD79B','MS4A1'),group.by = 'major_celltype',pt.size = 0)
VlnPlot(data.merge,features = c('CD79A','CD79B','MS4A1'),group.by = 'seurat_clusters',pt.size = 0)
dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/NK.pdf',width=16,height = 16)
plot_density(data.merge,features = c('NKG7','CST7','KLRB1','GNLY'),joint = F,reduction = 'umap')
plot_density(data.merge,features = c('NKG7','CST7','KLRB1','GNLY'),joint = F,reduction = 'tsne')

dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/plasma.pdf',width=16,height = 16)
plot_density(data.merge,features = c('JCHAIN','SDC1','MZB1','IGHG1'),joint = F,reduction = 'umap')
plot_density(data.merge,features = c('JCHAIN','SDC1','MZB1','IGHG1'),joint = F,reduction = 'tsne')

dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/macro.pdf',width=16,height = 16)
plot_density(data.merge,features = c('CD14','CD68','LYZ','C1QA'),joint = F,reduction = 'umap')
plot_density(data.merge,features = c('CD14','CD68','LYZ','C1QA'),joint = F,reduction = 'tsne')


dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/DC.pdf',width=16,height = 16)
plot_density(data.merge,features = c('CD1C','CLEC9A','BIRC3','LAMP3'),joint = F,reduction = 'umap')
plot_density(data.merge,features = c('CD1C','CLEC9A','BIRC3','LAMP3'),joint = F,reduction = 'tsne')


dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/cluster6_CD8CTL.pdf',width = 16,height = 16)
plot_density(data.merge,features = c('CD8A','CD8B','NKG7','GZMK'),reduction = 'umap',joint = F)
plot_density(data.merge,features = c('CD8A','CD8B','NKG7','GZMK'),reduction = 'tsne',joint = F)
dev.off()
####
pdf('./2.Cluster/Annotate/res_0.4/T_featureplot.pdf',width = 16,height = 16)
FeaturePlot(data.merge,features = c('CD3D','CD3E','CD3G','CD7'),reduction = 'umap',order = T)
dev.off()
####
pdf('2.Cluster/Annotate/res_0.4/imm.pdf',width = 8,height = 8)
plot_density(data.merge,features = c('PTPRC'),reduction = 'umap',joint = F)
plot_density(data.merge,features = c('PTPRC'),reduction = 'tsne',joint = F)
VlnPlot(data.merge,features = 'PTPRC',group.by = 'major_celltype',pt.size = 0)
dev.off()
####
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')

#### Cell type specific gene####
data.merge <- readRDS('2.Cluster/data.merge.pro.rds')
Idents(data.merge) <- data.merge$major_celltype
idents <- as.character(levels(data.merge))
cellType.all.markers <- FindAllMarkers(data.merge, 
                                       group.by = "major_celltype", 
                                       logfc.threshold = 0.25, 
                                       min.pct = 0.25, 
                                       test.use = "MAST", 
                                       latent.vars = "orig.ident")##同时包含上调和下调基因
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.all.markers$cluster == x)
  DEGs <- cellType.all.markers[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/Annotate/res_0.4/celltype.all.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.all.markers, file = "2.Cluster/Annotate/res_0.4/cellType.all.DEGs.rds")

#require logfc.threshold = 0.25 & p_val_adj < 0.05
cellType.sig.DEGs <- cellType.all.markers %>% filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.sig.DEGs$cluster == x)
  DEGs <- cellType.sig.DEGs[index,]
  DEGs <- DEGs %>% arrange(desc(avg_log2FC))
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/Annotate/res_0.4/cellType.sig.pos.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.sig.DEGs, file = "2.Cluster/Annotate/res_0.4/cellType.sig.pos.DEGs.rds")
top.genes <- cellType.sig.DEGs %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
top.genes <- top.genes[order(top.genes$cluster),]
data.merge <- ScaleData(data.merge,features = c(VariableFeatures(data.merge),unique(top.genes$gene)))
pdf("2.Cluster/Annotate/res_0.4/cellType.topgenes.pdf",width = 15)
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend() + scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()
source('./function/do_heatmap.R')
pdf('./2.Cluster/Annotate/res_0.4/cellType.topgenes.mean.pdf')
DoHeatmap_my(data.merge,features = unique(top.genes$gene),slot = 'scale.data',cluster_cols = F,color = colorRampPalette(Palettes[['blueRed']])(100),angle_col = 45,
             assay = 'RNA',group.by = 'major_celltype',cluster_rows = F)

dev.off()

####cell lineage ratio####
cellType.ratio <- as.data.frame(table(data.merge$major_celltype))
cellType.ratio$Type <- rep("Lymphoid cells", nrow(cellType.ratio))
idx <- which(cellType.ratio$Var1 %in% c("Myeloid"))
cellType.ratio$Type[idx] <- "Myeloid cells"
idx <- which(cellType.ratio$Var1 %in% c("Smooth","Endo",'Fibro'))
cellType.ratio$Type[idx] <- "Stromal cells"
group.ratio <- tapply(cellType.ratio$Freq, cellType.ratio$Type, sum)
group.ratio <- data.frame(Type = names(group.ratio), Ratio = group.ratio)
labs <- paste0(group.ratio$Type, " (", round(group.ratio$Ratio/sum(group.ratio$Ratio), 4)*100, "%)")
pdf("2.Cluster/Annotate/res_0.4/lineage.ratio.pdf")
p <- ggpie(group.ratio, "Ratio", label = labs,
           fill = "Type", color = "white", lab.pos = "out",
           palette = "jama")
print(p)
dev.off()

#### The proportion of various cell types in each patient####
library(plotly)
pdf("2.Cluster/Annotate/res_0.4/cellType.patient.ratio.pdf")
patients <- gsub(x = data.merge@meta.data$sample,pattern = '[A-Z]',"")
data.merge <- AddMetaData(data.merge,metadata = patients,col.name = 'patients')
index <- unique(data.merge$patients)
res <- lapply(index, function(x){
  a <- subset(data.merge, subset = patients==x)
  cellType.ratio <- as.data.frame(table(as.character(a$major_celltype)))
  group.ratio <- cellType.ratio$Freq/sum(cellType.ratio$Freq)
  group.ratio <- data.frame(Type = cellType.ratio$Var1, Ratio = group.ratio)
  group.ratio$Type = paste(group.ratio$Type,'(',round(group.ratio$Ratio, 4)*100, "%",')')
  group.ratio$labs = ""
  p <- ggpie(group.ratio, "Ratio", label = "labs",
             fill = "Type", color = "white", lab.pos = "out",palette = Palettes[['mycols_8_2']])
  p <- ggpar(p,title = x)
  print(p)
})
dev.off()
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')
##########
