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
plan("multisession", workers = 1) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
library(plot1cell)
#### Harmony corrected result####
####54136cells
data.merge <- readRDS('./2.Cluster/filtered/data.merge.filtered.pc20.seed45.rds')
DefaultAssay(data.merge) <- "RNA"
####res = 0.4####
dir.create('./2.Cluster/Annotate')
data.merge@meta.data$seurat_clusters <- data.merge@meta.data$RNA_snn_res.0.4
length(unique(data.merge$seurat_clusters))
table(data.merge$seurat_clusters)
pdf("2.Cluster/Annotate/cluster.pdf")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "seurat_clusters",cols = Palettes[['circus']])+NoLegend()
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "seurat_clusters",cols = Palettes[['circus']])+NoLegend()
DimPlot(object = data.merge,reduction = 'umap',label = TRUE,group.by = "seurat_clusters",cols = Palettes[['circus']],split.by = 'orig.ident',ncol = 3) + NoLegend()
DimPlot(object = data.merge,reduction = 'tsne',label = TRUE,group.by = "seurat_clusters",cols = Palettes[['circus']],split.by = 'orig.ident',ncol = 3) + NoLegend()
dev.off()

pdf("2.Cluster/Annotate/cluster.umap.pdf",width = 8,height = 8)
SCpubr::do_DimPlot(data.merge,group.by = 'seurat_clusters',reduction = 'umap',label = T,font.size = 15,pt.size = .3,border.size = 0,plot.axes = T) + 
  scale_color_manual(values = Palettes[['circus']]) + NoLegend()
dev.off()
## Plot the ratio of each cell type and sample situation##
pdf("2.Cluster/Annotate/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()
pdf("2.Cluster/Annotate/sample.cell.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "seurat_clusters", id.vars2 = "orig.ident",color.len = Palettes[['mycols_15']])
dev.off()


expMatrix <- GetAssayData(data.merge, slot = "scale.data")
highVariableGenes <- VariableFeatures(data.merge)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, data.merge$seurat_clusters, mean)
  return(mean.value)
})

pdf('2.Cluster/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
               N63+tmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2
library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("2.Cluster/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')
#### Differential expression##
data.merge <- readRDS('./2.Cluster/data.merge.rds')
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "2.Cluster/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "2.Cluster/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "2.Cluster/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("2.Cluster/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(subset(data.merge,downsample = 50), features = unique(top5.genes$gene), size = 2,group.bar = T) + NoLegend()
dev.off()

####The expression of the classic marker##
data.merge <- readRDS('./2.Cluster/data.merge.rds')
pdf('./2.Cluster/Annotate/observe.expression.T.B.pdf',width = 10,height = 8)
FeaturePlot(data.merge,features = c('MS4A1','CD19','CD79A','CD79B'),reduction ='umap',order = T,min.cutoff = 1.3,cols = Palettes[['greyMagma']])
FeaturePlot(data.merge,features = c('CD7','CD3D','CD3E','CD3G'),reduction ='umap',order = T,min.cutoff = 1.3,cols = Palettes[['greyMagma']])
dev.off()


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
names(annotation.colors) <- unique(cell.type.markers$Celltype)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Celltype, 
                                               levels = unique(cell.type.markers$Celltype)),
                                 col = list(Type = annotation.colors),show_annotation_name = F)
pdf("2.Cluster/Annotate/cluster.signature.expression.pdf")
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
                  cycling = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Cycling'],
                  Plasma = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Plasma'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Myeloid = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Myeloid'],
                  Epi = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Epi'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo'],
                  Fibro = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Fibro'],
                  Smooth = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Smooth'])
library(GEB)
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "RNA_snn_res.0.4",dot.scale = 4)
dev.off()

saveRDS(data.merge,file = './2.Cluster/data.merge.rds')
#####annotated the cell type####
#major celltype####
data.merge <- readRDS('./2.Cluster/data.merge.rds')
table(Idents(data.merge))
table(data.merge$seurat_clusters)
major_celltype <- data.merge@meta.data$seurat_clusters
major_celltype <- gsub("^0$", "T", major_celltype)
major_celltype <- gsub("^1$", "B", major_celltype)
major_celltype <- gsub("^2$", "T", major_celltype)
major_celltype <- gsub("^3$", "B", major_celltype)
major_celltype <- gsub("^4$", "T", major_celltype)
major_celltype <- gsub("^5$", "T", major_celltype)
major_celltype <- gsub("^6$", "T", major_celltype)
major_celltype <- gsub("^7$", "T", major_celltype) 
major_celltype <- gsub("^8$", "B", major_celltype)
major_celltype <- gsub("^9$", "Myeloid", major_celltype)
major_celltype <- gsub("^10$", "NK", major_celltype)
major_celltype <- gsub("^11$", "Myofibroblast", major_celltype)
major_celltype <- gsub("^12$", "LAMP3+cDC", major_celltype)
major_celltype <- gsub("^13$", "pDC", major_celltype)
major_celltype <- gsub("^14$", "Endo", major_celltype)
table(major_celltype)
data.merge <- AddMetaData(data.merge, major_celltype, col.name = "major_celltype")
table(data.merge$major_celltype)

pdf('./2.Cluster/Annotate/Imm.observe.pdf')
VlnPlot(data.merge,features = 'PTPRC',group.by = 'major_celltype',pt.size = 0,cols = Palettes[['circus']])
dev.off()

pdf('./2.Cluster/Annotate/Imm.dimplot.pdf',width = 8,height = 8)

data.merge$immune <- ifelse(data.merge$large_annotation=='Immune','Immune','Non-immune')
SCpubr::do_DimPlot(data.merge,group.by = 'immune',reduction = 'umap',border.size = 0,plot.axes = T) + 
  theme(legend.text = element_text(size = 15))
dev.off()

data.merge@meta.data$large_annotation <- ifelse(data.merge@meta.data$major_celltype %in% c('T','B','NK','Myeloid','LAMP3+cDC','pDC'),'Immune',
                                                data.merge@meta.data$major_celltype)
table(data.merge$large_annotation)
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')
####use plot1cell to visualize####
library(plot1cell)
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
table(Idents(data.merge))
circ_data <- prepare_circlize_data(data.merge, scale = 0.8)
cluster_colors<- Palettes[['circus']][1:length(unique(Idents(data.merge)))]
group_colors<-rand_color(length(names(table(data.merge$group))))
rep_colors<-rand_color(length(names(table(data.merge$orig.ident))))
pdf('2.Cluster/Annotate/circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "sample",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

pdf('./2.Cluster/Annotate/celltype_umap_clean.pdf',6,6)
DimPlot(data.merge,group.by = 'major_celltype',cols = Palettes[['circus']],label = F) + NoLegend()
dev.off()
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

write.xlsx(predicted_major,file = './2.Cluster/Annotate/celltypist_predicted.major.xlsx',rowNames = T)

View(data.merge@meta.data)
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')


####observe the annotated results####
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
data.merge@meta.data$major_celltype <- factor(data.merge@meta.data$major_celltype,levels = c(
  'T','B','NK','Myeloid','LAMP3+cDC','pDC','Myofibroblast','Endo'
),ordered = T)
source('./function/colorPalettes.R')
source('./function/do_dimplot.R')
pdf("2.Cluster/Annotate/cellType.pro.pdf")
length(unique(data.merge$major_celltype))
DimPlot(object = data.merge, reduction = 'tsne',label = FALSE, group.by = "large_annotation",cols = Palettes[['mycols_4']])
DimPlot(object = data.merge, reduction = 'umap',label = FALSE, group.by = "large_annotation",cols = Palettes[['mycols_4']])
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "major_celltype",cols = Palettes[['circus']])+NoLegend()
DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "major_celltype",cols = Palettes[['circus']])
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "major_celltype",cols = Palettes[['circus']])+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "major_celltype",cols = Palettes[['circus']])
dev.off()

pdf('./2.Cluster/Annotate/cellType.pro.umap.pdf',width = 9,height = 8)
DoDimplot(data.merge,groupby = 'major_celltype',colors = Palettes[['circus']]) + NoLegend()
dev.off()

library(SCpubr)
cols <- c("#b60e10","#233341","#f4b9b7","#6098b6","#987143","#d6828c","#c9ad99","#3e9896")
names(cols) <- levels(data.merge)
pdf('./2.Cluster/Annotate/celltype.umap.pro.pdf',8,8)
SCpubr::do_DimPlot(data.merge,reduction = 'umap',group.by = 'major_celltype',label = T,colors.use = cols,
                  plot_density_contour = T,pt.size = .6,contour.color = 'black')
dev.off()
Q##Plot--- celltype marker plot
pdf("2.Cluster/Annotate/cellType.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "major_celltype", angle = 60,color.len = Palettes[['circus']])
dev.off()

pdf("2.Cluster/Annotate/cellType.sample.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "major_celltype", id.vars2 = "orig.ident", angle = 60,color.len = Palettes[['circus']])
dev.off()

cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres_v2.txt", header = T, stringsAsFactors = F, sep = "\t")
exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(
                  T = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='T'],
                  B = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='B'],
                  NK = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='NK'],
                  Myeloid = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Myeloid'],
                  'LAMP3+cDC' = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='LAMP3+cDC'],
                  pDC = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='pDC'],
                  Myofibroblast = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Myofibroblast'],
                  Endo = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Celltype=='Endo']
                  )
pdf("2.Cluster/Annotate/cellType.marker.dotplot.pdf", height = 5,width = 18)
Idents(data.merge) <- data.merge$major_celltype
DotPlot(data.merge,features = gene_list,scale.by = 'size',col.min = 1) + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 12,angle = 25,hjust = 1))
dev.off()

Idents(data.merge) <- data.merge$major_celltype
exp.matrix <- GetAssayData(data.merge, slot = "data")
cell.type.markers_distinct$Celltype <- factor(cell.type.markers_distinct$Celltype,levels = levels(data.merge),ordered = T)
cell.type.markers_distinct <- cell.type.markers_distinct[order(cell.type.markers_distinct$Celltype),]
index <- match(cell.type.markers_distinct$Gene[!is.na(cell.type.markers_distinct$Celltype)], rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$major_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('2.Cluster/Annotate/cellType.marker.heatmap.pdf',width = 22,height = 12)
Heatmap((cluster.score.normailzed), 
        width = unit(28, "cm"), height = unit(8, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

group <- ifelse(grepl(pattern = 'T',data.merge@meta.data$sample),'tumor','normal')
data.merge <- AddMetaData(data.merge,metadata = group,col.name = 'group')
table(data.merge$group)
pdf("2.Cluster/Annotate/cellType.group.ratio.pdf", height = 10, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "group", id.vars2 = "major_celltype",angle = 60) + 
  scale_fill_manual(values = c("#FF4040", "#1874CD")) + theme(text = element_text(size = 24))
dev.off()

pdf("2.Cluster/Annotate/group.cellratio.pdf", height = 10, width = 8)
ratio.plot(seurat.object = data.merge, id.vars1 = "major_celltype", id.vars2 = "group", angle = 60,color.len = Palettes[['circus']])
dev.off()
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')


Idents(data.merge) <- data.merge$major_celltype
table(Idents(data.merge))
# 每个细胞亚群抽1/20
allCells=names(Idents(data.merge))
allType = levels(Idents(data.merge))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(data.merge)== x ]
  if(length(cgCells) > 1000){
    cg=sample(cgCells,size = floor(length(cgCells)/20))}
  else{
    cg = cgCells
  }
  cg
}))

use_colors <- Palettes[['circus']][1:length(unique(data.merge$major_celltype))]
cg_sce = data.merge[,allCells %in% choose_Cells]
table(Idents(cg_sce))
library(ggSankeyGrad)
tmp <- as.data.frame(table(cg_sce$major_celltype,cg_sce$sample))
cols1 <- unlist(lapply(use_colors,function(i){
  rep(i,6)
}))
use_colors <- Palettes[['circus']][1:length(unique(data.merge$sample))]
cols2 <- rep(use_colors,8)
tmp$cols1 = cols1
tmp$cols2 = cols2
pdf('./2.Cluster/Annotate/sample.celltype.sankeyplot.pdf')
with(tmp,ggSankeyGrad(c1 = Var1,c2 = Var2,values = tmp$Freq,col1 = cols1,col2 = cols2,
             label = T,alpha = 0.6,padding = 30))
dev.off()

tmp <- as.data.frame(table(cg_sce$major_celltype,cg_sce$group))
use_colors <- c("#003f5c","#58508d","#8a508f","#bc5090","#de5a79","#ff6361","#ff8531","#ffa600")

cols1 <- unlist(lapply(use_colors,function(i){
  rep(i,2)
}))
cols2 <- cols1
cols2 <- rep(c("#de5357","#0087b0"),8)
tmp$cols1 = cols1
tmp$cols2 = cols2
tmp$Freq[16] <- 15
pdf('./2.Cluster/Annotate/group.celltype.sankeyplot.pdf')
with(tmp,ggSankeyGrad(c1 = Var1,c2 = Var2,values = tmp$Freq,col1 = cols1,col2 = cols2,
                      label = T,alpha = 0.9,padding = 100,label_size = 16))
dev.off()
pdf('2.Cluster/Annotate/dimplot.group.pdf')
DimPlot(data.merge,group.by = 'group',reduction = 'umap',cols = Palettes[['circus']],label = F)
DimPlot(data.merge,group.by = 'group',reduction = 'tsne',cols = Palettes[['circus']],label = F)
dev.off()
pdf('2.Cluster/Annotate/dimplot.group.split.pdf',width = 12,height = 6)
DimPlot(data.merge,group.by = 'major_celltype',reduction = 'umap',cols = Palettes[['circus']],label = F,split.by = 'group')
dev.off()

cols <- c("#fbee5a","#f6cc5e","#f1a962","#ec8766","#e7656a","#e2436e")
names(cols)=unique(data.merge$sample)
pdf('2.Cluster/Annotate/dimplot.sample.pdf',6,6)
SCpubr::do_DimPlot(data.merge,reduction = 'umap',group.by = 'sample',border.size = 1,label = F,colors.use = cols,
                   pt.size = 5,raster = T,raster.dpi = 2048)
dev.off()
#### Cell type specific gene####
data.merge <- readRDS('2.Cluster/data.merge.pro.rds')
Idents(data.merge) <- data.merge$major_celltype
idents <- as.character(levels(data.merge))
cellType.all.markers <- FindAllMarkers(data.merge, 
                                       group.by = "major_celltype", 
                                       logfc.threshold = 0.25, 
                                       min.pct = 0.25,
                                       )##同时包含上调和下调基因
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.all.markers$cluster == x)
  DEGs <- cellType.all.markers[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/Annotate/celltype.all.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.all.markers, file = "2.Cluster/Annotate/cellType.all.DEGs.rds")

#require logfc.threshold = 0.25 & p_val_adj < 0.05
cellType.sig.DEGs <- cellType.all.markers %>% filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.sig.DEGs$cluster == x)
  DEGs <- cellType.sig.DEGs[index,]
  DEGs <- DEGs %>% arrange(desc(avg_log2FC))
  return(DEGs)
})
write.xlsx(saveFormat, file = "2.Cluster/Annotate/cellType.sig.pos.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.sig.DEGs, file = "2.Cluster/Annotate/cellType.sig.pos.DEGs.rds")
top.genes <- cellType.sig.DEGs %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
top.genes <- top.genes[order(top.genes$cluster),]
data.merge <- ScaleData(data.merge,features = c(VariableFeatures(data.merge),unique(top.genes$gene)))
pdf("2.Cluster/Annotate/cellType.topgenes.pdf",width = 15)
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend() + scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()
source('./function/do_heatmap.R')
pdf('./2.Cluster/Annotate/cellType.topgenes.mean.pdf')
DoHeatmap_my(data.merge,features = unique(top.genes$gene),slot = 'scale.data',cluster_cols = F,color = colorRampPalette(Palettes[['blueRed']])(100),angle_col = 45,
             assay = 'RNA',group.by = 'major_celltype',cluster_rows = F)

dev.off()

#### The proportion of various cell types in each patient####
library(plotly)
pdf("2.Cluster/Annotate/cellType.patient.ratio.pdf")
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
             fill = "Type", color = "white", lab.pos = "out",palette = Palettes[['mycols_8']])
  p <- ggpar(p,title = x)
  print(p)
})
dev.off()
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')
####test of celltype percentage between groups####
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
df <- as.data.frame(table(data.merge$major_celltype,data.merge$sample))
df$total <- apply(df,1,function(x){sum(df$Freq[df$Var2==x[2]])})
df$percent <- round(df$Freq/df$total,4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('./2.Cluster/Annotate/celltype.percentage.group.pdf',height = 6,width = 12)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',color = 'group',palette = 'npg',width = 0,line.size = .4,line.color = 'black') +
    scale_y_continuous(breaks = scales::pretty_breaks(n=5)) + 
    theme(axis.text.x = element_text(size = 15)) + 
    labs(title = i,x = '',y = '% of all cells') + NoLegend()
  return(res)
})
plot_grid(plotlist = res,ncol = 4)
dev.off()

library(speckle)
metadata <- data.merge@meta.data
res <- propeller(clusters = metadata$major_celltype,sample = metadata$sample,group = metadata$group)
write.csv(res,file = './2.Cluster/Annotate/celltype.propeller.test.csv')
saveRDS(data.merge,file = './2.Cluster/data.merge.pro.rds')
