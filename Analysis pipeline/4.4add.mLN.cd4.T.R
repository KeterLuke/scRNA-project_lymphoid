####add mLN of GSE131907 CD4T to merged CD4T
####load packages####
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
getwd()
setwd('../')
getwd()
source('./function/colorPalettes.R')
########
dir.create('./4.4add.mLN.cd4.T')
integrated.cd4 <- readRDS('./4.3add.normal.cd4.T/CD4.merge.rds')
table(integrated.cd4$paper)
####process GSE131907 mLn data and annotated major celltypes####
####read mLn samples of GSE131907
library(data.table)
obj=readRDS("./raw/extra data/normal_lym/GSE131907_Lung_Cancer_raw_UMI_matrix.rds/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
ana=fread("./raw/extra data/normal_lym/GSE131907_Lung_Cancer_cell_annotation.txt/GSE131907_Lung_Cancer_cell_annotation.txt")
table(ana$Sample_Origin)
ana=ana[ana$Sample_Origin=="mLN",]
table(ana$Sample)
mat=obj[,ana$Index]
colnames(mat)=paste0(colnames(mat),"-",'GSE131907_mLN')
GSE131907_mLN=CreateSeuratObject(counts = mat, project = 'GSE131907_mLN', min.cells = 3, min.features = 200)
GSE131907_mLN$source="GSE131907_mLN"
GSE131907_mLN <- GSE131907_mLN %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(GSE131907_mLN,harmony = F)
GSE131907_mLN  <- FindNeighbors(object = GSE131907_mLN , dims = 1:13, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
GSE131907_mLN  <- FindClusters(object = GSE131907_mLN , resolution = set.resolutions, verbose = FALSE)
GSE131907_mLN  <- RunUMAP(GSE131907_mLN , dims = 1:13)
pdf(file = "4.4add.mLN.cd4.T/GSE131907_mLN_PCA-test.pdf")
clustree(GSE131907_mLN)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = GSE131907_mLN, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.1
saveRDS(GSE131907_mLN,file = '4.4add.mLN.cd4.T//GSE131907_mLN.rds')
####remove doublet
source('./function/doubletDetect.R')
dir.create('./4.4add.mLN.cd4.T/doublet')
####GSE131907_nLN
sce <- readRDS('./4.4add.mLN.cd4.T/GSE131907_mLN.rds')
pdf("4.4add.mLN.cd4.T/doublet/GSE131907_mLN_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:13, doublet.rate = 0.1, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.1
saveRDS(sce,file = '4.4add.mLN.cd4.T/doublet/GSE131907_mLN.RDS')

mLn <- readRDS('./4.4add.mLN.cd4.T/doublet/GSE131907_mLN.RDS')
table(mLn$Doublet)
mLn <- subset(mLn,Doublet=='Singlet')
mLn <- mLn %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(mLn,harmony = F)
pc <- 55
mLn <- FindNeighbors(mLn,reduction = 'pca',dims = 1:pc,verbose = F)
set.resolutions <- seq(0.1,1,0.1)
mLn <- FindClusters(mLn,resolution = set.resolutions,verbose = F)
mLn <- RunUMAP(mLn,reduction = 'pca',dims = 1:pc,verbose = F)
pdf('./4.4add.mLN.cd4.T/mLn.res.observe.pdf')
clustree(mLn)
DimPlot(mLn,reduction = 'umap',group.by  = 'source')

sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = mLn, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()

####select res 0.3
mLn$seurat_clusters <- mLn$RNA_snn_res.0.3
table(mLn$seurat_clusters)
Idents(mLn) <- mLn$seurat_clusters
pdf('./4.4add.mLN.cd4.T/marker.observe.pdf')
FeaturePlot(mLn,features = c('CD3D','CD3E','CD79A','MS4A1'),order = T,min.cutoff = 1)
VlnPlot(mLn,features = c('CD3D','CD3E','CD79A','MS4A1'),pt.size = 0,ncol = 2)
FeaturePlot(mLn,features = c('PTPRC','EPCAM','CD3D','MS4A1'),order = T,min.cutoff = 1)
dev.off()
pdf('./4.4add.mLN.cd4.T/mLN.dimplot.pdf')
DimPlot(mLn,group.by = 'seurat_clusters',reduction = 'umap',label = T)
DimPlot(mLn,group.by = 'source',reduction = 'umap')
dev.off()
#####expression of classical markers
cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")
df <- data.frame(Celltype = 'Epi',Gene = c('EPCAM','KRT13','KRT18','KRT19'))
cell.type.markers <- rbind(cell.type.markers,df)
exp.matrix <- GetAssayData(mLn, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, mLn$seurat_clusters, mean)
  return(a)
})
library(vegan)
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
pdf("4.4add.mLN.cd4.T/cluster.signature.expression.pdf")
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
a <- mLn
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
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "seurat_clusters",dot.scale = 4)
dev.off()

plan("multisession",workers = 12)
cluster.pos.markers <- FindAllMarkers(mLn, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
plan('multisession',workers = 1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.4add.mLN.cd4.T/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.4add.mLN.cd4.T/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.4add.mLN.cd4.T/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)
####annotate
table(mLn$seurat_clusters)
major_celltype <- mLn@meta.data$seurat_clusters
major_celltype <- gsub("^0$", "B", major_celltype)
major_celltype <- gsub("^1$", "T", major_celltype)
major_celltype <- gsub("^2$", "Myeloid", major_celltype)
major_celltype <- gsub("^3$", "T", major_celltype)
major_celltype <- gsub("^4$", "Myeloid", major_celltype)
major_celltype <- gsub("^5$", "NK", major_celltype)
major_celltype <- gsub("^6$", "Epi", major_celltype)
major_celltype <- gsub("^7$", "Epi", major_celltype) 
major_celltype <- gsub("^8$", "Myeloid", major_celltype)
major_celltype <- gsub("^9$", "Cycling", major_celltype)
major_celltype <- gsub("^10$", "Stromal", major_celltype)
major_celltype <- gsub("^11$", "T", major_celltype)
major_celltype <- gsub("^12$", "Plasma", major_celltype)
major_celltype <- gsub("^13$", "Epi", major_celltype)
major_celltype <- gsub("^14$", "Epi", major_celltype)
major_celltype <- gsub("^15$", "Myeloid", major_celltype)
major_celltype <- gsub("^16$", "Epi", major_celltype)
major_celltype <- gsub("^17$", "Epi", major_celltype)
major_celltype <- gsub("^18$", "Mast", major_celltype)
major_celltype <- gsub("^19$", "Epi", major_celltype)
major_celltype <- gsub("^20$", "pDC", major_celltype)


table(major_celltype)
mLn <- AddMetaData(mLn, major_celltype, col.name = "major_celltype")
table(mLn$major_celltype)
saveRDS(mLn,file = './4.4add.mLN.cd4.T/annotated.GSE131907.mLN.rds')

mLn <- readRDS('./4.4add.mLN.cd4.T/annotated.GSE131907.mLN.rds')
pdf('./4.4add.mLN.cd4.T/mLn.major_celltype.pdf')
DimPlot(mLn,group.by = 'major_celltype',reduction = 'umap',label = T)
dev.off()

table(mLn$major_celltype)
genelist=c('MS4A1','CD19','CD79A','CD79B','MKI67','TOP2A','STMN1','KIAA0101','EPCAM','KRT18','SFTPD','CEACAM6',
           'TPSB2','TPSAB1','CPA3','MS4A2','APOC1','LYZ','CD14','CD68','NKG7','KLRB1',
           'GNLY','TRDC','LILRA4','IRF7','CXCR3','TCF4','MZB1','IGHG1','IGHG2','IGHG3','COL11A1','COL3A1','MYL9','TUBB2B',
           'CD3D','CD3E','CD3G','CD7')
mLn$major_celltype <- factor(mLn$major_celltype,levels = c('B','Cycling','Epi','Mast','Myeloid','NK','pDC','Plasma','Stromal','T'))
pdf("4.4add.mLN.cd4.T/major.cellType.marker.dotplot.pdf", height = 5,width = 12)
Idents(mLn) <- mLn$major_celltype
DotPlot(mLn,features = genelist,scale.by = 'size',col.min = 0) + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 8,angle = 25))
dev.off()

Idents(mLn) <- mLn$major_celltype
exp.matrix <- GetAssayData(mLn, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, mLn$major_celltype, mean)
  return(a)
})
library(vegan)
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('4.4add.mLN.cd4.T/major.cellType.marker.heatmap.pdf')
Heatmap(t(cluster.score.normailzed), 
        column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =15),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

table(mLn$major_celltype)
sub.T <- subset(mLn,major_celltype=='T')
FeaturePlot(sub.T,features = c('CD3D','CD3E','CD3G','MS4A1'),order = T,min.cutoff = 0)
sub.T <- DietSeurat(sub.T,assays = 'RNA')
saveRDS(sub.T,file = './4.4add.mLN.cd4.T/mLn.sub.T.rds')
saveRDS(mLn,file = './4.4add.mLN.cd4.T/annotated.GSE131907.mLN.rds')

#####process mLN subT and divided into CD4/CD8####
dir.create('./4.4add.mLN.cd4.T/sub.T')
sub.T <- readRDS('./4.4add.mLN.cd4.T/mLn.sub.T.rds')
sub.T <- sub.T %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(sub.T,harmony = F)
pc = 15
sub.T <- FindNeighbors(sub.T,reduction = 'pca',dims = 1:pc,verbose = F)
set.resolutions <- seq(0.1,1,0.1)
sub.T <- FindClusters(sub.T,resolution = set.resolutions,verbose = F)
sub.T <- RunUMAP(sub.T,reduction = 'pca',dims = 1:pc,verbose = F)
pdf('./4.4add.mLN.cd4.T/sub.T/mLn.sub.T.res.observe.pdf')
clustree(sub.T)
DimPlot(sub.T,group.by = 'source',label = T,reduction = 'umap')
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.T, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
####select res 0.4
sub.T$seurat_clusters <- sub.T$RNA_snn_res.0.4
table(sub.T$seurat_clusters)
pdf('./4.4add.mLN.cd4.T/sub.T/mLn.sub.T.dimplot.pdf')
DimPlot(sub.T,group.by = 'seurat_clusters',label = T,reduction = 'umap')
dev.off()

pdf('./4.4add.mLN.cd4.T/sub.T/marker.observe.pdf')
FeaturePlot(sub.T,reduction = 'umap',features = c('CD4','CD8A','CD8B','MS4A1'),order = T)
VlnPlot(sub.T,group.by = 'seurat_clusters',features = c('CD4','CD8A','CD8B'))
GEB::DoHeatmap_my(sub.T,features = c('CD4','CD8A','CD8B'),group.by = 'seurat_clusters')
genelist <- c('CD4','CD40LG','CD8A','CD8B','CD3D','CD3E','CD3G','NKG7','KLRB1','TRGV9','TRDV2','TRDC')
exp.matrix <- GetAssayData(sub.T, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.T$seurat_clusters, mean)
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

plan("multisession",workers = 12)
cluster.pos.markers <- FindAllMarkers(sub.T, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1)
plan('multisession',workers = 1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.4add.mLN.cd4.T/sub.T/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.4add.mLN.cd4.T/sub.T/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.4add.mLN.cd4.T/sub.T/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

FeaturePlot(sub.T,reduction = 'umap',features = c('CD3D','CD3E','CD3G'),order = T,min.cutoff = 0)

####annotate
table(sub.T$seurat_clusters)
CD4.CD8 <- sub.T@meta.data$seurat_clusters
CD4.CD8 <- gsub("^0$", "CD8", CD4.CD8)
CD4.CD8 <- gsub("^5$", "CD8", CD4.CD8)
CD4.CD8 <- gsub("^8$", "NKT", CD4.CD8)
CD4.CD8 <- gsub("^10$", "CD8", CD4.CD8)
CD4.CD8 <- gsub("^11$", "NKT", CD4.CD8)
CD4.CD8 <- gsub('^[0-9]','CD4',CD4.CD8)

table(CD4.CD8)
sub.T <- AddMetaData(sub.T, CD4.CD8, col.name = "CD4.CD8")
table(sub.T$CD4.CD8)
pdf('./4.4add.mLN.cd4.T/sub.T/celltype.pdf')
DimPlot(sub.T,group.by = 'CD4.CD8',label = T,reduction = 'umap')
dev.off()

genelist <- c('CD4','CD40LG','CD8A','CD8B','CD3D','CD3E','CD3G','NKG7')
sub.T <- ScaleData(sub.T,features = c(VariableFeatures(sub.T),genelist))
pdf('./4.4add.mLN.cd4.T/sub.T/celltype.marker.observe.pdf')
GEB::DoHeatmap_my(sub.T,group.by = 'CD4.CD8',slot = 'scale.data',features = genelist)
dev.off()
saveRDS(sub.T,file = './4.4add.mLN.cd4.T/sub.T/annotated.mLn.sub.T.rds')

CD4 <- subset(sub.T,CD4.CD8=='CD4')
CD4 <- DietSeurat(CD4,assays = 'RNA')
saveRDS(CD4,file = './4.4add.mLN.cd4.T/sub.T/mLn.T.CD4.rds')

CD8 <- subset(sub.T,CD4.CD8=='CD8')
CD8 <- DietSeurat(CD8,assays = 'RNA')
saveRDS(CD8,file = './4.4add.mLN.cd4.T/sub.T/mLn.T.CD8.rds')

####merge CD4T of mLn and nLn#####
dir.create('./4.4add.mLN.cd4.T/sub.T/integrate.CD4')
mLn.CD4 <- readRDS('./4.4add.mLN.cd4.T/sub.T/mLn.T.CD4.rds')
nLn.CD4 <- readRDS('./4.3add.normal.cd4.T/CD4.merge.rds')
nLn.CD4 <- DietSeurat(nLn.CD4,assays = 'RNA')
mLn.CD4$paper <- mLn.CD4$source
table(mLn.CD4$paper)
table(nLn.CD4$paper)
CD4.merge <- merge(mLn.CD4,nLn.CD4)
table(CD4.merge$source)
CD4.merge$origin <- CD4.merge$source
CD4.merge$origin[grep(CD4.merge$origin,pattern = 'r[0-9]')] <- 'roider_nln'
table(CD4.merge$origin)
CD4.merge <- CD4.merge %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
####cca
split.object <- SplitObject(CD4.merge,split.by = 'origin')
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
DimPlot(CD4.merge,reduction = 'pca',split.by = 'origin',group.by = 'origin')
DimPlot(integrated.data,reduction = 'pca',split.by = 'origin',group.by = 'origin')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(integrated.data,harmony = F)
DefaultAssay(integrated.data) <- 'integrated'
pc <- 55
integrated.data <- FindNeighbors(integrated.data,dims = 1:pc,reduction = 'pca',verbose = F)
set.resolutions <- seq(0.1,1,0.1)
integrated.data <- FindClusters(integrated.data,resolution = set.resolutions,verbose = F)
integrated.data <- RunUMAP(integrated.data,reduction = 'pca',dims = 1:pc,verbose = F)
pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/cca.integrated.CD4.pdf')
clustree(integrated.data)
DimPlot(integrated.data,reduction = 'umap',split.by = 'origin',ncol = 3,group.by = 'origin')
DimPlot(integrated.data,group.by = 'origin',reduction = 'umap')
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = integrated.data, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x))
  print(p)
})
dev.off()

# ####harmony
# library(harmony)
# harmony.data <- RunHarmony(CD4.merge,group.by.vars = 'origin',verbose = F)
# DimPlot(CD4.merge,split.by = 'origin',reduction = 'pca',group.by = 'origin')
# DimPlot(harmony.data,split.by = 'origin',reduction = 'harmony',group.by = 'origin')
# DimPlot(harmony.data,reduction = 'harmony',group.by = 'origin')
# source('./function/do_Elbow_quantitative_modified_by_lhc.R')
# do_Elbow_quantitative(harmony.data,harmony = T)
# pc = 50
# harmony.data <- FindNeighbors(harmony.data,dims = 1:pc,reduction = 'harmony',verbose = F)
# set.resolutions <- seq(0.1,1.5,0.1)
# harmony.data <- FindClusters(harmony.data,resolution = set.resolutions,verbose = F)
# harmony.data <- RunUMAP(harmony.data,reduction = 'harmony',dims = 1:pc,verbose = F)
# pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/harmony.integrated.CD4.pdf')
# clustree(harmony.data)
# DimPlot(harmony.data,reduction = 'umap',split.by = 'origin',ncol = 2,group.by = 'origin')
# DimPlot(harmony.data,group.by = 'origin',reduction = 'umap')
# sce.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = harmony.data, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
#   print(p)
# })
# dev.off()
# 
# library(Nebulosa)
# plot_density(harmony.data,'MIR155HG')
# harmony.data$seurat_clusters <- harmony.data$RNA_snn_res.0.9
# table(harmony.data$seurat_clusters)
# DotPlot(harmony.data,features = c('MIR155HG','LAG3','C1orf228'),group.by = 'seurat_clusters')
# 
# harmony.data$origin <- harmony.data$source
# harmony.data$origin[grep(harmony.data$origin,pattern = 'r[0-9]')] <- 'roider_nln'
# table(harmony.data$origin)
# metadata <- harmony.data@meta.data
# 
# df1 <- data.frame(origin = rep('a',5),cluster = 0,CD4 = 0,percentage = 0)
# df1$origin <- names(table(metadata$origin))
# df1$CD4 <- as.numeric(table(metadata$origin))
# df1$cluster <- as.numeric(table(metadata$origin[metadata$seurat_clusters=='24']))
# df1$percentage <- round(df1$cluster/df1$CD4,3)*100




DefaultAssay(integrated.data) <- 'RNA'
integrated.data <- integrated.data %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
library(Nebulosa)
plot_density(integrated.data,'MIR155HG')
integrated.data$seurat_clusters <- integrated.data$integrated_snn_res.0.5
table(integrated.data$seurat_clusters)
DotPlot(integrated.data,features = c('MIR155HG','LAG3','C1orf228'),group.by = 'seurat_clusters',col.min = 0)
####select res 0.5
table(integrated.data$seurat_clusters)
integrated.data$seurat_clusters <- factor(integrated.data$seurat_clusters,levels = as.character(0:21),ordered = T)
Idents(integrated.data) <- integrated.data$seurat_clusters
table(Idents(integrated.data))
pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/CD4.dimplot.pdf')
DimPlot(integrated.data,group.by = 'seurat_clusters',cols = sample(Palettes[['mycols_28']],size = 22))
DimPlot(integrated.data,cells = WhichCells(integrated.data,idents = '11'),pt.size = .5,group.by = 'origin',cols = Palettes[['mycols_8']])
dev.off()


table(integrated.data$origin)
metadata <- integrated.data@meta.data

df1 <- data.frame(origin = rep('a',5),cluster = 0,CD4 = 0,percentage = 0)
df1$origin <- names(table(metadata$origin))
df1$CD4 <- as.numeric(table(metadata$origin))
df1$cluster <- as.numeric(table(metadata$origin[metadata$seurat_clusters=='11']))
df1$percentage <- round(df1$cluster/df1$CD4,3)*100

index <- which(!is.na(metadata$group))
metadata$group[index] <- paste('thiswork_',metadata$group[index])
metadata$group[-index] <- metadata$paper[-index]
table(metadata$group)
df2 <- data.frame(origin = rep('a',2),cluster = 0,CD4 = 0,percentage = 0)
metadata_2 <- metadata[metadata$source=='thiswork',]
df2$origin <- names(table(metadata_2$group))
df2$CD4 <- as.numeric(table(metadata_2$group))
df2$cluster <- as.numeric(table(metadata_2$group[metadata_2$seurat_clusters=='11']))
df2$percentage <- round(df2$cluster/df2$CD4,3)*100



df <- rbind(df1,df2)
pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/paper.texpre.percentage.pdf')
ggbarplot(df,x = 'origin',y = 'percentage',fill = 'origin',palette = 'npg',label = T,title = 'MIR155HG+ Tex-pre') + 
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(y = '% of CD4 T cells')
dev.off()

pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/marker.texpre.pdf')
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

pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/marker.texpre.dotplot.pdf',width = 10,height = 5)
DotPlot(integrated.data,features = c('MIR155HG','LAG3','TSHZ2','C1orf228'),group.by = 'seurat_clusters',scale.by = 'size',col.max = 2) + scale_color_gradientn(colors = Palettes[['blueRed']]) +coord_flip()
dev.off()


selected.cluster <- subset(integrated.data,seurat_clusters=='11')
selected.cluster$cluster <- '11'
source('./function/ratio_plot.R')
pdf('./4.4add.mLN.cd4.T/sub.T/integrate.CD4/texpre.ratio.pdf')
ratio.plot(seurat.object = selected.cluster, id.vars1 = "source", id.vars2 = "cluster", angle = 60,color.len = Palettes[['summerNight']])
dev.off()

table(Idents(integrated.data))
plan("multisession",workers = 11)
cluster.pos.markers <- FindAllMarkers(integrated.data, only.pos = TRUE, group.by = "seurat_clusters",test.use = 'MAST',latent.vars = 'source',
                                      min.diff.pct = 0.1,logfc.threshold = 0.1)
plan('multisession',workers = 1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.4add.mLN.cd4.T/sub.T/integrate.CD4/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.4add.mLN.cd4.T/sub.T/integrate.CD4/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.4add.mLN.cd4.T/sub.T/integrate.CD4/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

saveRDS(integrated.data,file = './4.4add.mLN.cd4.T/sub.T/integrate.CD4/integrated.CD4.rds')
