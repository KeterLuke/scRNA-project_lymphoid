####process and integrate external data(normal lymphoid tissue: 5 samples)
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
options(future.globals.maxSize = 128000 * 1024^2)
setwd('../')
getwd()
set.seed(123)
dir.create("./4.1process.external.data")
name=("./4.1process.external.data/")
source('./function/colorPalettes.R')

####1 GSE131907
##input normal reference ln data
library(data.table)
obj=readRDS("./raw/extra data/normal_lym/GSE131907_Lung_Cancer_raw_UMI_matrix.rds/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
ana=fread("./raw/extra data/normal_lym/GSE131907_Lung_Cancer_cell_annotation.txt/GSE131907_Lung_Cancer_cell_annotation.txt")
table(ana$Sample_Origin)
ana=ana[ana$Sample_Origin=="nLN",]
table(ana$Sample)
mat=obj[,ana$Index]

colnames(mat)=paste0(colnames(mat),"-",'GSE131907_nLN')
GSE131907_nLN=CreateSeuratObject(counts = mat, project = 'GSE131907_nLN', min.cells = 3, min.features = 200)
GSE131907_nLN$source="GSE131907_nLN"
GSE131907_nLN <- GSE131907_nLN %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(GSE131907_nLN,harmony = F)
GSE131907_nLN  <- FindNeighbors(object = GSE131907_nLN , dims = 1:11, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
GSE131907_nLN  <- FindClusters(object = GSE131907_nLN , resolution = set.resolutions, verbose = FALSE)
GSE131907_nLN  <- RunUMAP(GSE131907_nLN , dims = 1:11)
pdf(file = "4.1process.external.data/GSE131907_nLN_PCA-test.pdf")
clustree(GSE131907_nLN)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = GSE131907_nLN, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.1
saveRDS(GSE131907_nLN,file = '4.1process.external.data/GSE131907_nLN.rds')
####2 roider et al 2021
read_fun=function(name,dir){
  # name='r2_nln';dir="./input/validation_data/sc_dat/roider/rLN2/"
  library(Matrix)
  barcode.path <- paste0(dir, "barcodes.tsv")
  features.path <- paste0(dir, "features.tsv")
  matrix.path <- paste0(dir, "matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)

  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  colnames(mat)=paste0(colnames(mat),"-",name)
  obj=CreateSeuratObject(counts = mat, project = name, min.cells = 3, min.features = 200)
  obj$source=name
  return(obj)
}

r1_nln=read_fun('r1_nln',"./raw/extra data/roider/rLN1/")
r2_nln=read_fun('r2_nln',"./raw/extra data/roider/rLN2/")
r3_nln=read_fun('r3_nln',"./raw/extra data/roider/rLN3/")
####r1_nln
r1_nln <- r1_nln %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(r1_nln,harmony = F)
r1_nln  <- FindNeighbors(object = r1_nln , dims = 1:11, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
r1_nln  <- FindClusters(object = r1_nln , resolution = set.resolutions, verbose = FALSE)
r1_nln  <- RunUMAP(r1_nln , dims = 1:11)
pdf(file = "4.1process.external.data/r1_nln_PCA-test.pdf")
clustree(r1_nln)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = r1_nln, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.4
saveRDS(r1_nln,file = '4.1process.external.data/r1_nln.rds')
####r2_nln
r2_nln <- r2_nln %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(r2_nln,harmony = F)
r2_nln  <- FindNeighbors(object = r2_nln , dims = 1:9, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
r2_nln  <- FindClusters(object = r2_nln , resolution = set.resolutions, verbose = FALSE)
r2_nln  <- RunUMAP(r2_nln , dims = 1:9)
pdf(file = "4.1process.external.data/r2_nln_PCA-test.pdf")
clustree(r2_nln)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = r2_nln, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.1
saveRDS(r2_nln,file = '4.1process.external.data/r2_nln.rds')
####r3_nln
r3_nln <- r3_nln %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(r3_nln,harmony = F)
r3_nln  <- FindNeighbors(object = r3_nln , dims = 1:13, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
r3_nln  <- FindClusters(object = r3_nln , resolution = set.resolutions, verbose = FALSE)
r3_nln  <- RunUMAP(r3_nln , dims = 1:13)
pdf(file = "4.1process.external.data/r3_nln_PCA-test.pdf")
clustree(r3_nln)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = r3_nln, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.1
saveRDS(r3_nln,file = '4.1process.external.data/r3_nln.rds')
###3 gse182436
mat=read.table("./raw/extra data/1/GSE182434_raw_count_matrix.txt/GSE182434_raw_count_matrix.txt",header = T,row.names = 1)
inf=fread("./raw/extra data/1/GSE182434_cell_annotation.txt/GSE182434_cell_annotation.txt")
table(inf$TumorNormal)
cells=inf$ID[inf$TumorNormal=='Normal']
mat=mat[,cells]
colnames(mat)=paste0(colnames(mat),"-",'gse182436_n')
obj_gse182436=CreateSeuratObject(counts = mat, project = 'gse182436_n', min.cells = 3, min.features = 200)
obj_gse182436$source='obj_gse182436_nln'
table(obj_gse182436$source)
obj_gse182436 <- obj_gse182436 %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(obj_gse182436,harmony = F)
obj_gse182436  <- FindNeighbors(object = obj_gse182436 , dims = 1:12, verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
obj_gse182436  <- FindClusters(object = obj_gse182436 , resolution = set.resolutions, verbose = FALSE)
obj_gse182436  <- RunUMAP(obj_gse182436 , dims = 1:12)
pdf(file = "4.1process.external.data/obj_gse182436_PCA-test.pdf")
clustree(obj_gse182436)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = obj_gse182436, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
###res:0.2
saveRDS(obj_gse182436,file = '4.1process.external.data/obj_gse182436.rds')
####remove doublet####
rm(list = ls())
source('./function/doubletDetect.R')
dir.create('./4.1process.external.data/doublet')
####GSE131907_nLN
sce <- readRDS('./4.1process.external.data/GSE131907_nLN.rds')
pdf("4.1process.external.data/doublet/GSE131907_nLN_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:11, doublet.rate = 0.1, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.1
saveRDS(sce,file = '4.1process.external.data/doublet/GSE131907_nLN.RDS')
####r1_nln
sce <- readRDS('./4.1process.external.data/r1_nln.rds')
pdf("4.1process.external.data/doublet/r1_nln_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:11, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.4", sct = F)
dev.off()
##select:rate = 0.02476
saveRDS(sce,file = '4.1process.external.data/doublet/r1_nln.RDS')
####r2_nln
sce <- readRDS('./4.1process.external.data/r2_nln.rds')
pdf("4.1process.external.data/doublet/r2_nln_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:9, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.016936
saveRDS(sce,file = '4.1process.external.data/doublet/r2_nln.RDS')
####r3_nln
sce <- readRDS('./4.1process.external.data/r3_nln.rds')
pdf("4.1process.external.data/doublet/r3_nln_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:13, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.022792
saveRDS(sce,file = '4.1process.external.data/doublet/r3_nln.RDS')
####obj_gse182436
sce <- readRDS('./4.1process.external.data/obj_gse182436.rds')
pdf("4.1process.external.data/doublet/obj_gse182436_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:12, doublet.rate = 0.01, annotation = "RNA_snn_res.0.2", sct = F)
dev.off()
##select:rate = 0.01
saveRDS(sce,file = '4.1process.external.data/doublet/obj_gse182436.RDS')

####merge external data
GSE131907_nLN <- readRDS('./4.1process.external.data/doublet/GSE131907_nLN.RDS')
GSE131907_nLN <- subset(GSE131907_nLN,Doublet=='Singlet')
r1_nln <- readRDS('./4.1process.external.data/doublet/r1_nln.RDS')
r1_nln <- subset(r1_nln,Doublet=='Singlet')
r2_nln <- readRDS('./4.1process.external.data/doublet/r2_nln.RDS')
r2_nln <- subset(r2_nln,Doublet=='Singlet')
r3_nln <- readRDS('./4.1process.external.data/doublet/r3_nln.RDS')
r3_nln <- subset(r3_nln,Doublet=='Singlet')
obj_gse182436 <- readRDS('./4.1process.external.data/doublet/obj_gse182436.RDS')
obj_gse182436 <- subset(obj_gse182436,Doublet=='Singlet')

files=c('GSE131907_nLN','r1_nln','r2_nln','r3_nln',
        'obj_gse182436')
i=1
if(i==1){
  print(i)
  name1=files[1]
  name2=files[2]
  merged_seurat=merge(get(name1),y=get(name2))
} 

for(j in 3:length(files)){
  print(j)
  name3=files[j]
  merged_seurat=merge( merged_seurat,y=get(name3)) }
table(merged_seurat$source)
table(merged_seurat$Doublet)
merged_seurat <- merged_seurat%>%NormalizeData()%>%
  FindVariableFeatures(nfeatures = 3000)%>%
  ScaleData()%>%
  RunPCA(npcs = 100)
####harmonny
library(harmony)
harmony_seurat <- RunHarmony(merged_seurat,group.by.vars = 'source',verbose = F)
DimPlot(harmony_seurat,reduction = 'pca',split.by = 'source')
DimPlot(harmony_seurat,reduction = 'harmony',split.by = 'source')
DimPlot(harmony_seurat,reduction = 'harmony',group.by = 'source')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(harmony_seurat,harmony = T)
harmony_seurat <- FindNeighbors(harmony_seurat,dims = 1:35,reduction = 'harmony',verbose = F)
set.resolutions <- seq(0.1,1,0.1)
harmony_seurat <- FindClusters(harmony_seurat,resolution = set.resolutions,verbose = F)
harmony_seurat <- RunUMAP(harmony_seurat,reduction = 'harmony',dims = 1:35,verbose = F)
pdf(file = "4.1process.external.data/externaldata.integrated.harmony.pdf",width = 8)
clustree(harmony_seurat)
DimPlot(harmony_seurat,reduction = 'umap',split.by = 'source',ncol = 3)
DimPlot(harmony_seurat,reduction = 'umap',group.by  = 'source')

sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = harmony_seurat, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
saveRDS(harmony_seurat,file = './4.1process.external.data/harmony.integrated.external.data.rds')
####cca
split_seurat <- SplitObject(merged_seurat, split.by = "source")
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
intregrated_objects=pca_cca_pipline(split_seurat)
intregrated_objects <- ScaleData(intregrated_objects)
intregrated_objects <- RunPCA(intregrated_objects,npcs = 100)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(intregrated_objects,harmony = F)
DimPlot(intregrated_objects,reduction = 'pca',split.by = 'source')
DimPlot(intregrated_objects,reduction = 'pca',group.by = 'source')
intregrated_objects  <- FindNeighbors(object = intregrated_objects , dims = 1:20,reduction = 'pca', verbose = FALSE)
set.resolutions <- seq(0.1,1,0.1)
intregrated_objects  <- FindClusters(object = intregrated_objects , resolution = set.resolutions,verbose = FALSE)
intregrated_objects  <- RunUMAP(intregrated_objects , reduction = 'pca',dims = 1:20)
pdf(file = "4.1process.external.data/externaldata.integrated.cca.pdf",width = 8)
clustree(intregrated_objects)
DimPlot(intregrated_objects,reduction = 'umap',split.by = 'source',ncol = 3)
DimPlot(intregrated_objects,reduction = 'umap',group.by = 'source')
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = intregrated_objects, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x))
  print(p)
})
dev.off()

DefaultAssay(intregrated_objects) <- 'RNA'
intregrated_objects <- intregrated_objects %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
FeaturePlot(intregrated_objects,features = c('CD3D','CD3E','MS4A1','CD79A'),order = T,min.cutoff = 1)
DimPlot(intregrated_objects,group.by = 'source',split.by = 'source',ncol = 3)
DimPlot(intregrated_objects,group.by = 'source')
FeaturePlot(intregrated_objects,features = c('CD8A','CD8B'),order = T,min.cutoff = 1)
saveRDS(intregrated_objects,file = './4.1process.external.data/cca.integrated.external.data.rds')

####annotate celltype
####select: cca integrated results
####res:0.3
external.data <- readRDS('./4.1process.external.data/cca.integrated.external.data.rds')
external.data$seurat_clusters <- external.data$integrated_snn_res.0.3
pdf('./4.1process.external.data/externaldata.dimplot.pdf')
DimPlot(external.data,group.by = 'seurat_clusters',reduction = 'umap',label = T)
DimPlot(external.data,group.by = 'source',reduction = 'umap')
DimPlot(external.data,group.by = 'source',split.by = 'source',reduction = 'umap',ncol = 3)
dev.off()
FeaturePlot(external.data,features = c('EPCAM','PTPRC','CD3D','MS4A1'),order = T,min.cutoff = 0)
FeaturePlot(external.data,features = c('MKI67','TOP2A'),order = T)
Idents(external.data) <- external.data$seurat_clusters
table(Idents(external.data))
plan("multisession",workers = 12)
cluster.pos.markers <- FindAllMarkers(external.data, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
plan('multisession',workers = 1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "4.1process.external.data/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "4.1process.external.data/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "4.1process.external.data/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#####expression of classical markers
cell.type.markers <- read.table(file = "2.Cluster/Annotate/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(external.data, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, external.data$seurat_clusters, mean)
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
pdf("4.1process.external.data/cluster.signature.expression.pdf")
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
a <- external.data
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
####annotate
table(external.data$seurat_clusters)
major_celltype <- external.data@meta.data$seurat_clusters
major_celltype <- gsub("^0$", "T", major_celltype)
major_celltype <- gsub("^1$", "B", major_celltype)
major_celltype <- gsub("^2$", "B", major_celltype)
major_celltype <- gsub("^3$", "T", major_celltype)
major_celltype <- gsub("^4$", "T", major_celltype)
major_celltype <- gsub("^5$", "T", major_celltype)
major_celltype <- gsub("^6$", "Plasma", major_celltype)
major_celltype <- gsub("^7$", "NK", major_celltype) 
major_celltype <- gsub("^8$", "B", major_celltype)
major_celltype <- gsub("^9$", "Cycling", major_celltype)
major_celltype <- gsub("^10$", "T", major_celltype)
major_celltype <- gsub("^11$", "Myeloid", major_celltype)
major_celltype <- gsub("^12$", "B", major_celltype)
major_celltype <- gsub("^13$", "pDC", major_celltype)
# major_celltype <- gsub("^14$", "pDC", major_celltype)
# major_celltype <- gsub("^15$", "T", major_celltype)
# major_celltype <- gsub("^16$", "LAMP3+DC", major_celltype)


table(major_celltype)
external.data <- AddMetaData(external.data, major_celltype, col.name = "major_celltype")
table(external.data$major_celltype)

pdf('./4.1process.external.data/celltype.pdf')
DimPlot(external.data,group.by = 'major_celltype',label = T,reduction = 'umap')
DimPlot(external.data,group.by = 'major_celltype',split.by = 'source',reduction = 'umap',ncol = 3)
dev.off()

pdf('./4.1process.external.data/observe.marker.pdf')
FeaturePlot(external.data,features = c('CD3D','CD3E','MS4A1','CD79A'),order = T,min.cutoff = 1)
VlnPlot(external.data,features = c('CD3D','CD3E','MS4A1','CD79A'),group.by = 'seurat_clusters',ncol = 2,pt.size = 0)
dev.off()

genelist=c('MS4A1','CD19','CD79A','CD79B','MKI67','TOP2A','STMN1','TUBB',
           'APOC1','LYZ','CST3','APOE','NKG7','KLRB1',
           'GNLY','TRDC','LILRA4','IRF7','FCER1G','TCF4','SDC1','MZB1','IGHG1','IGHG3',
           'CD3D','CD3E','CD3G','TRBC1')
external.data$major_celltype <- factor(external.data$major_celltype,levels = c('B','Cycling','Myeloid','NK','pDC','Plasma','T'))
pdf("4.1process.external.data/cellType.marker.dotplot.pdf", height = 5,width = 12)
Idents(external.data) <- external.data$major_celltype
DotPlot(external.data,features = genelist,scale.by = 'size',col.min = 0) + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 8,angle = 25))
dev.off()

Idents(external.data) <- external.data$major_celltype
exp.matrix <- GetAssayData(external.data, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, external.data$major_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('4.1process.external.data/cellType.marker.heatmap.pdf')
Heatmap(t(cluster.score.normailzed), 
       column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =15),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()
pdf('./4.1process.external.data/celltype.marker.observe.pdf')
VlnPlot(external.data,features = c('CD3D','CD3E','MS4A1','CD79A'),pt.size = 0,ncol = 2)
dev.off()
saveRDS(external.data,file = './4.1process.external.data/externaldata.annotated.rds')

table(external.data$major_celltype)
T_sub <- subset(external.data,subset = major_celltype=="T")
FeaturePlot(T_sub,features = c('CD3D','CD3E','CD3G','CD7'),order = T,min.cutoff = 0)
T_sub <- DietSeurat(T_sub,assays = 'RNA')
saveRDS(T_sub,file = './4.1process.external.data/external.T.sub.rds')
dir.create('./4.2process.external.T.sub/')
