

#' @description: annotate T subclusters

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
sub.scRNA <- subset(data.merge, subset = major_celltype == 'T')
rm(data.merge)
DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")####remove SCT assay, scale.data and dimension reduction 
# observe the batch effect
# sub.scRNA <- SCTransform(sub.scRNA, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
# sub.scRNA <- RunPCA(sub.scRNA, npcs = 30, verbose = FALSE) %>%
# RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
# FindClusters(resolution = seq(0.2, 1.2, by = 0.1), verbose = FALSE)

# split the dataset into a list of seurat objects by samples####
sub.list <- SplitObject(sub.scRNA, split.by = "orig.ident")
dir.create('./3.1Subcluster_T')
sub.list.Standard <- variableFeatureSelection(seurat.lists = sub.list, method = "Standard", nfeatures = 3000)
saveRDS(sub.list.Standard, file = "3.1Subcluster_T/sub.list.Standard.3000.rds")

sub.list.SCT <- variableFeatureSelection(seurat.lists = sub.list, method = "SCT", nfeatures = 3000,return.only.var.genes = T,vars.to.regress = c("nCount_RNA", "percent.mt"))
saveRDS(sub.list.SCT, file = "3.1Subcluster_T/sub.list.SCT.3000.rds")

## Remove previous clustering results
sub.scRNA <- merge(sub.list.SCT[[1]], y = sub.list.SCT[2:length(sub.list.SCT)], project = "T cells")
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(sub.scRNA@meta.data))
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-index]
View(sub.scRNA@meta.data)

# assay=SCT
DefaultAssay(sub.scRNA) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = sub.list.SCT, nfeatures = 3000)
VariableFeatures(sub.scRNA) <- seurat.features.SCT
# assay=RNA
DefaultAssay(sub.scRNA) <- "RNA"
seurat.features.RNA <- SelectIntegrationFeatures(object.list = sub.list.Standard, nfeatures = 3000)
VariableFeatures(sub.scRNA) <- seurat.features.RNA
sub.scRNA <- NormalizeData(sub.scRNA, verbose = FALSE)
sub.scRNA <- ScaleData(sub.scRNA, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(sub.scRNA))
####SCT.Harmony.PC20####
pdf("3.1Subcluster_T/SCT.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()

####SCT.Harmony.PC30####
pdf("3.1Subcluster_T/SCT.Harmony.PC30.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()

####SCT.Harmony.PC40####
pdf("3.1Subcluster_T/SCT.Harmony.PC40.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 40, npcs = 100)
dev.off()

####SCT.Harmony.PC45####
pdf("3.1Subcluster_T/SCT.Harmony.PC45.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 45, npcs = 100)
dev.off()

####SCT.Harmony.PC50####
pdf("3.1Subcluster_T/SCT.Harmony.PC50.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 50, npcs = 100)
dev.off()


####RNA.Harmony.PC7####
pdf("3.1Subcluster_T/RNA.Harmony.PC7.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 7, npcs = 100)
dev.off()

####RNA.Harmony.PC20####
pdf("3.1Subcluster_T/RNA.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()

####RNA.Harmony.PC25####
pdf("3.1Subcluster_T/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()

####RNA.Harmony.PC27####
pdf("3.1Subcluster_T/RNA.Harmony.PC27.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 27, npcs = 100)
dev.off()

####RNA.Harmony.PC30####
pdf("3.1Subcluster_T/RNA.Harmony.PC30.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()

####RNA.Harmony.PC40####
pdf("3.1Subcluster_T/RNA.Harmony.PC40.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 40, npcs = 100)
dev.off()

####RNA.Harmony.PC45####
pdf("3.1Subcluster_T/RNA.Harmony.PC45.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 45, npcs = 100)
dev.off()

####RNA.Harmony.PC50####
pdf("3.1Subcluster_T/RNA.Harmony.PC50.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 50, npcs = 100)
dev.off()

##select RNA.Harmony.PC25
pdf("3.1Subcluster_T/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()

sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$RNA_snn_res.0.7
table(sub.scRNA.harmony$seurat_clusters)
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
DefaultAssay(sub.scRNA.harmony) <- "RNA"
##remove mixed clusters####
sub.scRNA.harmony <- subset(sub.scRNA.harmony,subset = seurat_clusters %in% c(0:9,11:12,14))
sub.scRNA.harmony <- NormalizeData(sub.scRNA.harmony, verbose = FALSE)
sub.scRNA.harmony <- FindVariableFeatures(sub.scRNA.harmony,nfeatures = 3000)
sub.scRNA.harmony <- ScaleData(sub.scRNA.harmony, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = row.names(sub.scRNA.harmony@assay$RNA@data))
sub.scRNA.harmony <- RunPCA(sub.scRNA.harmony,verbose = F,npcs = 100)
sub.scRNA.harmony <- RunUMAP(sub.scRNA.harmony,verbose = F,reduction = 'pca',dims = 1:20)
DimPlot(sub.scRNA.harmony,reduction = 'umap',group.by = 'orig.ident')
options(timeout = 1600000000)
####perform integration and reduction####
sub.scRNA.harmony <- RunHarmony(sub.scRNA.harmony,reduction = 'pca',verbose = F,group.by.vars = 'orig.ident')
sub.scRNA.harmony <- RunUMAP(sub.scRNA.harmony,verbose = F,reduction = 'harmony',dims = 1:20)
sub.scRNA.harmony <- RunTSNE(sub.scRNA.harmony,verbose = F,reduction = 'harmony',dims = 1:20)
DimPlot(sub.scRNA.harmony,reduction = 'umap',group.by = 'orig.ident')
sub.scRNA.harmony <- FindNeighbors(sub.scRNA.harmony,verbose = F,reduction = 'harmony',dims = 1:20)
sub.scRNA.harmony <- FindClusters(sub.scRNA.harmony,resolution = seq(0.1, 1, by = 0.1),verbose = F)
pdf('./3.1Subcluster_T/RNA.Harmony.PC25.pro.PC20.pdf')
set.resolutions <- seq(0.1, 1, by = 0.1)
clustree(sub.scRNA.harmony)
DimPlot(sub.scRNA.harmony,reduction = 'umap',group.by = 'orig.ident')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.scRNA.harmony, reduction = 'umap',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
DimPlot(sub.scRNA.harmony,reduction = 'tsne',group.by = 'orig.ident')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.scRNA.harmony, reduction = 'tsne',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
dev.off()
####select resolution:0.5
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$RNA_snn_res.0.5
table(sub.scRNA.harmony$seurat_clusters)
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
table(Idents(sub.scRNA.harmony))
table(sub.scRNA.harmony$celltypist_minor,sub.scRNA.harmony$seurat_clusters)
saveRDS(sub.scRNA.harmony, file = "3.1Subcluster_T/sub.scRNA.harmony.rds")

sub.T <- readRDS('./3.1Subcluster_T/sub.scRNA.harmony.rds')
length(unique(sub.T$seurat_clusters))
dir.create('./3.1Subcluster_T/Annotate')
pdf('./3.1Subcluster_T/Annotate/cluster.pdf')
DimPlot(sub.T,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
DimPlot(sub.T,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
dev.off()
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.rds')
#### classical marker#### 
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
pdf('./3.1Subcluster_T/Annotate/classical.marker.cd4.pdf')
plot_density(sub.T,features  = c('CD4','CCR7','LEF1','TCF7'),reduction = 'umap',joint = F)##Naive
plot_density(sub.T,features  = c('CD4','CCR7','IL2','SELL'),reduction = 'umap',joint = F)##Naive/Tcm
plot_density(sub.T,features  = c('CD4','CCR7','ANXA1','ANXA2'),reduction = 'umap',joint = F)##Tcm
plot_density(sub.T,features  = c('CD4','CD40LG','ANXA1','S100A4'),reduction = 'umap',joint = F)##Tem
plot_density(sub.T,features  = c('CD4','CTLA4','FOXP3','IL2RA'),reduction = 'umap',joint = F)##Treg
plot_density(sub.T,features  = c('CD4','PDCD1','TIGIT','LAG3'),reduction = 'umap',joint = F)##Tex
dev.off()
pdf('./3.1Subcluster_T/Annotate/classical.marker.cd8.pdf')
plot_density(sub.T,features  = c('CD8A','CD8B','CCR7','LEF1'),reduction = 'umap',joint = F)##Naive
plot_density(sub.T,features  = c('CD8A','CD8B','IFNG','GNLY'),reduction = 'umap',joint = F)##CTL
dev.off()
##cluster similarity##
expMatrix <- GetAssayData(sub.T, slot = "scale.data")
highVariableGenes <- VariableFeatures(sub.T)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, sub.T$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.1Subcluster_T/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
cor <- (cor - min(cor))/(max(cor)-min(cor))
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'ward.D2')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.1Subcluster_T/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

##expression of functional markers
cell.type.markers <- read.table(file = "3.1Subcluster_T/Annotate/T_functional markers.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(sub.T, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.T$seurat_clusters, mean)
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
pdf("3.1Subcluster_T/Annotate/cluster.functional.signature.pdf",height = 10)
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

a <- sub.T
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
cell.type.markers_distinct <- cell.type.markers %>% distinct(Gene,.keep_all = T)
gene_list <- list(Naive = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Type=='Naive'],
                  Cytotoxic_effector = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Type=='Cytotoxic and effector'],
                  Memory = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Type=='Memory'],
                  Treg = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Type=='Treg'],
                  Exhausted = cell.type.markers_distinct$Gene[cell.type.markers_distinct$Type=='Exhausted'])
library(GEB)
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "seurat_clusters",dot.scale = 3)
dev.off()
#### Differential expression####
sub.T$seurat_clusters <- sub.T$RNA_snn_res.0.5
Idents(sub.T) <- sub.T$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.T, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sub.T <- ScaleData(sub.T,features = c(VariableFeatures(sub.T),unique(top5.genes$gene)))
pdf("3.1Subcluster_T/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(sub.T, features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['mycols_12']],group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

top10.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
sub.T <- ScaleData(sub.T,features = c(VariableFeatures(sub.T),unique(top10.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('./3.1Subcluster_T/Annotate/cluster.top10genes.mean.pdf',height = 12)
DoHeatmap_my(sub.T,features = unique(top10.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
####CD4 CED8
library(reticulate)
py_config()
scanpy <- import('scanpy')
pandas <- import('pandas')
numpy <- import('numpy')
plt <- import('matplotlib.pyplot')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(sub.T[['RNA']]@counts)))),
                       obs = pandas$DataFrame(sub.T@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(sub.T[['RNA']]@counts),
                                                         row.names = rownames(sub.T[['RNA']]@counts)))
)
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
scanpy$plotting$violin(adata,groupby = 'seurat_clusters',scale = 'width',keys = c('CD4'),stripplot = F,show = F)
plt$savefig('./3.1Subcluster_T/Annotate/CD4.pdf', bbox_inches = 'tight')

scanpy$plotting$violin(adata,groupby = 'seurat_clusters',scale = 'width',keys = c('CD8A','CD8B'),stripplot = F,show = F)
plt$savefig('./3.1Subcluster_T/Annotate/CD8.pdf', bbox_inches = 'tight')

pdf('3.1Subcluster_T/Annotate/CD4_CD8.pdf')
FeaturePlot(sub.T,reduction = 'umap',features = c('CD4','CD8A','CD8B'),order = T,cols = Palettes[['greyMagma']],min.cutoff = .5)
FeaturePlot(sub.T,reduction = 'tsne',features = c('CD4','CD8A','CD8B'),order = T,cols = Palettes[['greyMagma']],min.cutoff = .5)
dev.off()
####celltypist predict####
celltypist <- import('celltypist')
predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- sub.T@meta.data$seurat_clusters
sub.T <- AddMetaData(sub.T,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(sub.T$celltypist_minor,sub.T$seurat_clusters)

write.xlsx(predicted_minor,file = '3.1Subcluster_T/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(sub.T,file = '3.1Subcluster_T/sub.T.rds')

####enrichment####
cluster.DE <- FindAllMarkers(sub.T,only.pos = F,group.by = 'seurat_clusters',min.diff.pct = 0.1,test.use = 'MAST',
                             latent.vars = 'orig.ident')
idents <- levels(sub.T)
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "3.1Subcluster_T/Annotate/cluster.all.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "3.1Subcluster_T/Annotate/cluster.all.DE.rds")

##-- Functional enrichment analysis
DEGs <- lapply(saveFormat, function(x){
  x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.4) %>% arrange(desc(avg_log2FC))
  return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "3.1Subcluster_T/Annotate/clusters.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "3.1Subcluster_T/Annotate/clusters.DEGs.rds")
####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
DEGs <- readRDS('./3.1Subcluster_T/Annotate/clusters.DEGs.rds')
idents <- levels(sub.T)
DEGs <- lapply(DEGs, function(x){
  return(x$gene)
})
source("function/clusterProfiler.enricher.R")
####enricher
pdf('./3.1Subcluster_T/Annotate/enrichment/msigDB/cluster.clusterProfiler.enricher.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB",saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/enrichment/msigDB/'),
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
write.xlsx(res, file = "3.1Subcluster_T/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.rds')
####enrichgo
pdf('./3.1Subcluster_T/Annotate/enrichment/GO/cluster.clusterProfiler.enrichgo.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "GO",GO.ont = 'BP',filename = 'enrichgoResult.csv',
                                 saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/enrichment/GO/'),
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
write.xlsx(res, file = "3.1Subcluster_T/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.rds')
#### show cluster profiler result####
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
res <- lapply(1:12,function(i){
  pro <- read.xlsx('./3.1Subcluster_T/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.selected.xlsx',sheet = i)
  return(pro)
})
names(res) <- 0:11
library(ggthemes)
enrich.pathways <- unlist(lapply(res,function(i){
  id <- i$Description
  return(id)
}))
enrich.pathways <- unique(enrich.pathways)       
enrich <- lapply(res, function(x){
  idx <- which(x$Description %in% enrich.pathways)
  return(x[idx, c("Description", "p.adjust", "Count")])
})                    
names(enrich) <- names(res)
enrich <- lapply(names(enrich), function(x){
  enrich[[x]]$Cluster <- rep(x, nrow(enrich[[x]]))
  return(enrich[[x]])
})  
enrich.res <- Reduce(function(x,y) rbind(x,y), enrich)
rownames(enrich.res) <- NULL

enrich.res$p.adjust <- -log10(enrich.res$p.adjust)
enrich.res$wrap <- wrapText(enrich.res$Description, 45)
enrich.res$Cluster <- factor(enrich.res$Cluster,levels = c(0:11),ordered = T)
enrich.res$p.adjust[enrich.res$p.adjust >= 12] <- 12
pdf("3.1Subcluster_T/Annotate/enrichment/GO/cluterProfiler_go.enrichment.pdf",height = 15)
p1 <- ggplot(enrich.res, 
             aes(x = Cluster, 
                 y = wrap, 
                 size = Count, 
                 color = p.adjust,
                 fill = p.adjust,))
p2 <- p1 + guides(color='none') + geom_point(shape = 21) + theme_few() + scale_color_gradient(low = "#FF9375", high = "red") +
  scale_fill_gradient(low = "#FF9375", high = "red",breaks=c(2,5,8,12),limits = c(1.3,14))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 5)) + xlab("") + ylab("") + labs(fill='-log10(FDR)')
print(p2)

#heatmap
pathway.p <- enrich.res[, c("wrap", "Cluster", "p.adjust")]
pathway.data <- pivot_wider(pathway.p, names_from = "Cluster", values_from = "p.adjust")
pathway.data <- as.data.frame(pathway.data)
rownames(pathway.data) <- pathway.data$wrap
pathway.data <- as.data.frame(pathway.data[,-1])
cols <- colorRamp2(c(3,10), c("#ffff66", "red"))
Heatmap(pathway.data, name = "-log10(FDR)", show_column_dend = F, col = cols,show_row_dend = F,  
        na_col = "#f2f3f4", border = "grey", border_gp = gpar(col = "grey"), rect_gp = gpar(col = "grey"),
        width = unit(5, "cm"), height = unit(24, "cm"), cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 10))
dev.off()

####GSVA####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
Idents(sub.T) <- sub.T$seurat_clusters
av <- AverageExpression(sub.T,assays = 'RNA',group.by = 'seurat_clusters',slot = 'counts')
av <- av[[1]]
head(av)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(GSVA)
####gsvago
# gene_sets <- msigdbr(species = 'Homo sapiens',category = 'C5')
imm_sets <- data.table::fread('./3.1Subcluster_T/Annotate/enrichment/immune_related.txt')
# rn <- gene_sets$gs_name
# rn <- unlist(lapply(rn,function(i){
#   i <- gsub('GOBP_',"",i)
#   i <- gsub("_"," ",i)
#   return(i)
# }))
# rn <- tolower(rn)
# gene_sets$gs_name <- rn
# View(gene_sets)
# ####set the results of goenrichment as background
# go <- readRDS('./3.1Subcluster_T/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.rds')
# go_terms <- unique(unlist(lapply(go,function(i){
#  go <- i$Description
#   return(go)
# })))
# gene_sets <- gene_sets[gene_sets$gs_name%in%go_terms,]

gs=split(imm_sets$`Gene name`,imm_sets$`GO term name`)
gs = lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds,goId)
  {GeneSet(geneIds,geneIdType=EntrezIdentifier(),collectionType=GOCollection(goId),setName=goId)
},gs, names(gs)))
gsc
geneset <- gsc
X <- av
library(GSVA)
es.max <- gsva(X,geneset,
                mx.diff=F, verbose=T,parallel.sz=11)
es <- es.max
df <- do.call(rbind,lapply(1:ncol(es),function(i){
  dat <- data.frame(
    go = rownames(es),
    cluster = colnames(es)[i],
    sd.1 = es[,i],
    sd.2 = apply(es[,-i],1,median)
  )
}))
df$fc <- df$sd.1-df$sd.2
top <- df%>% group_by(cluster) %>% top_n(5,fc)
pdf('./3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaGo.pdf',width = 8)
pheatmap(es[unique(top$go),],show_rownames = T,clustering_method = 'ward.D2',fontsize = 5,fontsize_col = 8,angle_col = '45',main = 'GSVA_clusters')
dev.off()
saveRDS(es.max,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaGO.rds')
idents <- levels(sub.T)
df$cluster <- factor(df$cluster,levels = levels(sub.T))
res <- split(df,df$cluster)
saveFormat <- lapply(res,function(i){
  i <- i %>% dplyr::filter(fc > 0) %>% dplyr::arrange(desc(fc))
  return(i)
})
write.xlsx(saveFormat,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaGo.xlsx',sheetName = idents)
####gsvakegg
gene_sets <- msigdbr(species = 'Homo sapiens',category = 'C2',subcategory = 'KEGG')
rn <- gene_sets$gs_name
rn <- unlist(lapply(rn,function(i){
  i <- gsub('KEGG_',"",i)
  i <- gsub("_"," ",i)
  return(i)
}))
rn <- tolower(rn)
gene_sets$gs_name <- rn
View(gene_sets)
gs=split(gene_sets$gene_symbol,gene_sets$gs_name)
gs = lapply(gs, unique)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
gsc <- GeneSetCollection(mapply(function(geneIds,keggId)
{GeneSet(geneIds,geneIdType=EntrezIdentifier(),collectionType=KEGGCollection(keggId),setName=keggId)
},gs, names(gs)))
gsc
geneset <- gsc
X <- av
es.max <- gsva(X,geneset,
               mx.diff=F,verbose=T,parallel.sz=11)
es <- es.max
df <- do.call(rbind,lapply(1:ncol(es),function(i){
  dat <- data.frame(
    kegg = rownames(es),
    cluster = colnames(es)[i],
    sd.1 = es[,i],
    sd.2 = apply(es[,-i],1,median)
  )
}))
df$fc <- df$sd.1 - df$sd.2
top <- df%>% group_by(cluster) %>% top_n(8,fc)
pdf('./3.1Subcluster_T/Annotate/enrichment/GSVA/gsvakegg.pdf')
pheatmap(es[unique(top$kegg),],show_rownames = T,fontsize = 6,fontsize_col = 8,clustering_method = 'ward.D',angle_col = '45',main = 'GSVA_clusters')
dev.off()
saveRDS(es,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvakegg.rds')
idents <- levels(sub.T)
df$cluster <- factor(df$cluster,levels = levels(sub.T))
res <- split(df,df$cluster)
saveFormat <- lapply(res,function(i){
  i <- i %>% dplyr::filter(fc > 0) %>% dplyr::arrange(desc(fc))
  return(i)
})
write.xlsx(saveFormat,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvakegg.xlsx',sheetName = idents)
####gsvareactome
gene_sets <- msigdbr(species = 'Homo sapiens',category = 'C2',subcategory = 'REACTOME')
rn <- gene_sets$gs_name
rn <- unlist(lapply(rn,function(i){
  i <- gsub('REACTOME_',"",i)
  i <- gsub("_"," ",i)
  return(i)
}))
rn <- tolower(rn)
gene_sets$gs_name <- rn
View(gene_sets)
gs=split(gene_sets$gene_symbol,gene_sets$gs_name)
gs = lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds,roId)
{GeneSet(geneIds,geneIdType=EntrezIdentifier(),collectionType=GOCollection(roId,evidenceCode = 'IC'),setName=roId)
},gs, names(gs)))
gsc
geneset <- gsc
X <- av
es.max <- gsva(X,geneset,
               mx.diff=F,verbose=T,parallel.sz=11)
es <- es.max
df <- do.call(rbind,lapply(1:ncol(es),function(i){
  dat <- data.frame(
    ro = rownames(es),
    cluster = colnames(es)[i],
    sd.1 = es[,i],
    sd.2 = apply(es[,-i],1,median)
  )
}))
df$fc <- df$sd.1 - df$sd.2
top <- df%>% group_by(cluster) %>% top_n(10,fc)
pdf('./3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaReactome.pdf')
pheatmap(es[unique(top$ro),],show_rownames = T,fontsize = 6,fontsize_col = 8,clustering_method = 'ward.D',angle_col = '45',main = 'GSVA_clusters')
dev.off()
saveRDS(es,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaReactome.rds')
idents <- levels(sub.T)
df$cluster <- factor(df$cluster,levels = levels(sub.T))
res <- split(df,df$cluster)
saveFormat <- lapply(res,function(i){
  i <- i %>% dplyr::filter(fc > 0) %>% dplyr::arrange(desc(fc))
  return(i)
})
write.xlsx(saveFormat,file = './3.1Subcluster_T/Annotate/enrichment/GSVA/gsvaReactome.xlsx',sheetName = idents)
####add module score####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
Idents(sub.T) <- sub.T$seurat_clusters
signature <- list(Naive = c('CCR7','LEF1','TCF7','SELL'),
                  Cytotoxic = c('IFNG','GNLY','GZMA','GZMB','GZMH','NKG7'),
                  Memory = c('IL2','IL7R','FOS','S100A4','ANXA1','CD40LG'),
                  Treg = c('FOXP3','IL2RA','BATF','IKZF2'),
                  Exhausted = c('CTLA4','PDCD1','TIGIT','LAG3','HAVCR2'))
sce <- AddModuleScore(sub.T,features = signature,name = names(signature))
pdf('./3.1Subcluster_T/Annotate/Signature_score.pdf')
VlnPlot(sce,features = c('Naive1','Cytotoxic2','Memory3','Treg4','Exhausted5'),ncol = 2,pt.size = 0,cols = Palettes[['mycols_12']])
FeaturePlot(sce,features = c('Naive1','Cytotoxic2','Memory3','Treg4','Exhausted5'),reduction = 'umap',ncol = 2,order = T,min.cutoff = 0.4,cols = Palettes[['greyMagma']])
dev.off()
sub.T <- sce
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.rds')
####GSEA####
dir.create('./3.1Subcluster_T/Annotate/enrichment/GSEA')
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
gmtfile = './3.1Subcluster_T/Annotate/enrichment/GSEA/c5.go.bp.v2022.1.Hs.symbols.gmt'
geneset <- read.gmt(gmtfile)
deg <- readRDS('./3.1Subcluster_T/Annotate/cluster.all.DE.rds')
deg <- split(deg,deg$cluster)
res <- lapply(deg,function(i){
  genelist <- i$avg_log2FC
  names(genelist) <- toupper(i$gene)
  genelist <- sort(genelist,decreasing = T)
  result <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.9,verbose = F)
  cat(paste0('GSEA of cluster',unique(i$cluster),'_finished','\n'))
  gsea_result_df <- result@result
  return(gsea_result_df)
})
saveRDS(res,file = './3.1Subcluster_T/Annotate/enrichment/GSEA/GSEA_result.rds')
idents <- names(res)
write.xlsx(res,file = './3.1Subcluster_T/Annotate/enrichment/GSEA/GSEA_result.xlsx',sheetName = idents)
res_sig <- lapply(res,function(i){
  res <- i %>% dplyr::filter(pvalue < 0.05&p.adjust < 0.25)
  return(res)
})
write.xlsx(res_sig,file = './3.1Subcluster_T/Annotate/enrichment/GSEA/GSEA_result.sig.xlsx',sheetName = idents)
####annotate celltype####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
pdf('./3.1Subcluster_T/Annotate/cluster0_6.pdf')
plot_density(sub.T,features = c('CD55'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster1.pdf')
plot_density(sub.T,features = c('IL2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster2.pdf')
plot_density(sub.T,features = c('ANXA1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster3.pdf')
plot_density(sub.T,features = c('CCL4'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster5.pdf')
plot_density(sub.T,features = c('RTKN2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster6.pdf')
plot_density(sub.T,features = c('NELL2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster7.pdf')
plot_density(sub.T,features = c('HSPA1A'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster8.pdf')
plot_density(sub.T,features = c('IL2RA'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster9.pdf')
plot_density(sub.T,features = c('PDCD1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster10.pdf')
plot_density(sub.T,features = c('IFI44L'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/cluster11.pdf')
plot_density(sub.T,features = c('HAVCR2'),reduction = 'umap')
dev.off()
####annotate
minor_celltype <- sub.T@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "CD4Tn_CD55", minor_celltype)
minor_celltype <- gsub("^1$", "CD4Tcm_IL2", minor_celltype)
minor_celltype <- gsub("^2$", "CD4Tem_ANXA1", minor_celltype)
minor_celltype <- gsub("^3$", "CD8Tem_CCL4", minor_celltype)
minor_celltype <- gsub("^4$", "CD4Tun", minor_celltype)
minor_celltype <- gsub("^5$", "CD4Treg_RTKN2", minor_celltype)
minor_celltype <- gsub("^6$", "CD8Tn_CD55", minor_celltype)
minor_celltype <- gsub("^7$", "CD4Tstd_HSPA1A", minor_celltype)
minor_celltype <- gsub("^8$", "CD4Treg_IL2RA", minor_celltype)
minor_celltype <- gsub("^9$", "CD4Tex_PDCD1", minor_celltype)
minor_celltype <- gsub("^10$", "CD4Tif1_IFI44L", minor_celltype)
minor_celltype <- gsub("^11$", "CD8Tex_HAVCR2", minor_celltype)

table(minor_celltype)
sub.T <- AddMetaData(sub.T,minor_celltype,col.name = 'minor_celltype')
table(sub.T$minor_celltype)
sub.T$minor_celltype <- factor(sub.T$minor_celltype,levels = names(table(sub.T$minor_celltype)))
Idents(sub.T) <- sub.T$minor_celltype
DefaultAssay(sub.T) <- "RNA"

source('./function/do_dimplot.R')
pdf("3.1Subcluster_T/Annotate/cellType.pdf",height = 10,width = 12)
DoDimplot(sub.T,groupby = 'minor_celltype',colors = Palettes[['mycols_12_2']])
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.T, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = sub.T, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_12_2']])
ratio.plot(seurat.object = sub.T, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.T,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_12_2']])
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.T,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_12_2']],ncol = 2)
dev.off()

saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.pro.rds')

####expression of celltype markers####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.pro.rds')
genelist <- c('CD4','CD8A','CD8B','CCR7','LEF1','TCF7','CD55','NELL2','IL2','CD40LG','IL7R','ANXA1','S100A4','FOS','JUN',
              'HLA-DRB1','HLA-DPB1','CCL4','CCL5','GZMA','GNLY','GZMK','GZMB','GZMH','NKG7','RTKN2','FOXP3','IKZF2','MAF','IL2RA','BATF',
              'ICA1','TOX2','PDCD1','CTLA4','TIGIT','LYST','LGALS3','HAVCR2','HSPA1A','HSPA1B','HSPD1','HSPE1','DNAJB1','DNAJB4',
              'IFI44L','IFI6','IFI44','IFIT1','IFIT3')
sub.T <- ScaleData(sub.T,features = c(VariableFeatures(sub.T),genelist))
sub.T@meta.data$celltype <- sub.T@meta.data$minor_celltype
sub.T@meta.data$celltype <- factor(sub.T@meta.data$celltype,levels = c('CD4Tn_CD55','CD8Tn_CD55','CD4Tcm_IL2','CD4Tem_ANXA1',
                                                               'CD8Tem_CCL4','CD4Treg_RTKN2',
                                                               'CD4Treg_IL2RA','CD4Tex_PDCD1','CD8Tex_HAVCR2',
                                                               'CD4Tstd_HSPA1A','CD4Tif1_IFI44L'),ordered = T)
Idents(sub.T) <- sub.T$celltype
exp.matrix <- GetAssayData(sub.T, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.T$celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('./3.1Subcluster_T/Annotate/celltype_expression.pdf',width = 10,height = 18)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(30, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

# sub.T@meta.data$minor_celltype <- factor(sub.T@meta.data$minor_celltype,levels = c('CD4Tn_CCR7','CD8Tn_CCR7','CD4Tcm_IL2','CD4Tem_ANXA1',
#                                                                        'CD8Tem_GZMK','CD4Treg_RTKN2',
#                                                                        'CD4Treg_IL2RA','CD4Tex_PDCD1','CD8Tex_HAVCR2',
#                                                                        'CD4Tstd_HSPA1A','CD4Tif1_IFI44L','CD4Tun'),ordered = T)
# Idents(sub.T) <- sub.T$minor_celltype
# cluster.exp <- apply(gene.matrix, 1, function(x){
#   a <- tapply(x, sub.T$minor_celltype, mean)
#   return(a)
# })
# cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
# 
# pdf('./3.1Subcluster_T/Annotate/celltype_expression_2.pdf',width = 10,height = 18)
# Heatmap(t(cluster.score.normailzed), 
#         width = unit(15, "cm"), height = unit(20, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
#         show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6),rect_gp = gpar(col = 'white'),
#         column_names_gp = gpar(fontsize = 8),
#         heatmap_legend_param = list(
#           title = "Expression", at = c(0, 1), 
#           labels = c("min", "max")),col = Palettes[['blueRed']])
# dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_expression.dotplot.pdf',width = 15,height = 8)
Idents(sub.T) <- sub.T$celltype
DotPlot(subset(sub.T,subset = minor_celltype!='CD4Tun'),features = genelist,cols = c('grey','red'),scale.by = 'size',col.min = 0,scale.min = -5) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(sub.T) <- sub.T$minor_celltype
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.pro.rds')
