#' @description: annotate cycling subclusters

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
sub.scRNA <- subset(data.merge, subset = major_celltype == 'Cycling')
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
dir.create('./3.3Subcluster_cycling')
sub.list.Standard <- variableFeatureSelection(seurat.lists = sub.list, method = "Standard", nfeatures = 3000)
saveRDS(sub.list.Standard, file = "3.3Subcluster_cycling/sub.list.Standard.3000.rds")

sub.list.SCT <- variableFeatureSelection(seurat.lists = sub.list, method = "SCT", nfeatures = 3000,return.only.var.genes = T,vars.to.regress = c("nCount_RNA", "percent.mt"))
saveRDS(sub.list.SCT, file = "3.3Subcluster_cycling/sub.list.SCT.3000.rds")

## Remove previous clustering results
sub.scRNA <- merge(sub.list.SCT[[1]], y = sub.list.SCT[2:length(sub.list.SCT)], project = "Cycling cells")
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
####RNA.Harmony.PC25####
pdf("3.3Subcluster_cycling/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()

####RNA.Harmony.PC20####
pdf("3.3Subcluster_cycling/RNA.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()

####select: RNA.Harmony.PC20####
pdf("3.3Subcluster_cycling/RNA.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.3Subcluster_cycling/sub.scRNA.harmony.rds')

sub.cycling = readRDS('./3.3Subcluster_cycling/sub.scRNA.harmony.rds')
####select res:0.3####
sub.cycling$seurat_clusters = sub.cycling$RNA_snn_res.0.3
Idents(sub.cycling)=sub.cycling$seurat_clusters
dir.create('3.3Subcluster_cycling/Annotate')
pdf('3.3Subcluster_cycling/Annotate/cluster.pdf')
DimPlot(sub.cycling,group.by = 'seurat_clusters',cols = Palettes[['mycols_4']],reduction = 'umap',label = T)
DimPlot(sub.cycling,group.by = 'seurat_clusters',cols = Palettes[['mycols_4']],reduction = 'tsne',label = T)
dev.off()

####observe T&B cell markers expression####
pdf('./3.3Subcluster_cycling/Annotate/observe_expression.pdf')
FeaturePlot(sub.cycling,features = c('CD3D','CD3E','MS4A1','CD79A','CD79B','CD19'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']],min.cutoff = 1)
VlnPlot(sub.cycling,features = c('CD3D','CD3E','CD3G','MS4A1','CD79A','CD79B'),group.by = 'seurat_clusters')
FeaturePlot(sub.cycling,features = c('LRMP','SUGCT','AICDA','RGS13'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']],min.cutoff = 1)
FeaturePlot(sub.cycling,features = c('IL7R','CD4','CD8A','CD8B'),reduction = 'umap',order = T,cols = Palettes[['greyMagma']],min.cutoff = 1)
dev.off()

saveRDS(sub.cycling,file = './3.3Subcluster_cycling/sub.cycling.rds')

####differential expression gene####
sub.cycling <- readRDS('./3.3Subcluster_cycling/sub.cycling.rds')
Idents(sub.cycling) <- sub.cycling$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.cycling, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.25,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.3Subcluster_cycling/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.3Subcluster_cycling/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.3Subcluster_cycling/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sub.cycling <- ScaleData(sub.cycling,features = c(VariableFeatures(sub.cycling),unique(top5.genes$gene)))
pdf("3.3Subcluster_cycling/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(sub.cycling, features = unique(top5.genes$gene), size = 2,group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

top8.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
sub.cycling <- ScaleData(sub.cycling,features = c(VariableFeatures(sub.cycling),unique(top8.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('3.3Subcluster_cycling/Annotate/cluster.top8genes.mean.pdf',height = 12)
DoHeatmap_my(sub.cycling,features = unique(top8.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(sub.cycling,file = './3.2Subcluster_B/sub.cycling.rds')

####celltypist predict####
sub.cycling <- readRDS('./3.3Subcluster_cycling/sub.cycling.rds')
library(reticulate)
py_config()
scanpy <- import('scanpy')
pandas <- import('pandas')
numpy <- import('numpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(sub.cycling[['RNA']]@counts)))),
                       obs = pandas$DataFrame(sub.cycling@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(sub.cycling[['RNA']]@counts),
                                                         row.names = rownames(sub.cycling[['RNA']]@counts)))
)
adata
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)

predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- sub.cycling@meta.data$seurat_clusters
sub.cycling <- AddMetaData(sub.cycling,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(sub.cycling$celltypist_minor,sub.cycling$seurat_clusters)
write.xlsx(predicted_minor,file = '3.3Subcluster_cycling/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(sub.cycling,file = './3.2Subcluster_B/sub.cycling.rds')
####annotate cell type####
sub.cycling <- readRDS('./3.3Subcluster_cycling/sub.cycling.rds')
minor_celltype <- sub.cycling@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "Cycling Tm_IL7R", minor_celltype)
minor_celltype <- gsub("^1$", "Cycling Bgc_LRMP", minor_celltype)
minor_celltype <- gsub("^2$", "Undefined", minor_celltype)
minor_celltype <- gsub("^3$", "Cycling B_BANK1", minor_celltype)
table(minor_celltype)
sub.cycling <- AddMetaData(sub.cycling,minor_celltype,col.name = 'minor_celltype')
table(sub.cycling$minor_celltype)
sub.cycling$minor_celltype <- factor(sub.cycling$minor_celltype,levels = names(table(sub.cycling$minor_celltype)))
Idents(sub.cycling) <- sub.cycling$minor_celltype
DefaultAssay(sub.cycling) <- "RNA"

source('./function/do_dimplot.R')
pdf("3.3Subcluster_cycling/Annotate/cellType.pdf",width = 8,height = 6)
DoDimplot(sub.cycling,groupby = 'minor_celltype',colors = Palettes[['mycols_4']])
dev.off()

pdf('3.3Subcluster_cycling/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.cycling, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = sub.cycling, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_4']])
ratio.plot(seurat.object = sub.cycling, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('3.3Subcluster_cycling/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.cycling,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_4']])
dev.off()

pdf('3.3Subcluster_cycling/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.cycling,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_4']],ncol = 2)
dev.off()

saveRDS(sub.cycling,file = '3.3Subcluster_cycling/sub.cycling.pro.rds')
####expression of celltype markers####
sub.cycling <- readRDS('./3.3Subcluster_cycling/sub.cycling.pro.rds')
genelist <- c('CD3D','CD3E','CD3G','TRBC1','TRBC2','IL7R','S100A4','ANXA1','BANK1','PLAC8','MS4A1','CD79A','CD79B',
              'LRMP','RGS13','AICDA','SUGCT','FABP5','HSPE1','NME1','HSP90AB1')
sub.cycling <- ScaleData(sub.cycling,features = c(VariableFeatures(sub.cycling),genelist))
sub.cycling@meta.data$celltype <- sub.cycling@meta.data$minor_celltype
sub.cycling@meta.data$celltype <- factor(sub.cycling@meta.data$celltype,levels = c('Cycling Tm_IL7R','Cycling B_BANK1','Cycling Bgc_LRMP',
                                                                       'Undefined'),ordered = T)
Idents(sub.cycling) <- sub.cycling$celltype
exp.matrix <- GetAssayData(sub.cycling, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.cycling$celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform standardization for all clusters per gene
pdf('3.3Subcluster_cycling/Annotate/celltype_expression.pdf',width = 12,height = 22)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(24, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

Idents(sub.cycling) <- sub.cycling$celltype
pdf('3.3Subcluster_cycling/Annotate/celltype_expression.dotplot.pdf',width = 10)
DotPlot(sub.cycling,features = genelist,cols = c('grey','red'),col.min = 0,scale.by = 'size') + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(sub.cycling) <- sub.cycling$minor_celltype
saveRDS(sub.cycling,file = '3.3Subcluster_cycling/sub.cycling.pro.rds')
