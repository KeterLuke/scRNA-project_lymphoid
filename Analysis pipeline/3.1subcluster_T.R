

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
data.merge <- readRDS("./2.Cluster/data.merge.pro.rds")
table(Idents(data.merge))
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
rm(sub.list,sub.list.SCT,sub.list.Standard)
####correct batch effect####
# ####SCT.Harmony.PC20####
# pdf("3.1Subcluster_T/SCT.Harmony.PC20.Integration.pdf")
# sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
# dev.off()
# 
# ####SCT.Harmony.PC30####
# pdf("3.1Subcluster_T/SCT.Harmony.PC30.Integration.pdf")
# sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
# dev.off()
# 
# ####SCT.Harmony.PC40####
# pdf("3.1Subcluster_T/SCT.Harmony.PC40.Integration.pdf")
# sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 40, npcs = 100)
# dev.off()
# 
# ####SCT.Harmony.PC45####
# pdf("3.1Subcluster_T/SCT.Harmony.PC45.Integration.pdf")
# sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 45, npcs = 100)
# dev.off()
# 
# ####SCT.Harmony.PC50####
# pdf("3.1Subcluster_T/SCT.Harmony.PC50.Integration.pdf")
# sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "SCT", set.resolutions = seq(0.1, 1, by = 0.1), PC = 50, npcs = 100)
# dev.off()

####RNA.Harmony.PC20####
pdf("3.1Subcluster_T/RNA.Harmony.PC20.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC20.rds')
####RNA.Harmony.PC25####
pdf("3.1Subcluster_T/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC25.rds')
####RNA.Harmony.PC30####
pdf("3.1Subcluster_T/RNA.Harmony.PC30.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC30.rds')


##select RNA.Harmony.PC20 res:0.6
sub.T <- readRDS('./3.1Subcluster_T/RNA.Harmony.Integration.PC20.rds')
sub.T$seurat_clusters <- sub.T$RNA_snn_res.0.6
length(unique(sub.T$seurat_clusters))
table(sub.T$seurat_clusters)
dir.create('./3.1Subcluster_T/Annotate')
pdf('./3.1Subcluster_T/Annotate/cluster.pdf')
DimPlot(sub.T,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_15']])
DimPlot(sub.T,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_15']])
dev.off()
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.rds')
#### classical marker#### 
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
pdf('./3.1Subcluster_T/Annotate/cd4.cd8.observe.pdf')
FeaturePlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0,ncol = 2)
VlnPlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),group.by = 'seurat_clusters',ncol = 2)
dev.off()

pdf('./3.1Subcluster_T/Annotate/B.Observe.pdf')
FeaturePlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0,ncol = 2)
VlnPlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),group.by = 'seurat_clusters',ncol = 2)
dev.off()

genelist <- c('CD4','CD8A','CD8B','IFNG','GZMK')
Idents(sub.T) <- sub.T$seurat_clusters
exp.matrix <- GetAssayData(sub.T, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.T$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform standardization for all clusters per gene
pdf('./3.1Subcluster_T/Annotate/marker_expression.pdf',width = 10,height = 5)
Heatmap(t(cluster.score.normailzed), 
       column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Expression", at = c(0,1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()
##cluster similarity#####
expMatrix <- GetAssayData(sub.T, slot = "scale.data")
highVariableGenes <- VariableFeatures(sub.T)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, sub.T$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.1Subcluster_T/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.1Subcluster_T/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
#### Differential expression####
Idents(sub.T) <- sub.T$seurat_clusters
table(Idents(sub.T))
cluster.pos.markers <- FindAllMarkers(sub.T, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)
####remove cluster11(express B markers)
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
sub.T <- subset(sub.T,seurat_clusters != '11')
table(sub.T$seurat_clusters)
####rerun integration####
sub.T <- NormalizeData(sub.T)
sub.T <- FindVariableFeatures(sub.T,nfeatures = 3000)
sub.T <- ScaleData(sub.T,features = VariableFeatures(sub.T),vars.to.regress = c("nCount_RNA", "percent.mt"))
pdf("3.1Subcluster_T/RNA.Harmony.PC20.pro.PC20.pdf")
sub.T.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.T, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(sub.T.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC20.pro.PC20.rds')

pdf("3.1Subcluster_T/RNA.Harmony.PC20.pro.PC25.pdf")
sub.T.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.T, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()
saveRDS(sub.T.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC20.pro.PC25.rds')

pdf("3.1Subcluster_T/RNA.Harmony.PC20.pro.PC30.pdf")
sub.T.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.T, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()
saveRDS(sub.T.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC20.pro.PC30.rds')
####select:RNA harmony pc30,res:0.6####
sub.T <- readRDS('./3.1Subcluster_T/RNA.Harmony.Integration.PC20.pro.PC30.rds')
sub.T$seurat_clusters <- sub.T$RNA_snn_res.0.6
table(sub.T$seurat_clusters)
pdf('./3.1Subcluster_T/Annotate/cluster.pdf')
DimPlot(sub.T,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
DimPlot(sub.T,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
dev.off()
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.reclustered.rds')
#### classical marker#### 
sub.T <- readRDS('./3.1Subcluster_T/sub.T.reclustered.rds')
pdf('./3.1Subcluster_T/Annotate/cd4.cd8.observe.pdf')
FeaturePlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0,ncol = 2)
VlnPlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),group.by = 'seurat_clusters',ncol = 2)
dev.off()

pdf('./3.1Subcluster_T/Annotate/B.Observe.pdf')
FeaturePlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0,ncol = 2)
VlnPlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),group.by = 'seurat_clusters',ncol = 2)
dev.off()

genelist <- c('CD4','CD8A','CD8B','IFNG','GZMK')
Idents(sub.T) <- sub.T$seurat_clusters
exp.matrix <- GetAssayData(sub.T, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.T$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform standardization for all clusters per gene
pdf('./3.1Subcluster_T/Annotate/marker_expression.pdf',width = 10,height = 5)
Heatmap(t(cluster.score.normailzed), 
        column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Expression", at = c(0,1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()
##cluster similarity#####
expMatrix <- GetAssayData(sub.T, slot = "scale.data")
highVariableGenes <- VariableFeatures(sub.T)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, sub.T$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.1Subcluster_T/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.1Subcluster_T/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
#### Differential expression####
Idents(sub.T) <- sub.T$seurat_clusters
table(Idents(sub.T))
cluster.pos.markers <- FindAllMarkers(sub.T, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)
####annotate cd4 and cd8 clusters####
T_major <-  sub.T@meta.data$seurat_clusters
T_major <- gsub("^3$", "CD8", T_major)
T_major <- gsub("^6$", "CD8", T_major)
T_major <- gsub("^11$", "CD8", T_major)
T_major[T_major!='CD8'] <- 'CD4'
table(T_major)
sub.T <- AddMetaData(sub.T,T_major,col.name = 'T_major')
table(sub.T$T_major)
sub.T$T_major <- factor(sub.T$T_major,levels = names(table(sub.T$T_major)))
Idents(sub.T) <- sub.T$T_major
DefaultAssay(sub.T) <- "RNA"

pdf("3.1Subcluster_T/Annotate/cellType.pdf",height = 10,width = 12)
DimPlot(sub.T,group.by = 'T_major',reduction = 'umap',label = T)
DimPlot(sub.T,group.by = 'T_major',reduction = 'tsne',label = T)
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.T, id.vars1 = "sample", id.vars2 = "T_major", angle = 60)
ratio.plot(seurat.object = sub.T, id.vars1 = "T_major", id.vars2 = "sample", angle = 60)
ratio.plot(seurat.object = sub.T, id.vars1 = "group", id.vars2 = "T_major", angle = 60)
ratio.plot(seurat.object = sub.T, id.vars1 = "T_major", id.vars2 = "group", angle = 60)
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.T,reduction = 'umap',group.by = 'T_major',split.by = 'group')
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.T,reduction = 'umap',group.by = 'T_major',split.by = 'sample',ncol = 2)
dev.off()

saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.reclustered.pro.rds')
####test of celltype percentage between groups####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.reclustered.pro.rds')
df <- as.data.frame(table(sub.T$T_major,sub.T$sample))
df$percent <- round(df$Freq/ncol(sub.T),4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('./3.1Subcluster_T/Annotate/celltype.percentage.group.pdf')
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    stat_compare_means(comparisons = comparison,method = 'wilcox.test',paired = T) +
    theme_classic() +
    labs(title = i,x = 'Group',y = '% of all T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype.percentage.group.unpaired.pdf')
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggplot(tmp,aes(x = group,y = percent,fill = group)) + geom_boxplot(width = .4) + geom_point() + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of all T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('./3.1Subcluster_T/Annotate/celltype.percentage.group.barplot.pdf')
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggbarplot(tmp,x = 'group',y = 'percent',fill = 'group',add = 'mean_se') + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of all T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.reclustered.pro.rds')

####split T cells into CD4/CD8 subclusters####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.reclustered.pro.rds')
table(Idents(sub.T))
CD4 <- subset(sub.T,T_major == 'CD4')
CD8 <- subset(sub.T,T_major == 'CD8')
CD4 <- DietSeurat(CD4)
CD8 <- DietSeurat(CD8)
dir.create('./3.1Subcluster_T/Annotate/CD4')
dir.create('./3.1Subcluster_T/Annotate/CD8')
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.rds')
