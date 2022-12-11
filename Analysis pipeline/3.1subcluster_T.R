

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

## Remove previous clustering results
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(sub.scRNA@meta.data))
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-index]
View(sub.scRNA@meta.data)

DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- NormalizeData(sub.scRNA, verbose = FALSE)
sub.scRNA <- FindVariableFeatures(sub.scRNA,nfeatures = 3000,verbose = F)
sub.scRNA <- ScaleData(sub.scRNA, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(sub.scRNA))

####correct batch effect####
sub.scRNA <- RunPCA(sub.scRNA,features = VariableFeatures(sub.scRNA),npcs = 100,verbose = F)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(sub.scRNA,harmony = F)
####RNA.Harmony.PC15####
pdf("3.1Subcluster_T/RNA.Harmony.PC15.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 15, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC15.rds')
####RNA.Harmony.PC10####
pdf("3.1Subcluster_T/RNA.Harmony.PC10.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 10, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.1Subcluster_T/RNA.Harmony.Integration.PC10.rds')
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

##select RNA.Harmony.PC15 res:0.6
sub.T <- readRDS('./3.1Subcluster_T/RNA.Harmony.Integration.PC15.rds')
sub.T$seurat_clusters <- sub.T$RNA_snn_res.0.6
length(unique(sub.T$seurat_clusters))
table(sub.T$seurat_clusters)
dir.create('./3.1Subcluster_T/Annotate')
pdf('./3.1Subcluster_T/Annotate/cluster.pdf')
DimPlot(sub.T,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
DimPlot(sub.T,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']])
dev.off()
saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.rds')
#### classical marker#### 
sub.T <- readRDS('./3.1Subcluster_T/sub.T.rds')
pdf('./3.1Subcluster_T/Annotate/cd4.cd8.observe.pdf')
FeaturePlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0.5,ncol = 2)
VlnPlot(sub.T,features = c('CD4','CD8A','CD8B','TRDC'),group.by = 'seurat_clusters',ncol = 2,cols = Palettes[['mycols_12']])
dev.off()

pdf('./3.1Subcluster_T/Annotate/B.Observe.pdf')
FeaturePlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),reduction = 'umap',cols = Palettes[['greyMagma']],
            order = T,min.cutoff = 0,ncol = 2)
VlnPlot(sub.T,features = c('MS4A1','CD19','CD79A','CD79B'),group.by = 'seurat_clusters',ncol = 2,Palettes[['mycols_12']])
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

#### Differential expression####
Idents(sub.T) <- sub.T$seurat_clusters
table(Idents(sub.T))
cluster.pos.markers <- FindAllMarkers(sub.T, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

####annotate cd4 and cd8 clusters####
T_major <-  sub.T@meta.data$seurat_clusters
T_major <- gsub("^3$", "CD8", T_major)
T_major <- gsub("^8$", "CD8", T_major)
T_major[T_major!='CD8'] <- 'CD4'
table(T_major)
sub.T <- AddMetaData(sub.T,T_major,col.name = 'T_major')
table(sub.T$T_major)
sub.T$T_major <- factor(sub.T$T_major,levels = names(table(sub.T$T_major)))
Idents(sub.T) <- sub.T$T_major
DefaultAssay(sub.T) <- "RNA"

pdf("3.1Subcluster_T/Annotate/cellType.pdf",height = 10,width = 12)
DimPlot(sub.T,group.by = 'T_major',reduction = 'umap',label = T,cols = Palettes[['mycols_4']],label.size = 8)
DimPlot(sub.T,group.by = 'T_major',reduction = 'tsne',label = T,cols = Palettes[['mycols_4']],label.size = 8)
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

saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.pro.rds')
####cell marker expression####
Idents(sub.T) <- sub.T$T_major
celltypemarkers <- FindAllMarkers(sub.T,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
celltypemarkers.sig <- celltypemarkers[which(celltypemarkers$p_val_adj<0.05),]
genelist <- c('CD4','CD8A','CD8B','GZMK')
pdf('./3.1Subcluster_T/Annotate/celltype.marker.expression.violin.pdf')
VlnPlot(sub.T,features = genelist,group.by = 'T_major',stack = T,pt.size = 0,fill.by = 'ident',flip = T) + 
  ggsci::scale_fill_lancet()
dev.off()
####test of celltype percentage between groups####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.pro.rds')
df <- as.data.frame(table(sub.T$T_major,sub.T$sample))
df$total <- apply(df,1,function(x){sum(df$Freq[df$Var2==x[2]])})
df$percent <- round(df$Freq/df$total,4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('./3.1Subcluster_T/Annotate/celltype.percentage.group.pdf')
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    theme_classic() + scale_y_continuous(breaks = scales::pretty_breaks(n=10)) + 
    labs(title = i,x = 'Group',y = '% of all T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()


saveRDS(sub.T,file = './3.1Subcluster_T/sub.T.pro.rds')

####split T cells into CD4/CD8 subclusters####
sub.T <- readRDS('./3.1Subcluster_T/sub.T.reclustered.pro.rds')
table(Idents(sub.T))
CD4 <- subset(sub.T,T_major == 'CD4')
CD8 <- subset(sub.T,T_major == 'CD8')
CD4 <- DietSeurat(CD4,assays = 'RNA')
CD8 <- DietSeurat(CD8,assays = 'RNA')
dir.create('./3.1Subcluster_T/Annotate/CD4')
dir.create('./3.1Subcluster_T/Annotate/CD8')
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
saveRDS(CD8,file = './3.1Subcluster_T/Annotate/CD8/CD8.rds')
