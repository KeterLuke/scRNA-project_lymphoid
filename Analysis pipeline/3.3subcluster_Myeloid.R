#' @description: annotate myeloid subclusters

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
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
source('./function/variableFeatureSelection.R')
source('./function/Integrate_multisamples.R')
####primarily annotated results####
data.merge <- readRDS("./2.Cluster/data.merge.pro.rds")
sub.scRNA <- subset(data.merge, subset = major_celltype == 'Myeloid')
rm(data.merge)
DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- DietSeurat(sub.scRNA, assays = "RNA")####remove SCT assay, scale.data and dimension reduction 
## Remove previous clustering results
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(sub.scRNA@meta.data))
sub.scRNA@meta.data <- sub.scRNA@meta.data[,-index]
View(sub.scRNA@meta.data)

DefaultAssay(sub.scRNA) <- "RNA"
sub.scRNA <- NormalizeData(sub.scRNA, verbose = FALSE)
sub.scRNA <- FindVariableFeatures(sub.scRNA,nfeatures = 3000,verbose = F)
sub.scRNA <- ScaleData(sub.scRNA, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(sub.scRNA))
####correct batch effect####
dir.create('./3.3Subcluster_myeloid')
sub.scRNA <- RunPCA(sub.scRNA,features = VariableFeatures(sub.scRNA),npcs = 100,verbose = F)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(sub.scRNA,harmony = F)
####RNA.Harmony.PC15####
pdf("3.3Subcluster_myeloid/RNA.Harmony.PC15.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 15, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.3Subcluster_myeloid/RNA.Harmony.Integration.PC15.rds')
####RNA.Harmony.PC10####
pdf("3.3Subcluster_myeloid/RNA.Harmony.PC10.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 10, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.3Subcluster_myeloid/RNA.Harmony.Integration.PC10.rds')
####RNA.Harmony.PC5####
pdf("3.3Subcluster_myeloid/RNA.Harmony.PC5.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 5, npcs = 100)
dev.off()
saveRDS(sub.scRNA.harmony,file = './3.3Subcluster_myeloid/RNA.Harmony.Integration.PC5.rds')
####select RNA Harmony pc5 resolution:0.3####
sub.m <- readRDS('./3.3Subcluster_myeloid/RNA.Harmony.Integration.PC5.rds')
sub.m$seurat_clusters <- sub.m$RNA_snn_res.0.3
table(sub.m$seurat_clusters)
saveRDS(sub.m,file = './3.3Subcluster_myeloid/sub.m.rds')

Idents(sub.m) <- sub.m$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.m, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
table(sub.m$seurat_clusters)
length(unique(sub.m$seurat_clusters))
dir.create('./3.3Subcluster_myeloid/Annotate')
pdf('./3.3Subcluster_myeloid/Annotate/cluster.pdf')
DimPlot(sub.m,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']],shuffle = T)
DimPlot(sub.m,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']],shuffle = T)
dev.off()
saveRDS(sub.m,file = './3.3Subcluster_myeloid/sub.m.rds')
#### Differential expression####
sub.m <- readRDS('./3.3Subcluster_myeloid/sub.m.rds')
Idents(sub.m) <- sub.m$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.m, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.25,logfc.threshold = 0.25)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.3Subcluster_myeloid/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.3Subcluster_myeloid/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.3Subcluster_myeloid/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sub.m <- ScaleData(sub.m,features = c(VariableFeatures(sub.m),unique(top5.genes$gene)))
pdf("3.3Subcluster_myeloid/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(sub.m, features = unique(top5.genes$gene), size = 2,group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['blueRed']])
dev.off()

top8.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
sub.m <- ScaleData(sub.m,features = c(VariableFeatures(sub.m),unique(top8.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('./3.3Subcluster_myeloid/Annotate/cluster.top8genes.mean.pdf',height = 12)
DoHeatmap_my(sub.m,features = unique(top8.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(sub.m,file = './3.3Subcluster_myeloid/sub.m.rds')
####annotate celltype####
sub.m <- readRDS('./3.3Subcluster_myeloid/sub.m.rds')
####annotate
minor_celltype <- sub.m@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "cDC2_CD1C", minor_celltype)
minor_celltype <- gsub("^1$", "cDC1_CLEC9A", minor_celltype)
minor_celltype <- gsub("^2$", "Doublet", minor_celltype)
minor_celltype <- gsub("^3$", "Macro_C1QA", minor_celltype)



table(minor_celltype)
sub.m <- AddMetaData(sub.m,minor_celltype,col.name = 'minor_celltype')
table(sub.m$minor_celltype)
sub.m$minor_celltype <- factor(sub.m$minor_celltype,levels = names(table(sub.m$minor_celltype)))
Idents(sub.m) <- sub.m$minor_celltype
DefaultAssay(sub.m) <- "RNA"

source('./function/do_dimplot.R')
pdf("3.3Subcluster_myeloid/Annotate/cellType.pdf")
DoDimplot(sub.m,groupby = 'minor_celltype',colors = Palettes[['summerNight']])
dev.off()

pdf('./3.3Subcluster_myeloid/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.m, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = sub.m, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['summerNight']])
ratio.plot(seurat.object = sub.m, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('./3.3Subcluster_myeloid/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.m,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['summerNight']])
dev.off()

pdf('./3.3Subcluster_myeloid/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.m,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['summerNight']],ncol = 2)
dev.off()

saveRDS(sub.m,file = './3.3Subcluster_myeloid/sub.m.pro.rds')
####expression of celltype markers####
sub.m <- readRDS('./3.3Subcluster_myeloid/sub.m.pro.rds')
genelist <- c('IDO1','CLEC9A','XCR1','IRF8','CD1C','CD1E','CLEC10A','FCER1A','CD3D','CD3G','MS4A1','CD79A','CD14','CD68','C1QA','APOE'
)

exp.matrix <- GetAssayData(sub.m, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.m$minor_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('./3.3Subcluster_myeloid/Annotate/celltype_expression.pdf',width = 12,height = 22)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(24, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

Idents(sub.m) <- sub.m$minor_celltype
pdf('./3.3Subcluster_myeloid/Annotate/celltype_expression.dotplot.pdf',width = 10)
DotPlot(sub.m,features = genelist,scale.by = 'size') + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(sub.m) <- sub.m$minor_celltype
saveRDS(sub.m,file = './3.3Subcluster_myeloid/sub.m.pro.rds')
####test of celltype percentage between groups####
sub.m <- readRDS('./3.3Subcluster_myeloid/sub.m.pro.rds')
df <- as.data.frame(table(sub.m$minor_celltype,sub.m$sample))
df$total <- apply(df,1,function(x){sum(df$Freq[df$Var2==x[2]])})
df$percent <- round(df$Freq/df$total,4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('./3.3Subcluster_myeloid/Annotate/celltype.percentage.group.pdf',height = 12,width = 15)
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
