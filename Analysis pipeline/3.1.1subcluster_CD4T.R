#' @description: annotate CD4T subclusters

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
library(speckle)
library(cowplot)
set.seed(101)
library(future)
plan("multisession", workers = 1) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
source('./function/variableFeatureSelection.R')
source('./function/Integrate_multisamples.R')
####primarily annotated results####
CD4 <- readRDS("./3.1Subcluster_T/Annotate/CD4/CD4.rds")
table(Idents(CD4))
## Remove previous clustering results
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(CD4@meta.data))
CD4@meta.data <- CD4@meta.data[,-index]
View(CD4@meta.data)
# assay=RNA
DefaultAssay(CD4) <- "RNA"
CD4 <- NormalizeData(CD4, verbose = FALSE)
CD4 <- FindVariableFeatures(CD4,nfeatures = 3000,verbose = F)
CD4 <- ScaleData(CD4, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(CD4))
####correct batch effect####
####RNA.Harmony.PC15####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC15.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 15, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC15.rds')
####RNA.Harmony.PC20####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC20.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 20, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC20.rds')
####RNA.Harmony.PC25####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC25.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC25.rds')
####RNA.Harmony.PC12####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC12.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 12, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC12.rds')
####RNA.Harmony.PC30####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC30.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC30.rds')
####RNA.Harmony.PC35####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC35.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 35, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC35.rds')
##select RNA.Harmony.PC15 res:0.7####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC15.rds')
CD4$seurat_clusters <- CD4$RNA_snn_res.0.7
length(unique(CD4$seurat_clusters))
table(CD4$seurat_clusters)
dir.create('./3.1Subcluster_T/Annotate/CD4/Annotate')
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster.pdf')
DimPlot(CD4,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']],shuffle = T)
DimPlot(CD4,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_12']],shuffle = T)
dev.off()
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
#### Differential expression####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
Idents(CD4) <- CD4$seurat_clusters
table(Idents(CD4))
cluster.pos.markers <- FindAllMarkers(CD4, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1)
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
CD4 <- ScaleData(CD4,features = c(VariableFeatures(CD4),unique(top5.genes$gene)))
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(subset(CD4,downsample = 50), features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['mycols_12']],group.bar = T) + NoLegend()
dev.off()

top10.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
CD4 <- ScaleData(CD4,features = c(VariableFeatures(CD4),unique(top10.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/cluster.top10genes.mean.pdf',height = 12)
DoHeatmap_my(CD4,features = unique(top10.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
##cluster similarity#####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
expMatrix <- GetAssayData(CD4, slot = "scale.data")
highVariableGenes <- VariableFeatures(CD4)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, CD4$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
##expression of functional markers####
cell.type.markers <- read.table(file = "3.1Subcluster_T/Annotate/CD4/Annotate/CD4_markers.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(CD4, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, CD4$seurat_clusters, mean)
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
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/cluster.functional.signature.pdf",height = 10)
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
####add module score####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
Idents(CD4) <- CD4$seurat_clusters
signature <- list(Naive = c('CCR7','CD55','LEF1','TCF7'),
                  Memory = c('IL7R','CD40LG','IL2','FOS','ANXA1','S100A4','JUN'),
                  Treg = c('FOXP3','IL2RA','BATF','IKZF2'),
                  Exhaustion = c('PDCD1','TIGIT','LAG3','HAVCR2','TOX','TOX2','ITGAE','LYST','CTLA4'))
sce <- AddModuleScore(CD4,features = signature,name = names(signature))
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/Signature_score.pdf',width = 8,height = 12)
VlnPlot(sce,features = c('Naive1','Memory2','Treg3','Exhaustion4'),ncol = 2,pt.size = 0,cols = Palettes[['mycols_12']])
FeaturePlot(sce,features = c('Naive1','Memory2','Treg3','Exhaustion4'),reduction = 'umap',ncol = 2,order = T,min.cutoff = 0,cols = Palettes[['greyMagma']])
dev.off()
CD4 <- sce
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
####celltypist predict####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
library(reticulate)
py_config()
numpy <- import('numpy')
pandas <- import('pandas')
scanpy <- import('scanpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(CD4[['RNA']]@counts)))),
                       obs = pandas$DataFrame(CD4@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(CD4[['RNA']]@counts),
                                                         row.names = rownames(CD4[['RNA']]@counts)))
)
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- CD4@meta.data$seurat_clusters
CD4 <- AddMetaData(CD4,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(CD4$celltypist_minor,CD4$seurat_clusters)

write.xlsx(predicted_minor,file = '3.1Subcluster_T/Annotate/CD4/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(CD4,file = '3.1Subcluster_T/Annotate/CD4/CD4.rds')
####annotate celltype####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster0.pdf')
plot_density(CD4,features = c('CD55'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster1.pdf')
plot_density(CD4,features = c('ANXA1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster2.pdf')
plot_density(CD4,features = c('AIM1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster3.pdf')
plot_density(CD4,features = c('RTKN2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster4.pdf')
plot_density(CD4,features = c('FOS'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster5.pdf')
plot_density(CD4,features = c('HSPA1A'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster6.pdf')
plot_density(CD4,features = c('MAF'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster7.pdf')
plot_density(CD4,features = c('GIMAP5'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster8.pdf')
plot_density(CD4,features = c('PDCD1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster9.pdf')
plot_density(CD4,features = c('MIR155HG'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster10.pdf')
plot_density(CD4,features = c('EOMES'),reduction = 'umap')
dev.off()

####annotate
minor_celltype <- CD4@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "C0 CD4Tn-CD55", minor_celltype)
minor_celltype <- gsub("^1$", "C1 CD4Tm-ANXA1", minor_celltype)
minor_celltype <- gsub("^2$", "C2 CD4Tm-AIM1", minor_celltype)
minor_celltype <- gsub("^3$", "C3 CD4Treg-RTKN2", minor_celltype)
minor_celltype <- gsub("^4$", "C4 CD4Tm-FOS", minor_celltype)
minor_celltype <- gsub("^5$", "C5 CD4Tstd-HSPA1A", minor_celltype)
minor_celltype <- gsub("^6$", "C6 CD4Treg-MAF", minor_celltype)
minor_celltype <- gsub("^7$", "C7 CD4T-GIMAP5", minor_celltype)
minor_celltype <- gsub("^8$", "C8 CD4Tex-PDCD1", minor_celltype)
minor_celltype <- gsub("^9$", "C9 CD4Tex_pre-MIR155HG", minor_celltype)
minor_celltype <- gsub("^10$", "C10 CD4Tr1-EOMES", minor_celltype)


table(minor_celltype)
CD4 <- AddMetaData(CD4,minor_celltype,col.name = 'minor_celltype')
table(CD4$minor_celltype)
CD4$minor_celltype <- factor(CD4$minor_celltype,levels = c('C0 CD4Tn-CD55','C2 CD4Tm-AIM1','C1 CD4Tm-ANXA1','C4 CD4Tm-FOS',
                                                           'C3 CD4Treg-RTKN2','C6 CD4Treg-MAF',
                                                           'C8 CD4Tex-PDCD1','C10 CD4Tr1-EOMES','C9 CD4Tex_pre-MIR155HG','C5 CD4Tstd-HSPA1A','C7 CD4T-GIMAP5'))
Idents(CD4) <- CD4$minor_celltype
DefaultAssay(CD4) <- "RNA"

source('./function/do_dimplot.R') 
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/cellType.pdf",height = 15,width = 18)
DoDimplot(CD4,groupby = 'minor_celltype',colors = Palettes[['mycols_11']],pt.size = .8)
dev.off()

pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cellType.clean.pdf',width = 10,height = 8)
DimPlot(CD4,group.by = 'minor_celltype',reduction = 'umap',cols = Palettes[['mycols_11']],pt.size = .5)
DimPlot(CD4,group.by = 'minor_celltype',reduction = 'tsne',cols = Palettes[['mycols_11']],pt.size = .5)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = CD4, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = CD4, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_11']])
ratio.plot(seurat.object = CD4, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(CD4,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_11']])
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(CD4,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_11']],ncol = 2)
dev.off()

pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/C9.marker.pdf')
VlnPlot(CD4,features = c('MIR155HG','C1orf228','TSHZ2'),pt.size = 0,ncol = 2)
FeaturePlot(CD4,features = c('MIR155HG','C1orf228','TSHZ2'),order = T,reduction = 'umap',min.cutoff = 0,ncol = 2)
dev.off()

saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
####use plot1cell to visualize####
library(plot1cell)
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
table(Idents(CD4))
circ_data <- prepare_circlize_data(CD4, scale = 0.8 )
group_colors<-rand_color(length(names(table(CD4$group))))
rep_colors<-rand_color(length(names(table(CD4$orig.ident))))
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.3, col.use = Palettes[['mycols_11']] ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()
####expression of celltype markers####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
genelist <- c('CD55','CCR7','LEF1','TCF7','IL7R','CD40LG','CD69','AIM1','ANXA1','FOS','JUN','FOXP3','IL2RA','BATF','RTKN2','IKZF2','MAF',
              'CTLA4','TIGIT','PDCD1','TOX','TOX2','ITM2A','EOMES','HAVCR2','IL10','GZMK','GZMA','LGALS3','LYST','MIR155HG','LAG3','C1orf228','TSHZ2','REL','IRF4','ENO1','GBP2','HSPA1A','HSPA1B','HSPD1','HSPE1','DNAJB1','DNAJB4','GIMAP5','GIMAP4','GIMAP7','GIMAP1')
Idents(CD4) <- CD4$minor_celltype
exp.matrix <- GetAssayData(CD4, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, CD4$minor_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform standardization for all clusters per gene
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/celltype_expression.pdf',width = 10,height = 18)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(30, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Expression", at = c(0,1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/celltype_expression.dotplot.pdf',width = 18,height = 8)
Idents(CD4) <- CD4$minor_celltype
DotPlot(CD4,features = genelist,scale.by = 'size',col.min = -1.5) + scale_color_gradientn(colours = Palettes[['blueRed']]) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(CD4) <- CD4$minor_celltype
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
########test of celltype percentage between groups####
CD4 <- readRDS('3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
df <- as.data.frame(table(CD4$minor_celltype,CD4$sample))
df$total <- apply(df,1,function(x){sum(df$Freq[df$Var2==x[2]])})
df$percent <- round(df$Freq/df$total,4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype.percentage.group.pdf',height = 12,width = 15)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    theme_classic() +
    labs(title = i,x = 'Group',y = '% of CD4+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 4)
dev.off()

metadata <- CD4@meta.data
res <- propeller(clusters = metadata$minor_celltype,sample = metadata$sample,group = metadata$group)
plotCellTypeProps(clusters = metadata$minor_celltype,sample = metadata$sample)
write.csv(res,file = './3.1Subcluster_T/Annotate/CD4/celltype.propeller.test.csv')
saveRDS(CD4,file = '3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
