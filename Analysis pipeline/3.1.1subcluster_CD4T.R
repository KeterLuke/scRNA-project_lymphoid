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
####RNA.Harmony.PC40####
pdf("3.1Subcluster_T/Annotate/CD4/RNA.Harmony.PC40.Integration.pdf")
CD4.harmony <- Harmony.integration.reduceDimension(seurat.object = CD4, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 40, npcs = 100)
dev.off()
saveRDS(CD4.harmony,file = './3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC40.rds')
##select RNA.Harmony.PC15 res:0.4####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/RNA.Harmony.Integration.PC15.rds')
CD4$seurat_clusters <- CD4$RNA_snn_res.0.4
length(unique(CD4$seurat_clusters))
table(CD4$seurat_clusters)
dir.create('./3.1Subcluster_T/Annotate/CD4/Annotate')
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster.pdf')
DimPlot(CD4,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_9']])
DimPlot(CD4,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_9']])
dev.off()
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.rds')
#### Differential expression####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
Idents(CD4) <- CD4$seurat_clusters
table(Idents(CD4))
cluster.pos.markers <- FindAllMarkers(CD4, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
CD4 <- ScaleData(CD4,features = c(VariableFeatures(CD4),unique(top5.genes$gene)))
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(CD4, features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['mycols_9']],group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
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
signature <- list(Naive = c('CCR7','LEF1','CD55','TCF7'),
                  Memory = c('IL7R','FOS','CD40LG','JUN','ANXA1','S100A4'),
                  Treg = c('FOXP3','IL2RA','BATF','IKZF2','CTLA4'),
                  Exhaustion = c('PDCD1','TIGIT','LAG3','HAVCR2','TOX2','MAF','ITM2A','LGALS3'))
sce <- AddModuleScore(CD4,features = signature,name = names(signature))
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/Signature_score.pdf',width = 8,height = 12)
VlnPlot(sce,features = c('Naive1','Memory2','Treg3','Exhaustion4'),ncol = 2,pt.size = 0,cols = Palettes[['mycols_9']])
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
####enrichment####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
cluster.DE <- FindAllMarkers(CD4,only.pos = F,group.by = 'seurat_clusters',min.diff.pct = 0.1,test.use = 'MAST',
                             latent.vars = 'orig.ident')
idents <- levels(CD4)
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.all.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "3.1Subcluster_T/Annotate/CD4/Annotate/cluster.all.DE.rds")

##-- Functional enrichment analysis
DEGs <- lapply(saveFormat, function(x){
  x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.25) %>% arrange(desc(avg_log2FC))
  return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "3.1Subcluster_T/Annotate/CD4/Annotate/clusters.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "3.1Subcluster_T/Annotate/CD4/Annotate/clusters.DEGs.rds")
####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
DEGs <- readRDS('3.1Subcluster_T/Annotate/CD4/Annotate/clusters.DEGs.rds')
idents <- levels(CD4)
DEGs <- lapply(DEGs, function(x){
  return(x$gene)
})
source("function/clusterProfiler.enricher.R")
####enricher
dir.create('3.1Subcluster_T/Annotate/CD4/Annotate/enrichment')
dir.create('3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/msigDB')
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/msigDB/cluster.clusterProfiler.enricher.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "MsigDB",saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/msigDB'),
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
write.xlsx(res, file = "3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/msigDB/clusterProfiler.enricher.result.rds')
####enrichgo
dir.create('3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO')
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/cluster.clusterProfiler.enrichgo.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "GO",GO.ont = 'BP',filename = 'enrichgoResult.csv',
                                 saveDir = paste0(getwd(),'/3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/'),
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
write.xlsx(res, file = "3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.rds')
####show clusterprofiler result####
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
res <-  read.xlsx('./3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx',sheet = 6)
library(ggthemes)
enrich <- res[1:10, c("Description", "p.adjust", "Count",'GeneRatio')]
rownames(enrich) <- NULL

enrich$p.adjust <- -log10(enrich$p.adjust)
enrich$wrap <- wrapText(enrich$Description, 45)
enrich$p.adjust[enrich$p.adjust >= 12] <- 12
enrich <- enrich %>% arrange(Count)
enrich$wrap <- factor(enrich$wrap,ordered = T,levels = enrich$wrap)
rownames(enrich) <- NULL
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/enrichment/GO/cluterProfiler_go.enrichment_cluster5.pdf")
ggplot(data = enrich,
                aes(y = wrap,
                x = Count, 
                fill = p.adjust)) + geom_bar(stat = 'identity',position = 'dodge') + theme_classic() +
  scale_fill_gradient(low = "grey", high = "red",breaks=c(2,5,8,12),limits = c(1.3,14)) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1, size = 8), 
                 axis.text.y = element_text(size = 5)) + xlab("Count") + ylab("") + 
  labs(fill='-log10(p.adjust)',title = 'cluster 5')
dev.off()
####annotate celltype####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.rds')
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster0.pdf')
plot_density(CD4,features = c('CCR7'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster1.pdf')
plot_density(CD4,features = c('IL2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster2.pdf')
plot_density(CD4,features = c('ANXA1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster3.pdf')
plot_density(CD4,features = c('IKZF2'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster4.pdf')
plot_density(CD4,features = c('GIMAP1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster5.pdf')
plot_density(CD4,features = c('PDCD1'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster6.pdf')
plot_density(CD4,features = c('HSPA1A'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster7.pdf')
plot_density(CD4,features = c('IL2RA'),reduction = 'umap')
dev.off()
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cluster8.pdf')
plot_density(CD4,features = c('IFI6'),reduction = 'umap')
dev.off()
####annotate
minor_celltype <- CD4@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "C0 CD4Tn-CCR7", minor_celltype)
minor_celltype <- gsub("^1$", "C1 CD4Tm-IL2", minor_celltype)
minor_celltype <- gsub("^2$", "C2 CD4Tm-ANXA1", minor_celltype)
minor_celltype <- gsub("^3$", "C3 CD4Treg-IKZF2", minor_celltype)
minor_celltype <- gsub("^4$", "C4 CD4Tm-GIMAP1", minor_celltype)
minor_celltype <- gsub("^5$", "C5 CD4Tex-PDCD1", minor_celltype)
minor_celltype <- gsub("^6$", "C6 CD4Tstd_HSPA1A", minor_celltype)
minor_celltype <- gsub("^7$", "C7 CD4Treg-IL2RA", minor_celltype)
minor_celltype <- gsub("^8$", "C8 CD4Tm-IFI6", minor_celltype)

table(minor_celltype)
CD4 <- AddMetaData(CD4,minor_celltype,col.name = 'minor_celltype')
table(CD4$minor_celltype)
CD4$minor_celltype <- factor(CD4$minor_celltype,levels = c('C0 CD4Tn-CCR7','C1 CD4Tm-IL2','C2 CD4Tm-ANXA1','C4 CD4Tm-GIMAP1',
                                                           'C8 CD4Tm-IFI6','C3 CD4Treg-IKZF2','C7 CD4Treg-IL2RA',
                                                           'C5 CD4Tex-PDCD1','C6 CD4Tstd_HSPA1A'))
Idents(CD4) <- CD4$minor_celltype
DefaultAssay(CD4) <- "RNA"

source('./function/do_dimplot.R') 
pdf("3.1Subcluster_T/Annotate/CD4/Annotate/cellType.pdf",height = 10,width = 12)
DoDimplot(CD4,groupby = 'minor_celltype',colors = Palettes[['mycols_9']],pt.size = .5)
dev.off()

pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/cellType.clean.pdf',width = 10,height = 8)
DimPlot(CD4,group.by = 'minor_celltype',reduction = 'umap',cols = Palettes[['mycols_9']],pt.size = .5)
DimPlot(CD4,group.by = 'minor_celltype',reduction = 'tsne',cols = Palettes[['mycols_9']],pt.size = .5)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = CD4, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = CD4, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_9']])
ratio.plot(seurat.object = CD4, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(CD4,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_9']])
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(CD4,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_9']],ncol = 2)
dev.off()

saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
####expression of celltype markers####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
genelist <- c('CD55','CCR7','LEF1','TCF7','IL2','CD40LG','IL7R','FOS','JUN','ANXA1','S100A4','GIMAP1','GIMAP4',
              'GIMAP7','IFI6','IFI44L','IFI44','IFIT1','IFIT3','FOXP3','IKZF2','RTKN2','IL2RA','BATF','CTLA4','TIGIT','MAF','PDCD1',
              'ICA1','TOX2','HSPA1A','HSPA1B','HSPD1','HSPE1','DNAJB1','DNAJB4')
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
DotPlot(CD4,features = genelist,cols = c('grey','red'),scale.by = 'size',col.min = 0,scale.min = -5) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(CD4) <- CD4$minor_celltype
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
########test of celltype percentage between groups####
CD4 <- readRDS('3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
df <- as.data.frame(table(CD4$minor_celltype,CD4$sample))
df$percent <- round(df$Freq/ncol(CD4),4) * 100
df$group <- ifelse(grepl(df$Var2,pattern = 'T'),'tumor','normal')
idents <- unique(df$Var1)
pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype.percentage.group.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggpaired(tmp,x = 'group',y = 'percent',fill = 'group',palette = 'npg',width = .5,line.size = .1,line.color = 'grey') +
    stat_compare_means(comparisons = comparison,method = 'wilcox.test',paired = T) +
    theme_classic() +
    labs(title = i,x = 'Group',y = '% of CD4+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype.percentage.group.unpaired.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggplot(tmp,aes(x = group,y = percent,fill = group)) + geom_boxplot(width = .4) + geom_point() + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of CD4+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

pdf('3.1Subcluster_T/Annotate/CD4/Annotate/celltype.percentage.group.barplot.pdf',height = 18)
res <- lapply(idents,function(i){
  tmp <- df[df$Var1==i,]
  comparison = list(c('normal','tumor'))
  res <- ggbarplot(tmp,x = 'group',y = 'percent',fill = 'group',add = 'mean_se') + 
    stat_compare_means(comparisons = comparison,method = 'wilcox.test') +
    theme_classic() + ggsci::scale_fill_npg() + 
    labs(title = i,x = 'Group',y = '% of CD4+ T cells',fill = 'Group')
  return(res)
})
plot_grid(plotlist = res,ncol = 2)
dev.off()

saveRDS(CD4,file = '3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
####observe immune-exhaustion score in 3 Treg clusters####
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
Idents(CD4) <- CD4$minor_celltype
df <- FetchData(CD4,vars = c('Exhaustion4','minor_celltype'),cells = WhichCells(CD4,idents = c('C3 CD4Treg-IKZF2','C7 CD4Treg-IL2RA','C5 CD4Tex-PDCD1')))
table(df$minor_celltype)
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/Immune.exhaustion.score.compare.pdf')
ggviolin(df,x = 'minor_celltype',y = 'Exhaustion4',fill = 'minor_celltype',palette = 'npg',add = 'boxplot') + 
  stat_compare_means(comparisons = list(c('C3 CD4Treg-IKZF2','C7 CD4Treg-IL2RA'),
                                        c('C7 CD4Treg-IL2RA','C5 CD4Tex-PDCD1'),
                                        c('C3 CD4Treg-IKZF2','C5 CD4Tex-PDCD1')),
                     method = 'wilcox.test',label = 'p.signif') + ylab('Immune exhaustion score') + 
  stat_summary(fun.data = function(x) data.frame(y=-.8, label = paste("Mean=",round(mean(x),3))), geom="text")
dev.off()

df2 <- FetchData(CD4,vars = c('Exhaustion4','minor_celltype','group'),cells = WhichCells(CD4,idents = c('C3 CD4Treg-IKZF2','C7 CD4Treg-IL2RA','C5 CD4Tex-PDCD1')))
pdf('./3.1Subcluster_T/Annotate/CD4/Annotate/Immune.exhaustion.score.group.compare.pdf',width = 10)
ggviolin(df2,x = 'group',y = 'Exhaustion4',fill = 'minor_celltype',palette = 'npg',add = 'boxplot',facet.by = 'minor_celltype') + 
  stat_compare_means(comparisons = list(c('normal','tumor')),
                     method = 'wilcox.test',label = 'p.signif') + ylab('Immune exhaustion score') + 
  stat_summary(fun.data = function(x) data.frame(y=-.8, label = paste("Mean=",round(mean(x),3))), geom="text")
dev.off()
saveRDS(CD4,file = './3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
