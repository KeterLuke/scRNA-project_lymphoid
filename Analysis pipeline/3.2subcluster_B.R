#' @description: annotate B subclusters

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
sub.scRNA <- subset(data.merge, subset = major_celltype == 'B')
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
dir.create('./3.2Subcluster_B')
sub.list.Standard <- variableFeatureSelection(seurat.lists = sub.list, method = "Standard", nfeatures = 3000)
saveRDS(sub.list.Standard, file = "3.2Subcluster_B/sub.list.Standard.3000.rds")

sub.list.SCT <- variableFeatureSelection(seurat.lists = sub.list, method = "SCT", nfeatures = 3000,return.only.var.genes = T,vars.to.regress = c("nCount_RNA", "percent.mt"))
saveRDS(sub.list.SCT, file = "3.2Subcluster_B/sub.list.SCT.3000.rds")

## Remove previous clustering results
sub.scRNA <- merge(sub.list.SCT[[1]], y = sub.list.SCT[2:length(sub.list.SCT)], project = "B cells")
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
pdf("3.2Subcluster_B/RNA.Harmony.PC25.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 25, npcs = 100)
dev.off()

####RNA.Harmony.PC30####
pdf("3.2Subcluster_B/RNA.Harmony.PC30.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 30, npcs = 100)
dev.off()

####RNA.Harmony.PC35####
pdf("3.2Subcluster_B/RNA.Harmony.PC35.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 35, npcs = 100)
dev.off()

##select RNA.Harmony.PC35
pdf("3.2Subcluster_B/RNA.Harmony.PC35.Integration.pdf")
sub.scRNA.harmony <- Harmony.integration.reduceDimension(seurat.object = sub.scRNA, assay = "RNA", set.resolutions = seq(0.1, 1, by = 0.1), PC = 35, npcs = 100)
dev.off()

####select resolution:0.3
sub.scRNA.harmony$seurat_clusters <- sub.scRNA.harmony$RNA_snn_res.0.3
table(sub.scRNA.harmony$seurat_clusters)
Idents(sub.scRNA.harmony) <- sub.scRNA.harmony$seurat_clusters
table(Idents(sub.scRNA.harmony))
table(sub.scRNA.harmony$celltypist_minor,sub.scRNA.harmony$seurat_clusters)
saveRDS(sub.scRNA.harmony, file = "3.2Subcluster_B/sub.scRNA.harmony.rds")

sub.B <- readRDS('./3.2Subcluster_B/sub.scRNA.harmony.rds')
length(unique(sub.B$seurat_clusters))
dir.create('./3.2Subcluster_B/Annotate')
pdf('./3.2Subcluster_B/Annotate/cluster.pdf')
DimPlot(sub.B,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']])
DimPlot(sub.B,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']])
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')

#### Differential expression####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
Idents(sub.B) <- sub.B$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.B, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

####remove cluster5(express T&B cell markers while not express cycling-related genes)
sub.B <- subset(sub.B,subset = seurat_clusters!='5')
sub.B <- NormalizeData(sub.B, verbose = FALSE)
sub.B <- FindVariableFeatures(sub.B,nfeatures = 3000)
sub.B <- ScaleData(sub.B, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(sub.B))
sub.B <- RunPCA(sub.B,verbose = F,npcs = 100)
options(timeout = 1600000000)
####perform integration and reduction####
sub.B  <- RunHarmony(sub.B ,reduction = 'pca',verbose = F,group.by.vars = 'orig.ident')
sub.B  <- RunUMAP(sub.B ,verbose = F,reduction = 'harmony',dims = 1:35)
sub.B  <- RunTSNE(sub.B ,verbose = F,reduction = 'harmony',dims = 1:35)
sub.B  <- FindNeighbors(sub.B ,verbose = F,reduction = 'harmony',dims = 1:35)
sub.B  <- FindClusters(sub.B ,resolution = seq(0.1, 1, by = 0.1),verbose = F)
pdf('./3.2Subcluster_B/RNA.Harmony.PC35.pro.PC35.pdf')
set.resolutions <- seq(0.1, 1, by = 0.1)
clustree(sub.B)
DimPlot(sub.B ,reduction = 'umap',group.by = 'orig.ident')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.B , reduction = 'umap',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
DimPlot(sub.B ,reduction = 'tsne',group.by = 'orig.ident')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sub.B , reduction = 'tsne',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
dev.off()
####select resolution:0.2
sub.B$seurat_clusters <- sub.B$RNA_snn_res.0.2
table(sub.B$seurat_clusters)
Idents(sub.B) <- sub.B$seurat_clusters
table(Idents(sub.B))
length(unique(sub.B$seurat_clusters))
cluster.pos.markers <- FindAllMarkers(sub.B, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf('./3.2Subcluster_B/Annotate/cluster.pdf')
DimPlot(sub.B,reduction = 'umap',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']])
DimPlot(sub.B,reduction = 'tsne',group.by = 'seurat_clusters',label = T,cols = Palettes[['mycols_8']])
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
#### classical marker#### 
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
pdf('./3.2Subcluster_B/Annotate/classical.marker.B.pdf')
FeaturePlot(sub.B,features  = c('CD79A','CD79B','CD19'),reduction = 'umap',order = T,ncol = 2)
FeaturePlot(sub.B,features  = c('MS4A1'),reduction = 'umap',order = F,pt.size = .3)
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.Naive.pdf')
plot_density(sub.B,features  = c('MS4A1','IGHD','TCL1A','IL4R'),reduction = 'umap',joint = F)##Naive
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.Memory.pdf')
plot_density(sub.B,features  = c('MS4A1','CD27','AIM2','TNFRSF13B'),reduction = 'umap',joint = F)##Memory
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.Plasma.pdf')
plot_density(sub.B,features  = c('JCHAIN','SDC1','IGHA1','IGHG1'),reduction = 'umap',joint = F)##Plasma
dev.off()
pdf('./3.2Subcluster_B/Annotate/classical.marker.GC.pdf')
plot_density(sub.B,features  = c('AICDA','LRMP','RGS13','SUGCT'),reduction = 'umap',joint = F)##Germinal center B
dev.off()

##cluster similarity#####
expMatrix <- GetAssayData(sub.B, slot = "scale.data")
highVariableGenes <- VariableFeatures(sub.B)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, sub.B$seurat_clusters, mean)
  return(mean.value)
})

pdf('3.2Subcluster_B/Annotate/cor_cluster.pdf')
cor <- cor(t(meanExpCluster), method="spearman")
cor <- (cor - min(cor))/(max(cor)-min(cor))
pheatmap::pheatmap(cor,angle_col = 45,clustering_method = 'complete')
dev.off()

corrMatrix <- (1- cor(t(meanExpCluster), method="spearman"))/2

library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("3.2Subcluster_B/Annotate/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
###############celltypist predict########
library(reticulate)
py_config()
scanpy <- import('scanpy')
pandas <- import('pandas')
numpy <- import('numpy')
celltypist <- import('celltypist')
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(sub.B[['RNA']]@counts)))),
                       obs = pandas$DataFrame(sub.B@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(sub.B[['RNA']]@counts),
                                                         row.names = rownames(sub.B[['RNA']]@counts)))
)
adata
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)
predicted_minor <- as.data.frame(predictions$predicted_labels)
predicted_minor$seurat_cluster <- sub.B@meta.data$seurat_clusters
sub.B <- AddMetaData(sub.B,predicted_minor$majority_voting,col.name = 'celltypist_minor')
table(sub.B$celltypist_minor,sub.B$seurat_clusters)

write.xlsx(predicted_minor,file = '3.2Subcluster_B/Annotate/celltypist_predicted.minor.xlsx',rowNames = T)
saveRDS(sub.B,file = '3.2Subcluster_B/sub.B.rds')

####expression of functional markers####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
cell.type.markers <- read.table(file = "./3.2Subcluster_B/Annotate/B_functional markers.txt", header = T, stringsAsFactors = F, sep = "\t")
exp.matrix <- GetAssayData(sub.B, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.B$seurat_clusters, mean)
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
pdf("3.2Subcluster_B/Annotate/cluster.functional.signature.pdf",height = 10)
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
#### Differential expression####
Idents(sub.B) <- sub.B$seurat_clusters
cluster.pos.markers <- FindAllMarkers(sub.B, only.pos = TRUE, group.by = "seurat_clusters",
                                      min.diff.pct = 0.1,logfc.threshold = 0.1,test.use = 'MAST',latent.vars = 'orig.ident')
cluster.sig.markers <- cluster.pos.markers[which(cluster.pos.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top20 <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(cluster.sig.markers.top20, file = "3.2Subcluster_B/Annotate/cluster.sig.markers.top20.txt", col.names = T, row.names = F, sep = "\t", quote = F)

top5.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sub.B <- ScaleData(sub.B,features = c(VariableFeatures(sub.B),unique(top5.genes$gene)))
pdf("3.2Subcluster_B/Annotate/cluster.top5genes.pdf",width = 15,height = 15)
DoHeatmap(sub.B, features = unique(top5.genes$gene), size = 2,group.colors = Palettes[['mycols_12']],group.bar = T) + NoLegend() + 
  scale_fill_gradientn(colors = Palettes[['greenBlue']])
dev.off()

top8.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
sub.B <- ScaleData(sub.B,features = c(VariableFeatures(sub.B),unique(top8.genes$gene)))
source(('./function/do_heatmap.R'))
pdf('./3.2Subcluster_B/Annotate/cluster.top8genes.mean.pdf',height = 12)
DoHeatmap_my(sub.B,features = unique(top8.genes$gene),slot = 'scale.data',assay = 'RNA',cluster_cols = F,group.by = 'seurat_clusters',
             color = colorRampPalette(Palettes[['blueRed']])(100))
dev.off()
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
####GSVA####
dir.create('./3.2Subcluster_B/Annotate/enrichment')
dir.create('./3.2Subcluster_B/Annotate/enrichment/GSVA')
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
Idents(sub.B) <- sub.B$seurat_clusters
av <- AverageExpression(sub.B,assays = 'RNA',group.by = 'seurat_clusters',slot = 'counts')
av <- av[[1]]
head(av)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(GSVA)
####gsvago
# gene_sets <- msigdbr(species = 'Homo sapiens',category = 'C5',subcategory = 'GO:BP')
gene_sets <- data.table::fread('./3.1Subcluster_T/Annotate/enrichment/immune_related.txt')
# rn <- gene_sets$gs_name
# rn <- unlist(lapply(rn,function(i){
#   i <- gsub('GOBP_',"",i)
#   i <- gsub("_"," ",i)
#   return(i)
# }))
# rn <- tolower(rn)
# gene_sets$gs_name <- rn

gs=split(gene_sets$`Gene name`,gene_sets$`GO term name`)
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
top <- df%>% group_by(cluster) %>% top_n(8,fc)
pdf('./3.2Subcluster_B/Annotate/enrichment/GSVA/gsvaGo.pdf',width = 8)
pheatmap(es[unique(top$go),],show_rownames = T,clustering_method = 'complete',fontsize = 5,fontsize_col = 8,angle_col = '45',main = 'GSVA_clusters')
dev.off()
saveRDS(es.max,file = './3.2Subcluster_B/Annotate/enrichment/GSVA/gsvaGO.rds')
idents <- levels(sub.B)
df$cluster <- factor(df$cluster,levels = levels(sub.B))
res <- split(df,df$cluster)
saveFormat <- lapply(res,function(i){
  i <- i %>% dplyr::filter(fc > 0) %>% dplyr::arrange(desc(fc))
  return(i)
})
write.xlsx(saveFormat,file = './3.2Subcluster_B/Annotate/enrichment/GSVA/gsvaGo.xlsx',sheetName = idents)
#####enrichgo####
cluster.DE <- FindAllMarkers(sub.B,only.pos = F,group.by = 'seurat_clusters',min.diff.pct = 0.1,logfc.threshold = 0.1,
                             test.use = 'MAST',latent.vars = 'orig.ident')
idents <- levels(sub.B)
saveFormat <- lapply(idents, function(x){
  index <- which(cluster.DE$cluster == x)
  DEGs <- cluster.DE[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "3.2Subcluster_B/Annotate/cluster.all.DE.xlsx", sheetName = idents, rowNames = F)
saveRDS(cluster.DE, file = "3.2Subcluster_B/Annotate/cluster.all.DE.rds")

DEGs <- lapply(saveFormat, function(x){
  x <- x %>% filter(p_val_adj<0.05 & avg_log2FC>0.4) %>% arrange(desc(avg_log2FC))
  return(x)
})
names(DEGs) <- idents
write.xlsx(DEGs, file = "3.2Subcluster_B/Annotate/clusters.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(DEGs, file = "3.2Subcluster_B/Annotate/clusters.DEGs.rds")
idents <- levels(sub.B)
DEGs <- lapply(DEGs, function(x){
  return(x$gene)
})
source("function/clusterProfiler.enricher.R")

####enrichgo
dir.create('./3.2Subcluster_B/Annotate/enrichment/GO')
pdf('./3.2Subcluster_B/Annotate/enrichment/GO/cluster.clusterProfiler.enrichgo.pdf')
res <- lapply(names(DEGs), function(x){
  y <- DEGs[[x]]
  res <- cluterProfiler.enricher(gene = y, geneType = "SYMBOL", db.type = "GO",GO.ont = 'BP',filename = 'enrichgoResult.csv',
                                 saveDir = paste0(getwd(),'/3.2Subcluster_B/Annotate/enrichment/GO/'),
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
write.xlsx(res, file = "3.2Subcluster_B/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx", sheetName = idents, rowNames = F)
saveRDS(res,file = '3.2Subcluster_B/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.rds')
#### show cluster profiler result####
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}
res <- lapply(1:6,function(i){
  pro <- read.xlsx('./3.2Subcluster_B/Annotate/enrichment/GO/clusterProfiler.enrichgo.result.xlsx',sheet = i)
  return(pro)
})
names(res) <- 0:5
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
enrich.res <- enrich.res %>% group_by(Cluster) %>% top_n(10,wt = p.adjust)
enrich.res$wrap <- wrapText(enrich.res$Description, 45)
enrich.res$Cluster <- factor(enrich.res$Cluster,levels = c(0:5),ordered = T)
enrich.res$p.adjust[enrich.res$p.adjust >= 12] <- 12
pdf("3.2Subcluster_B/Annotate/enrichment/GO/cluterProfiler_go.enrichment.pdf",height = 15)
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
####add module score####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
Idents(sub.B) <- sub.B$seurat_clusters
signature <- list(Naive = c('YBX3','IGHD','TCL1A','FCER2','IL4R'),
                  Memory = c('CD27','TNFRSF13B','AIM2'),
                  Plasma = c('SDC1','MZB1','JCHAIN','XBP1','IGHG1','IGHA1','IGHG4'),
                  GC = c('AICDA','MKI67','KIAA0101','LRMP','SUGCT'))
sce <- AddModuleScore(sub.B,features = signature,name = names(signature))
pdf('./3.2Subcluster_B/Annotate/Signature_score.pdf')
VlnPlot(sce,features = c('Naive1','Memory2','Plasma3','GC4'),ncol = 2,pt.size = 0,cols = Palettes[['mycols_8']])
FeaturePlot(sce,features = c('Naive1','Memory2','Plasma3','GC4'),reduction = 'umap',ncol = 2,order = T,min.cutoff = 0,cols = Palettes[['greyMagma']],pt.size = 0.01)
dev.off()
sub.B <- sce
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.rds')
####GSEA####
dir.create('./3.2Subcluster_B/Annotate/enrichment/GSEA')
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
gmtfile = './3.1Subcluster_T/Annotate/enrichment/GSEA/c5.go.bp.v2022.1.Hs.symbols.gmt'
geneset <- read.gmt(gmtfile)
deg <- readRDS('./3.2Subcluster_B/Annotate/cluster.all.DE.rds')
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
saveRDS(res,file = './3.2Subcluster_B/Annotate/enrichment/GSEA/GSEA_result.rds')
idents <- names(res)
write.xlsx(res,file = './3.2Subcluster_B/Annotate/enrichment/GSEA/GSEA_result.xlsx',sheetName = idents)
res_sig <- lapply(res,function(i){
  res <- i %>% dplyr::filter(pvalue < 0.05&p.adjust < 0.25)
  return(res)
})
write.xlsx(res_sig,file = './3.2Subcluster_B/Annotate/enrichment/GSEA/GSEA_result.sig.xlsx',sheetName = idents)
####annotate celltype####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.rds')
pdf('./3.2Subcluster_B/Annotate/CD27.pdf')
FeaturePlot(sub.B,features = 'CD27',reduction = 'umap',order = T,cols = Palettes[['greyMagma']],pt.size = .3)
dev.off()
pdf('./3.2Subcluster_B/Annotate/TCL1A.pdf')
FeaturePlot(sub.B,features = 'TCL1A',reduction = 'umap',order = T,cols = Palettes[['greyMagma']],pt.size = .3,min.cutoff = .5)
dev.off()
pdf('./3.2Subcluster_B/Annotate/IGHD.pdf')
FeaturePlot(sub.B,features = 'IGHD',reduction = 'umap',order = T,cols = Palettes[['greyMagma']],pt.size = .3,min.cutoff = .5)
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster0.pdf')
plot_density(sub.B,features = c('TNFRSF13B'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster1.pdf')
plot_density(sub.B,features = c('IGHD'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster2.pdf')
plot_density(sub.B,features = c('FOSB'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster3.pdf')
plot_density(sub.B,features = c('HSPA1A'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster4.pdf')
plot_density(sub.B,features = c('MZB1'),reduction = 'umap')
dev.off()
pdf('./3.2Subcluster_B/Annotate/cluster5.pdf')
plot_density(sub.B,features = c('LRMP'),reduction = 'umap')
dev.off()

####annotate
minor_celltype <- sub.B@meta.data$seurat_clusters
minor_celltype <- gsub("^0$", "Bmem_TNFRSF13B", minor_celltype)
minor_celltype <- gsub("^1$", "Bn_IGHD", minor_celltype)
minor_celltype <- gsub("^2$", "Bmem_FOSB", minor_celltype)
minor_celltype <- gsub("^3$", "Bstd_HSPA1A", minor_celltype)
minor_celltype <- gsub("^4$", "Plasma_MZB1", minor_celltype)
minor_celltype <- gsub("^5$", "Bgc_LRMP", minor_celltype)

table(minor_celltype)
sub.B <- AddMetaData(sub.B,minor_celltype,col.name = 'minor_celltype')
table(sub.B$minor_celltype)
sub.B$minor_celltype <- factor(sub.B$minor_celltype,levels = names(table(sub.B$minor_celltype)))
Idents(sub.B) <- sub.B$minor_celltype
DefaultAssay(sub.B) <- "RNA"

source('./function/do_dimplot.R')
pdf("3.2Subcluster_B/Annotate/cellType.pdf",height = 10,width = 12)
DoDimplot(sub.B,groupby = 'minor_celltype',colors = Palettes[['mycols_6']])
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_ratio.pdf',height = 12)
ratio.plot(seurat.object = sub.B, id.vars1 = "sample", id.vars2 = "minor_celltype", angle = 60)
ratio.plot(seurat.object = sub.B, id.vars1 = "minor_celltype", id.vars2 = "sample", angle = 60,color.len = Palettes[['mycols_6']])
ratio.plot(seurat.object = sub.B, id.vars1 = "group", id.vars2 = "minor_celltype", angle = 60)
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_group.dimplot.pdf',width = 12,height = 5)
DimPlot(sub.B,reduction = 'umap',group.by = 'minor_celltype',split.by = 'group',cols = Palettes[['mycols_6']])
dev.off()

pdf('./3.2Subcluster_B/Annotate/celltype_sample.dimplot.pdf',width = 10,height = 12)
DimPlot(sub.B,reduction = 'umap',group.by = 'minor_celltype',split.by = 'sample',cols = Palettes[['mycols_6']],ncol = 2)
dev.off()

saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.pro.rds')
####expression of celltype markers####
sub.B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
genelist <- c('IGHD','TCL1A','FCER2','IL4R','YBX3','CD27','TNFRSF13B','AIM2','VIM','FOSB','FOS','IER2','HSPA1A',
              'HSPA1B','DNAJB1','HSPE1','HSPH1','HSPD1','MZB1','SDC1','JCHAIN','XBP1','IGHG1','IGHA1','IGHG4',
              'LRMP','RGS13','AICDA','MKI67','KIAA0101','SUGCT')
sub.B <- ScaleData(sub.B,features = c(VariableFeatures(sub.B),genelist))
sub.B@meta.data$celltype <- sub.B@meta.data$minor_celltype
sub.B@meta.data$celltype <- factor(sub.B@meta.data$celltype,levels = c('Bn_IGHD','Bmem_TNFRSF13B','Bmem_FOSB',
                                                                       'Bstd_HSPA1A','Plasma_MZB1','Bgc_LRMP'),ordered = T)
Idents(sub.B) <- sub.B$celltype
exp.matrix <- GetAssayData(sub.B, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, sub.B$celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
pdf('./3.2Subcluster_B/Annotate/celltype_expression.pdf',width = 12,height = 22)
Heatmap(t(cluster.score.normailzed), 
        width = unit(15, "cm"), height = unit(24, "cm"),column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

Idents(sub.B) <- sub.B$celltype
pdf('./3.2Subcluster_B/Annotate/celltype_expression.dotplot.pdf',width = 10)
DotPlot(sub.B,features = genelist,cols = c('grey','red'),scale.by = 'size',col.min = -0.4,scale.min = -5) + 
  theme(axis.text.x = element_text(size = 8)) + coord_flip()
dev.off()

Idents(sub.B) <- sub.B$minor_celltype
saveRDS(sub.B,file = './3.2Subcluster_B/sub.B.pro.rds')
