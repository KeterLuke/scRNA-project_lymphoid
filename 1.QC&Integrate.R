#' @description: quality control,remove doublets,remove batch effect and integrate samples
# this work_lympho 6 samples
#### load package (always load all the package before you run the code)####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(circlize)
library(vegan)
library(ggplot2)
library(future)
set.seed(101)
setwd('../')
getwd()
options(future.globals.maxSize = 80000*1024^2)
plan("multisession",workers = 2)
####read samples####
sce <- Read10X('./raw/filtered_feature_bc_matrix/')

####Preliminary filtration####
#min.cell >3 & min.features >200   
sce <- CreateSeuratObject(counts = sce,
                   min.cells = 3,
                   min.features = 200,
                   project = 'lympho')
index <- gsub(pattern = '[A-Z]*-',replacement = "",colnames(sce))
head(index)
table(index)
sce$sample <- ifelse(index=='1','C4_6',
                     ifelse(index=='2','N6_15',
                            ifelse(index=='3','N6_8',
                                   ifelse(index=='4','T4_6',
                                          ifelse(index=='5','T6_15',
                                                 ifelse(index=='6','T6_8',NA))))))
table(sce$sample)
sce$orig.ident <- sce$sample
table(sce$orig.ident)
View(sce@meta.data)

## Mitochondrial gene ratio
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
dir.create('./1.QualityControl')
pdf(file = "1.QualityControl/count.feature.mt.pdf",width = 10)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1

ggdensity(sce@meta.data, x = "nCount_RNA", title = "nCount")
ggdensity(sce@meta.data, x = "nFeature_RNA", title = "nFeature")
dev.off()

#### Detect the resolution parameters by each sample cluster ####
##After the parameters are determined, you can block them without executing [test]
sce_obj <- subset(sce, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 500 & nFeature_RNA < 6000) 
sce_obj <- SplitObject(sce_obj,split.by = 'orig.ident')
##C4_6
C4_6 <- sce_obj[['C4_6']]
saveRDS(C4_6,file = './C4_6_afterqc.RDS')
##T4_6
T4_6 <- sce_obj[['T4_6']]
saveRDS(T4_6,file = './T4_6_afterqc.RDS')
##N6_15
N6_15 <- sce_obj[['N6_15']]
saveRDS(N6_15,file = './N6_15_afterqc.RDS')
##T6_15
T6_15 <- sce_obj[['T6_15']]
saveRDS(T6_15,file = './T6_15_afterqc.RDS')
##N6_8
N6_8 <- sce_obj[['N6_8']]
saveRDS(N6_8,file = './N6_8_afterqc.RDS')
##T6_8
T6_8 <- sce_obj[['T6_8']]
saveRDS(T6_8,file = './T6_8_afterqc.RDS')

rm(list = ls())
source('./function/colorPalettes.R')
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
set.resolutions <- seq(0.1, 1, by = 0.1)

####C4_6 8049 samples
sce_obj <- readRDS('./C4_6_afterqc.RDS')
pdf(file = "1.QualityControl/C4_6_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    print(p)
  })

dev.off()
saveRDS(sce_obj,file = './C4_6_afterqc.RDS')

####T4_6 6731 samples
sce_obj <- readRDS('./T4_6_afterqc.RDS')
pdf(file = "1.QualityControl/T4_6_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T4_6_afterqc.RDS')

####N6_8 15174 samples
sce_obj <- readRDS('./N6_8_afterqc.RDS')
pdf(file = "1.QualityControl/N6_8_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './N6_8_afterqc.RDS')

####T6_8 11194 samples
sce_obj <- readRDS('./T6_8_afterqc.RDS')
pdf(file = "1.QualityControl/T6_8_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T6_8_afterqc.RDS')


####N6_15 8597 samples
sce_obj <- readRDS('./N6_15_afterqc.RDS')
pdf(file = "1.QualityControl/N6_15_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './N6_15_afterqc.RDS')

####T6_15 6661 samples
sce_obj <- readRDS('./T6_15_afterqc.RDS')
pdf(file = "1.QualityControl/T6_15_PCA-test.pdf",width =15)
sce_obj <- SCTransform(sce_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:40, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:40)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T6_15_afterqc.RDS')
rm(list = ls())


#### remove doublet####
library(DoubletFinder) # Require cleanup of low-quality cells in advance
source('./function/doubletDetect.R')

####C4_6 8049 samples
sce <- readRDS('./C4_6_afterqc.RDS')
pdf("1.QualityControl/doublet_high/C4_6_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.1", sct = T)
dev.off()
##select:rate = 0.064
saveRDS(sce,file = './doublet_high/C4_6_afterqc.RDS')

####T4_6 6731 samples
sce <- readRDS('./T4_6_afterqc.RDS')
pdf("1.QualityControl/doublet_high/T4_6_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.3", sct = T)
dev.off()
##select:rate = 0.054
saveRDS(sce,file = './doublet_high/T4_6_afterqc.RDS')

####N6_8 15174 samples
sce <- readRDS('./N6_8_afterqc.RDS')
pdf("1.QualityControl/doublet_high/N6_8_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.2", sct = T)
dev.off()
##select:rate = 0.121
saveRDS(sce,file = './doublet_high/N6_8_afterqc.RDS')

####T6_8 11194 samples
sce <- readRDS('./T6_8_afterqc.RDS')
pdf("1.QualityControl/doublet_high/T6_8_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.1", sct = T)
dev.off()
##select:rate = 0.089
saveRDS(sce,file = './doublet_high/T6_8_afterqc.RDS')


####N6_15 8597 samples
sce <- readRDS('./N6_15_afterqc.RDS')
pdf("1.QualityControl/doublet_high/N6_15_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.5", sct = T)
dev.off()
##select:rate = 0.069
saveRDS(sce,file = './doublet_high/N6_15_afterqc.RDS')

####T6_15 6661 samples
sce <- readRDS('./T6_15_afterqc.RDS')
pdf("1.QualityControl/doublet_high/T6_15_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:40, doublet.rate = ncol(sce)*8*1e-6, annotation = "SCT_snn_res.0.5", sct = T)
dev.off()
##select:rate = 0.053
saveRDS(sce,file = './doublet_high/T6_15_afterqc.RDS')

#####subset singlet####
rm(list = ls())
#C4_6
sce <- readRDS('./doublet_high/C4_6_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "./doublet_high/C4_6.single.RDS")

#T4_6
sce <- readRDS('./doublet_high/T4_6_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet_high/T4_6.single.RDS")

#N6_8
sce <- readRDS('./doublet_high/N6_8_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet_high/N6_8.single.RDS")

#T6_8
sce <- readRDS('./doublet_high/T6_8_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet_high/T6_8.single.RDS")

#N6_15
sce <- readRDS('./doublet_high/N6_15_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet_high/N6_15.single.RDS")

#T6_15
sce <- readRDS('./doublet_high/T6_15_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet_high/T6_15.single.RDS")

#### merge data and correct the batch effect####
rm(list = ls())
C4_6 <- readRDS('./doublet_high/C4_6.single.RDS')
DefaultAssay(C4_6) <- "RNA"
T4_6 <- readRDS('./doublet_high/T4_6.single.RDS')
DefaultAssay(T4_6) <- "RNA"
N6_8 <- readRDS('./doublet_high/N6_8.single.RDS')
DefaultAssay(N6_8) <- "RNA"
T6_8 <- readRDS('./doublet_high/T6_8.single.RDS')
DefaultAssay(T6_8) <- "RNA"
N6_15 <- readRDS('./doublet_high/N6_15.single.RDS')
DefaultAssay(N6_15) <- "RNA"
T6_15 <- readRDS('./doublet_high/T6_15.single.RDS')
DefaultAssay(T6_15) <- "RNA"


source('./function/variableFeatureSelection.R')
scelist <- list(C4_6 = C4_6,
                T4_6 = T4_6,
                N6_8 = N6_8,
                T6_8 = T6_8,
                N6_15 = N6_15,
                T6_15 = T6_15)
rm(list = setdiff(ls(),c('scelist','variableFeatureSelection')))

sce.list.Standard <- variableFeatureSelection(seurat.lists = scelist, method = "Standard", nfeatures = 3000)
saveRDS(sce.list.Standard, file = "./doublet_high/sce.list.Standard.3000.rds")

sce.list.SCT <- variableFeatureSelection(seurat.lists = scelist, method = "SCT", nfeatures = 3000,return.only.var.genes = T)
saveRDS(sce.list.SCT, file = "./doublet_high/sce.list.SCT.3000.rds")

rm(scelist)

## Remove previous clustering results and doublet detection results
data.merge <- merge(sce.list.SCT[[1]], y = sce.list.SCT[2:length(sce.list.SCT)],add.cell.ids = names(sce.list.SCT))
index <- match(paste0("SCT_snn_res.", seq(0.1,1, by=0.1)), colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
index <- match(grep('DF.classification',colnames(data.merge@meta.data),value = T),colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
index <- match(grep('pANN',colnames(data.merge@meta.data),value = T),colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
View(data.merge@meta.data)

# assay=SCT
DefaultAssay(data.merge) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = sce.list.SCT, nfeatures = 3000)
VariableFeatures(data.merge) <- seurat.features.SCT
# assay=RNA
DefaultAssay(data.merge) <- "RNA"
seurat.features.RNA <- SelectIntegrationFeatures(object.list = sce.list.Standard, nfeatures = 3000)
VariableFeatures(data.merge) <- seurat.features.RNA
data.merge <- NormalizeData(data.merge, verbose = FALSE)
data.merge <- ScaleData(data.merge, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(data.merge))
DefaultAssay(data.merge) <- "SCT"
saveRDS(data.merge, file = "./doublet_high/data.merge.rds")

#####Evaluation of cellcycle and patient bias####
rm(list = setdiff(ls(),'data.merge'))
DefaultAssay(data.merge) <- "SCT"
pdf("1.QualityControl/doublet_high/filtered.statistics.pdf")
source('./function/colorPalettes.R')
VlnPlot(object = data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", cols = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))], pt.size = 0)
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
# Draw the distribution of the number of samples
cell.number <- as.data.frame(table(data.merge$orig.ident))
pdf("1.QualityControl/doublet_high/sample.distribution.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))],
          sort.by.groups=FALSE, #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

#### Assess the cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.merge <- CellCycleScoring(data.merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
data.merge <- RunPCA(data.merge, features = c(s.genes, g2m.genes))
pdf("1.QualityControl/doublet_high/cellCycle.afterMerge.pdf")
DimPlot(data.merge, dims = c(1, 2), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(1, 3), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(2, 3), reduction = "pca", group.by = "Phase")
dev.off()
saveRDS(data.merge, file = "./doublet_high/data.merge.rds")
####
rm(list = ls())
data.merge <- readRDS('./doublet_high/data.merge.rds')
set.resolutions <- seq(0.1, 1, by = 0.1)
#### First observe whether the clustering effect will depend on the sample
a <- data.merge
a <- RunPCA(a, npcs = 50, verbose = T)
pdf("1.QualityControl/doublet_high/merge.observe.batch.pdf")
a <- FindNeighbors(a, dims = 1:40, verbose = T)
a <- FindClusters(object = a, resolution = set.resolutions, verbose = T) 
clustree(a)
a <- RunUMAP(a, dims = 1:40)
DimPlot(object = a, reduction = 'umap',label = TRUE, split.by = "orig.ident",ncol = 3)
DimPlot(object = a, reduction = 'umap',label = TRUE, split.by = 'Phase')
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})
dev.off()

######correct batch effect####
options(future.globals.maxSize = 80000*1024^2)
plan("multisession",workers = 4)
rm(list = ls())
data.merge <- readRDS(('./doublet_high/data.merge.rds'))
dir.create('./2.Cluster')

####harmony
DefaultAssay(data.merge) <- 'RNA'
data.merge <- RunPCA(data.merge, verbose = FALSE, npcs = 100)
#If dims.use is not specified, all PCs will be used by default
#assay.use defaults to RNA. If the SCTransform standardization method is used, you need to specify assay.use="SCT" in RunHarmony
#Compare the differences between the two modes and find that the default is better and more complete
data.merge <- RunHarmony(object = data.merge, group.by.vars = "orig.ident", assay.use='RNA', verbose = FALSE)
PC <- 30
data.merge <- RunUMAP(data.merge, reduction = "harmony", dims = 1:PC, verbose = FALSE)
data.merge <- RunTSNE(data.merge, reduction = "harmony", dims = 1:PC, verbose = FALSE)
data.merge <- FindNeighbors(data.merge, dims = 1:PC, reduction = "harmony", verbose = FALSE) #使用harmony替代PCA
set.resolutions = seq(0.1, 1, by = 0.1)
data.merge <- FindClusters(data.merge, resolution = set.resolutions, verbose = FALSE)
pdf('2.Cluster/doublet high/RNA.Harmony.Integration.pc30.pdf')
p <- clustree(data.merge)
print(p)
p <- DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
p <- DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge, reduction = 'tsne',label = TRUE, group.by = paste0('RNA', "_snn_res.", x)) + NoLegend()
  print(p)
})
dev.off()
saveRDS(data.merge,file = './doublet_high/RNA.Harmony.Integration.PC30.rds')




####sct
sce.list.SCT <- readRDS('./doublet_highs/sce.list.SCT.3000.rds')
seurat.features <- SelectIntegrationFeatures(object.list = sce.list.SCT, nfeatures = 3000)
sce.list.SCT <- PrepSCTIntegration(object.list = sce.list.SCT, anchor.features = seurat.features)
seurat.anchors <- FindIntegrationAnchors(object.list = sce.list.SCT, normalization.method = "SCT", anchor.features = seurat.features)
saveRDS(seurat.anchors,file = './doublet_high/sct_anchors.rds')

seurat_anchors <- readRDS('./doublet_high/sct_anchors.rds')
data.merge.seurat <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")           
DefaultAssay(data.merge.seurat) <- "integrated"
data.merge.seurat <- ScaleData(data.merge.seurat)
data.merge.seurat <- RunPCA(data.merge.seurat, verbose = FALSE, npcs = 100)
pdf('2.Cluster/doublet high/SCT.SCT.Integration.pc40.pdf')
PC <- 40
set.resolutions = seq(0.1, 1, by = 0.1)
data.merge.seurat <- RunUMAP(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
data.merge.seurat <- RunTSNE(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
data.merge.seurat <- FindNeighbors(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
data.merge.seurat <- FindClusters(data.merge.seurat, resolution = set.resolutions, verbose = FALSE)
p <- clustree(data.merge.seurat)
print(p)
p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
  print(p)
})
p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
  print(p)
})

dev.off()
saveRDS(data.merge.seurat, file = "./doublet_high/SCT.SCT.Integration.PC40.rds")

####cca
sce.list.Standard <- readRDS('./doublet_high/sce.list.Standard.3000.rds')
seurat.features <- SelectIntegrationFeatures(object.list = sce.list.Standard, nfeatures = 3000)
seurat.anchors <- FindIntegrationAnchors(object.list = sce.list.Standard, anchor.features = seurat.features)
saveRDS(seurat.anchors,file = './doublet_high/cca_anchors.rds')

seurat_anchors <- readRDS('./doublet_high/cca_anchors.rds')
data.merge.seurat <- IntegrateData(anchorset = seurat_anchors)           
DefaultAssay(data.merge.seurat) <- "integrated"
data.merge.seurat <- ScaleData(data.merge.seurat)
data.merge.seurat <- RunPCA(data.merge.seurat, verbose = FALSE, npcs = 100)
pdf('2.Cluster/doublet high/RNA.standard.Integration.pc30.pdf')

PC <- 30
data.merge.seurat <- RunUMAP(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
data.merge.seurat <- RunTSNE(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
data.merge.seurat <- FindNeighbors(data.merge.seurat, reduction = "pca", dims = 1:PC, verbose = FALSE)
set.resolutions = seq(0.1, 1, by = 0.1)
data.merge.seurat <- FindClusters(data.merge.seurat, resolution = set.resolutions, verbose = FALSE)
p <- clustree(data.merge.seurat)
print(p)
p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge.seurat, reduction = 'umap',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
  print(p)
})
p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
print(p)
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = data.merge.seurat, reduction = 'tsne',label = TRUE, group.by = paste0("integrated_snn_res.", x)) + NoLegend()
  print(p)
})

dev.off()
saveRDS(data.merge.seurat, file = "./doublet_high/RNA.cca.Integration.PC30.rds")
