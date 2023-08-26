# @description: quality control,remove doublets,remove batch effect and integrate samples
# this work_lympho 6 samples
#### load package (always load all the package before you run the code)####
library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
library(colourpicker)
library(cowplot)
plan("multisession", workers = 11) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
source('./function/ratio_plot.R')
####read samples####
n4_6 <- Read10X('./raw/N4_6/filtered_feature_bc_matrix/')
T4_6 <- Read10X('./raw/T4_6/filtered_feature_bc_matrix/')
N6_15 <- Read10X('./raw/N6_15/filtered_feature_bc_matrix/')
T6_15 <- Read10X('./raw/T6_15/filtered_feature_bc_matrix/')
N6_8 <- Read10X('./raw/N6_8/filtered_feature_bc_matrix/')
T6_8 <- Read10X('./raw/T6_8/filtered_feature_bc_matrix/')
####Preliminary filtration####
#min.cell >3 & min.features >200   
N4_6 <- CreateSeuratObject(counts = N4_6,
                          min.cells = 3,
                          min.features = 200,
                          project = 'N4_6')
T4_6 <- CreateSeuratObject(counts = T4_6,
                           min.cells = 3,
                           min.features = 200,
                           project = 'T4_6')
N6_15 <- CreateSeuratObject(counts = N6_15,
                           min.cells = 3,
                           min.features = 200,
                           project = 'N6_15')
T6_15 <- CreateSeuratObject(counts = T6_15,
                            min.cells = 3,
                            min.features = 200,
                            project = 'T6_15')
N6_8 <- CreateSeuratObject(counts = N6_8,
                            min.cells = 3,
                            min.features = 200,
                            project = 'N6_8')
T6_8 <- CreateSeuratObject(counts = T6_8,
                           min.cells = 3,
                           min.features = 200,
                           project = 'T6_8')
scelist <- list(N4_6 = N4_6,
                T4_6 = T4_6,
                N6_15 = N6_15,
                T6_15 = T6_15,
                N6_8 = N6_8,
                T6_8 = T6_8)
sce = merge(scelist[[1]], y = scelist[2:length(scelist)],add.cell.ids <- names(scelist))
sce$sample <- sce$orig.ident
table(sce$sample)
View(sce@meta.data)
####quality control####
sce$nUMI <- sce$nCount_RNA
sce$nGene <- sce$nFeature_RNA
sce$log10GenesPerUMI <- log10(sce$nGene)/log10(sce$nUMI)
## Mitochondrial gene ratio
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
## Ribosome gene ratio
sce[["percent.rb"]] <- PercentageFeatureSet(sce, pattern = "^RPL|^RPS")
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
## HBC gene ratio
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HBA|^HBB")
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
dir.create('./1.QualityControl')
pdf(file = "1.QualityControl/count.feature.mt.pdf",width = 10)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"),ncol = 3)
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + ploT1

ggdensity(sce@meta.data, x = "nCount_RNA", title = "nCount")
ggdensity(sce@meta.data, x = "nFeature_RNA", title = "nFeature") + geom_vline(xintercept = 500)
ggdensity(sce@meta.data, x = "log10GenesPerUMI", title = "log10GenesPerUMI") + geom_vline(xintercept = 0.8)
ggdensity(sce@meta.data, x = "percent.rb", title = "percent.rb") +  geom_vline(xintercept = 8)
dev.off()

#### Detect the resolution parameters by each sample cluster ####
##After the parameters are determined, you can block them without executing [test]
sce_obj <- subset(sce, subset = nFeature_RNA > 600 & percent.mt < 10 & nCount_RNA > 1000 &
                  nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.rb > 8 & percent.hb < 0.1 & 
                  log10GenesPerUMI > 0.8) 

pdf("1.QualityControl/filtered.statistics.pdf")
source('./function/colorPalettes.R')
VlnPlot(object = sce_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rb','percent.hb'), ncol = 2, group.by = "orig.ident", cols = Palettes$group_pal[1:length(unique(sce_obj@meta.data$orig.ident))], pt.size = 0)
FeatureScatter(object = sce_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = sce_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
sce_obj$mito.ratio <- sce_obj$percent.mt/100
source('./function/do_qc_plot.R')
do_qc_plot(sce_obj,dir = './1.QualityControl')
sce_obj <- SplitObject(sce_obj,split.by = 'orig.ident')
##N4_6
N4_6 <- sce_obj[['N4_6']]
saveRDS(N4_6,file = './N4_6_afterqc.RDS')
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

####N4_6 7923 samples
sce_obj <- readRDS('./N4_6_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:12, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
pdf(file = "1.QualityControl/N4_6_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:12)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './N4_6_afterqc.RDS')

####T4_6 6678 samples
sce_obj <- readRDS('./T4_6_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:9, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
pdf(file = "1.QualityControl/T4_6_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:9)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T4_6_afterqc.RDS')

####N6_8 15097 samples
sce_obj <- readRDS('./N6_8_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:15, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE) 
pdf(file = "1.QualityControl/N6_8_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:15)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './N6_8_afterqc.RDS')

####T6_8 11091 samples
sce_obj <- readRDS('./T6_8_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:17, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE)
pdf(file = "1.QualityControl/T6_8_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:17)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T6_8_afterqc.RDS')


####N6_15 8362 samples
sce_obj <- readRDS('./N6_15_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:12, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE)
pdf(file = "1.QualityControl/N6_15_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:12)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './N6_15_afterqc.RDS')

####T6_15 6547 samples
sce_obj <- readRDS('./T6_15_afterqc.RDS')
sce_obj <- NormalizeData(sce_obj,verbose = F)
sce_obj <- FindVariableFeatures(sce_obj,nfeatures = 3000,verbose = F)
sce_obj <- ScaleData(sce_obj,vars.to.regress = c("nCount_RNA", "percent.mt"),features = VariableFeatures(sce_obj),verbose = F)
sce_obj <- RunPCA(sce_obj, npcs = 100, verbose = FALSE)
do_Elbow_quantitative(sce_obj,harmony = F)
sce_obj  <- FindNeighbors(object = sce_obj , dims = 1:8, verbose = FALSE)
sce_obj  <- FindClusters(object = sce_obj , resolution = set.resolutions, verbose = FALSE)
pdf(file = "1.QualityControl/T6_15_PCA-test.pdf",width =15)
clustree(sce_obj)
sce_obj  <- RunUMAP(sce_obj , dims = 1:8)
sce.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = sce_obj, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})

dev.off()
saveRDS(sce_obj,file = './T6_15_afterqc.RDS')
rm(list = ls())


#### remove doublet####
library(DoubletFinder) # Require cleanup of low-quality cells in advance
source('./function/doubletDetect.R')
dir.create('./doublet')
####N4_6
sce <- readRDS('./N4_6_afterqc.RDS')
pdf("1.QualityControl/N4_6_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:12, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.044
saveRDS(sce,file = './doublet/N4_6_afterqc.RDS')

####T4_6
sce <- readRDS('./T4_6_afterqc.RDS')
pdf("1.QualityControl/T4_6_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:9, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.056
saveRDS(sce,file = './doublet/T4_6_afterqc.RDS')

####N6_8
sce <- readRDS('./N6_8_afterqc.RDS')
pdf("1.QualityControl/N6_8_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:15, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.136
saveRDS(sce,file = './doublet/N6_8_afterqc.RDS')

####T6_8
sce <- readRDS('./T6_8_afterqc.RDS')
pdf("1.QualityControl/T6_8_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:17, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.095
saveRDS(sce,file = './doublet/T6_8_afterqc.RDS')


####N6_15
sce <- readRDS('./N6_15_afterqc.RDS')
pdf("1.QualityControl/N6_15_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:12, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.070
saveRDS(sce,file = './doublet/N6_15_afterqc.RDS')

####T6_15
sce <- readRDS('./T6_15_afterqc.RDS')
pdf("1.QualityControl/T6_15_doublet.pdf")
sce <- doubletDetect(Seurat.object = sce, PCs = 1:8, doublet.rate = ncol(sce)*8*1e-6, annotation = "RNA_snn_res.0.1", sct = F)
dev.off()
##select:rate = 0.054
saveRDS(sce,file = './doublet/T6_15_afterqc.RDS')

#####subset singlet####
rm(list = ls())
#N4_6
sce <- readRDS('./doublet/N4_6_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "./doublet/N4_6.single.RDS")

#T4_6
sce <- readRDS('./doublet/T4_6_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet/T4_6.single.RDS")

#N6_8
sce <- readRDS('./doublet/N6_8_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet/N6_8.single.RDS")

#T6_8
sce <- readRDS('./doublet/T6_8_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet/T6_8.single.RDS")

#N6_15
sce <- readRDS('./doublet/N6_15_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet/N6_15.single.RDS")

#T6_15
sce <- readRDS('./doublet/T6_15_afterqc.RDS')
sce.sing <- subset(sce, subset = Doublet == "Singlet")
saveRDS(sce.sing, file = "doublet/T6_15.single.RDS")

#### merge data and correct the batch effect####
rm(list = ls())
N4_6 <- readRDS('./doublet/N4_6.single.RDS')
DefaultAssay(N4_6) <- "RNA"
T4_6 <- readRDS('./doublet/T4_6.single.RDS')
DefaultAssay(T4_6) <- "RNA"
N6_8 <- readRDS('./doublet/N6_8.single.RDS')
DefaultAssay(N6_8) <- "RNA"
T6_8 <- readRDS('./doublet/T6_8.single.RDS')
DefaultAssay(T6_8) <- "RNA"
N6_15 <- readRDS('./doublet/N6_15.single.RDS')
DefaultAssay(N6_15) <- "RNA"
T6_15 <- readRDS('./doublet/T6_15.single.RDS')
DefaultAssay(T6_15) <- "RNA"


source('./function/variableFeatureSelection.R')
scelist <- list(N4_6 = N4_6,
                T4_6 = T4_6,
                N6_8 = N6_8,
                T6_8 = T6_8,
                N6_15 = N6_15,
                T6_15 = T6_15)
rm(list = setdiff(ls(),c('scelist','variableFeatureSelection')))

sce.list.Standard <- variableFeatureSelection(seurat.lists = scelist, method = "Standard", nfeatures = 3000)
saveRDS(sce.list.Standard, file = "./sce.list.Standard.3000.rds")

## Remove previous clustering results and doublet detection results
data.merge <- merge(sce.list.Standard[[1]], y = sce.list.Standard[2:length(sce.list.Standard)])
index <- match(paste0("RNA_snn_res.", seq(0.1,1, by=0.1)), colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
index <- match(grep('DF.classification',colnames(data.merge@meta.data),value = T),colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
index <- match(grep('pANN',colnames(data.merge@meta.data),value = T),colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
View(data.merge@meta.data)

# assay=RNA
DefaultAssay(data.merge) <- "RNA"
seurat.features.RNA <- SelectIntegrationFeatures(object.list = sce.list.Standard, nfeatures = 3000)
VariableFeatures(data.merge) <- seurat.features.RNA
####remove cells coexpress cell type markers####
cells1 <- WhichCells(data.merge,expression = CD3D&CD3E&CD79A&MS4A1 > 0)
cells <- colnames(data.merge)
cells <- setdiff(cells,cells1)
data.merge <- subset(data.merge,cells = cells)
data.merge <- NormalizeData(data.merge, verbose = FALSE)
data.merge <- ScaleData(data.merge, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = VariableFeatures(data.merge))
saveRDS(data.merge, file = "./data.merge.rds")

#####Evaluation of cellcycle and patient bias####
rm(list = setdiff(ls(),'data.merge'))
DefaultAssay(data.merge) <- "RNA"
pdf("1.QualityControl/doublet.removed.statistics.pdf")
source('./function/colorPalettes.R')
VlnPlot(object = data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", cols = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))], pt.size = 0)
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
source('./function/do_qc_plot.R')
do_qc_plot(data.merge,dir = '1.QualityControl')
# Draw the distribution of the number of samples
cell.number <- as.data.frame(table(data.merge$orig.ident))
pdf("1.QualityControl/sample.distribution.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))],
          sort.by.groups=FALSE, #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()
#### Assess the cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.merge <- CellCycleScoring(data.merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
data.merge <- RunPCA(data.merge, features = c(s.genes, g2m.genes),verbose = F)
pdf("1.QualityControl/cellCycle.afterMerge.pdf")
DimPlot(data.merge, dims = c(1, 2), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(1, 3), reduction = "pca", group.by = "Phase")
DimPlot(data.merge, dims = c(2, 3), reduction = "pca", group.by = "Phase")
dev.off()
saveRDS(data.merge, file = "./data.merge.rds")
####
rm(list = ls())
data.merge <- readRDS('./data.merge.rds')
set.resolutions <- seq(0.1, 1, by = 0.1)
######correct batch effect####
options(future.globals.maxSize = 80000*1024^2)
plan("multisession",workers = 1)
rm(list = ls())
data.merge <- readRDS('./data.merge.rds')
dir.create('./2.Cluster')

####harmony
DefaultAssay(data.merge) <- 'RNA'
data.merge <- RunPCA(data.merge, verbose = FALSE, npcs = 100)
#If dims.use is not specified, all PCs will be used by default
#assay.use defaults to RNA. If the SCTransform standardization method is used, you need to specify assay.use="SCT" in RunHarmony
#Compare the differences between the two modes and find that the default is better and more complete
data.merge <- RunHarmony(object = data.merge, group.by.vars = "orig.ident", assay.use='RNA', verbose = FALSE)
source('./function/do_Elbow_quantitative_modified_by_lhc.R')
do_Elbow_quantitative(data.merge,harmony = T)
PC <- 35
data.merge <- RunUMAP(data.merge, reduction = "harmony", dims = 1:PC, verbose = FALSE)
data.merge <- RunTSNE(data.merge, reduction = "harmony", dims = 1:PC, verbose = FALSE)
data.merge <- FindNeighbors(data.merge, dims = 1:PC, reduction = "harmony", verbose = FALSE) #使用harmony替代PCA
set.resolutions = seq(0.1, 1, by = 0.1)
data.merge <- FindClusters(data.merge, resolution = set.resolutions, verbose = FALSE)
pdf('2.Cluster/RNA.Harmony.Integration.pc35.pdf')
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
saveRDS(data.merge,file = './RNA.Harmony.Integration.PC35.rds')

####select: RNA harmony pc35,res:0.2
data.merge <- readRDS('./RNA.Harmony.Integration.PC35.rds')
data.merge <- RunUMAP(data.merge,reduction = "harmony", dims = 1:35,seed.use = 50, verbose = FALSE)
data.merge <- RunTSNE(data.merge,reduction = "harmony", dims = 1:35,seed.use = 50, verbose = FALSE)
pdf('2.Cluster/RNA.Harmony.Integration.pc35.pro.pdf')
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
saveRDS(data.merge,file = './data.merge.pro.rds')
####The expression of the classic marker#####
data.merge <- readRDS('./data.merge.pro.rds')
pdf('./2.Cluster/observe.expression.T.B.pdf')
FeaturePlot(data.merge,features = c('MS4A1','CD19','CD79A','CD79B'),reduction ='umap',order = T,min.cutoff = 1)
FeaturePlot(data.merge,features = c('CD7','CD3D','CD3E','CD3G'),reduction ='umap',order = T,min.cutoff = 1)
dev.off()
data.merge$seurat_clusters <- data.merge$RNA_snn_res.0.2
data.merge$seurat_clusters <- factor(data.merge$seurat_clusters,levels = as.character(0:15))
table(data.merge$seurat_clusters)
cell.type.markers <- read.table(file = "2.Cluster/CellMarker_lowres.txt", header = T, stringsAsFactors = F, sep = "\t")

exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$seurat_clusters, mean)
  return(a)
})
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
pdf("2.Cluster/cluster.signature.expression.pdf")
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
a <- data.merge
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
DotPlot_ByColumnList(object = a,features = gene_list,group.by = "RNA_snn_res.0.4",dot.scale = 4)
dev.off()
####qc plot####
DefaultAssay(data.merge) <- "RNA"
pdf("1.QualityControl/doublet.removed.statistics.pdf",width = 9,height = 6)
source('./function/colorPalettes.R')
VlnPlot(object = data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample", cols = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))], pt.size = 0)
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
source('./function/do_qc_plot.R')
do_qc_plot(data.merge,dir = '1.QualityControl')
# Draw the distribution of the number of samples
cell.number <- as.data.frame(table(data.merge$orig.ident))
pdf("1.QualityControl/sample.distribution.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))],
          sort.by.groups=FALSE, #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

saveRDS(data.merge,file = './2.Cluster/filtered/data.merge.filtered.pc20.seed45.rds')
