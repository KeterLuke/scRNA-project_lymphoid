2####load packages####
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
library(monocle)
set.seed(101)
library(future)
plan("multisession", workers = 1) 
options(future.globals.maxSize = 80000 * 1024^2)
setwd('../')
getwd()
source('./function/colorPalettes.R')
library(plot1cell)
data.merge <- readRDS('./2.Cluster/data.merge.pro.rds')
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
CD8 <- readRDS('./3.1Subcluster_T/Annotate/CD8/CD8.pro.rds')
B <- readRDS('./3.2Subcluster_B/sub.B.pro.rds')
cd4_celltype <- FetchData(CD4,vars = 'minor_celltype')
cd4_celltype$minor_celltype <- as.character(cd4_celltype$minor_celltype)
cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4')] <- 'CD4T'
# cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4Tn')] <- 'CD4Tn'
# cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4Tm')] <- 'CD4Tm'
# cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4Treg')] <- 'CD4Treg'
# cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4Tst')] <- 'CD4Tstd'
# cd4_celltype$minor_celltype[grepl(cd4_celltype$minor_celltype,pattern = 'CD4Tex')] <- 'CD4Tex'
table(cd4_celltype$minor_celltype)

cd8_celltype <- FetchData(CD8,vars = 'minor_celltype')
cd8_celltype$minor_celltype <- as.character(cd8_celltype$minor_celltype)
cd8_celltype$minor_celltype[grepl(cd8_celltype$minor_celltype,pattern = 'CD8Tn')] <- 'CD8Tn'
cd8_celltype$minor_celltype[grepl(cd8_celltype$minor_celltype,pattern = 'CD8Te')] <- 'CD8Teff'
table(cd8_celltype$minor_celltype)

b_celltype <- FetchData(B,vars = 'minor_celltype')
b_celltype$minor_celltype <- as.character(b_celltype$minor_celltype)
b_celltype$minor_celltype[grepl(b_celltype$minor_celltype,pattern = 'Bpro')] <- 'Bpro'
b_celltype$minor_celltype[grepl(b_celltype$minor_celltype,pattern = 'Bgc')] <- 'Bgc'
b_celltype$minor_celltype[!b_celltype$minor_celltype%in%c('Bpro','Bgc')] <- 'B'
table(b_celltype$minor_celltype)

m_celltype <- FetchData(M,vars = 'minor_celltype')
m_celltype$minor_celltype <- as.character(m_celltype$minor_celltype)
m_celltype$minor_celltype[grepl(m_celltype$minor_celltype,pattern = 'Macro')] <- 'Macro'
m_celltype$minor_celltype[grepl(m_celltype$minor_celltype,pattern = 'cDC1')] <- 'cDC1'
m_celltype$minor_celltype[grepl(m_celltype$minor_celltype,pattern = 'cDC2')] <- 'cDC2'
m_celltype$minor_celltype[grepl(m_celltype$minor_celltype,pattern = 'Doublet')] <- NA

table(m_celltype$minor_celltype)

data.merge$minor_celltype<- data.merge$major_celltype
index <- match(rownames(cd4_celltype),rownames(data.merge@meta.data))
data.merge$minor_celltype <- as.character(data.merge$minor_celltype)
data.merge@meta.data$minor_celltype[index] <- cd4_celltype$minor_celltype
index <- match(rownames(cd8_celltype),rownames(data.merge@meta.data))
data.merge@meta.data$minor_celltype[index] <- cd8_celltype$minor_celltype
index <- match(rownames(b_celltype),rownames(data.merge@meta.data))
data.merge@meta.data$minor_celltype[index] <- b_celltype$minor_celltype
index <- match(rownames(m_celltype),rownames(data.merge@meta.data))
data.merge@meta.data$minor_celltype[index] <- m_celltype$minor_celltype

table(data.merge$minor_celltype,useNA = 'always')
data.merge <- data.merge[,!is.na(data.merge$minor_celltype)]
length(unique(data.merge$minor_celltype))
Idents(data.merge) <- data.merge$minor_celltype
dir.create('./5.merge.minor.major')
length(unique(data.merge$major_celltype))
pdf('./5.merge.minor.major/merge.annotation.pdf',width = 10,height = 8)
DimPlot(data.merge,reduction = 'umap',group.by = 'major_celltype',label = T,cols = Palettes[['mycols_8']])
DimPlot(data.merge,reduction = 'tsne',group.by = 'major_celltype',label = T,cols = Palettes[['mycols_8']])
DimPlot(data.merge,reduction = 'umap',group.by = 'minor_celltype',label = T,repel = T,cols = Palettes[['circus']])
DimPlot(data.merge,reduction = 'tsne',group.by = 'minor_celltype',label = T,repel = T,cols = Palettes[['circus']])
dev.off()
saveRDS(data.merge,file = './5.merge.minor.major/data.merge.pro.rds')

####use plot1cell to visualize####
library(plot1cell)
data.merge <- readRDS('./5.merge.minor.major/data.merge.pro.rds')
table(Idents(data.merge))
length(unique(Idents(data.merge)))
circ_data <- prepare_circlize_data(data.merge, scale = 0.8 )
cluster_colors<- Palettes[['circus']][1:length(unique(Idents(data.merge)))]
group_colors<-rand_color(length(names(table(data.merge$group))))
rep_colors<-rand_color(length(names(table(data.merge$orig.ident))))
pdf('5.merge.minor.major/circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

data.merge <- readRDS('5.merge.minor.major/data.merge.pro.rds')
CD4 <- readRDS('./3.1Subcluster_T/Annotate/CD4/CD4.pro.rds')
data.merge$minor_celltype <- as.character(data.merge$major_celltype)
data.merge$minor_celltype[data.merge$major_celltype=='T'] <- 'Other T'
index <- match(colnames(CD4[,CD4$minor_celltype=='C9 CD4Tex_pre-MIR155HG']),colnames(data.merge))
data.merge$minor_celltype[index] <- "Tex_pre-MIR155HG"
data.merge$minor_celltype[data.merge$major_celltype%in%c('Myeloid','LAMP3+cDC','pDC')] <- 'Myeloid'
table(data.merge$minor_celltype)
pdf('./5.merge.minor.major/cellmarker.observe.in.allcells.pdf',width = 10)
DotPlot(data.merge,group.by = 'minor_celltype',features = c('MIR155HG','LAG3','C1orf228','TSHZ2','CCR7','SESN3'),scale.by = 'size',col.min = -.5) + 
  scale_color_gradientn(colours = Palettes[['blueRed']])
genelist <- c('MIR155HG','LAG3','C1orf228','TSHZ2','CCR7','SESN3')
exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(genelist, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]

# Calculate the average expression of the cell type of each cluster
cluster.exp <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$minor_celltype, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.exp, "range", 2) ##perform 0-1 standardization for all clusters per gene
Heatmap(t(cluster.score.normailzed), 
        column_names_rot = 45,cluster_columns = F, cluster_rows = F, 
        show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 10),rect_gp = gpar(col = 'white'),
        column_names_gp = gpar(fontsize =12),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")),col = Palettes[['blueRed']])
dev.off()

