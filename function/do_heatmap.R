DoHeatmap_my=function(sce, features, slot="scale.data", assay="RNA", main="",cluster_cols,group.by="seurat_clusters",color, rotate=F,angle_col,cluster_rows){
  
  # 1. get df from Seurat: scaled gene expression
  .data = GetAssayData(	object = sce, assay = assay, slot = slot)[features,]
  .group=as.character( as.data.frame(sce@meta.data)[, group.by] )
  #dim(.data)
  #2. get mean by sample
  .data2=sapply( split( as.data.frame(t(.data)), .group), function(x){
    colMeans(x)
  })
  .data2 = .data2[,levels(sce)]
  # remove all 0 rows
  .data2=.data2[apply(.data2, 1, sd)>0,]
  
  #3. clean large outliers
  #.data2[.data2<cut_dims[1]] <-cut_dims[1]
  #.data2[.data2>cut_dims[2]] <-cut_dims[2]
  if(0){ # how is this fig?
    tmp=AverageExpression(sce, slot="data", group.by = group.by)$RNA[features, ]
    dim(tmp)
    pheatmap(log(tmp+1), border_color = NA,
             clustering_method = "ward.D2",
             scale = "row",
             main="try4"
             )
  }
  
  params2=list(mat=.data2,
               main=main,
               scale = 'row',
               clustering_method = "ward.D2",
               color = color,
               angle_col = angle_col,
               cluster_cols = cluster_cols,
               cluster_rows = cluster_rows)
    if(rotate){
    params2$mat=t(.data2)
    params2$scale="column"
  }
  library(pheatmap)
  do.call(pheatmap::pheatmap, params2)
}