####Using DoubletFinder function for Doublet detection

doubletDetect <- function(Seurat.object, PCs, doublet.rate = 0.076, annotation = "seurat_clusters", pN_value = 0.25, GT = FALSE, sct = FALSE){
  
  library(DoubletFinder) #Require cleanup of low-quality cells in advance
  T1.sample <- paramSweep_v3(Seurat.object, PCs = PCs, sct = sct)
  T1.gt.calls <- summarizeSweep(T1.sample, GT = GT)
  
  #Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep. 
  #Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions. 
  #If ground-truth doublet classifications are available, BCmvn is plotted along with mean ROC AUC for each pK.
  sweep.T1 <- find.pK(T1.gt.calls)
  pK_value <- as.numeric(as.character(sweep.T1$pK[sweep.T1$BCmetric == max(sweep.T1$BCmetric)])) #计算最优pK
  
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  plot(x = as.numeric(as.character(sweep.T1$pK)), y = sweep.T1$BCmetric, pch = 16,type="b", col = "blue",lty=1)
  abline(v=pK_value,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_value,max(sweep.T1$BCmetric),as.character(pK_value),pos = 4,col = "red")
  
  #potential doublet rate，
  nExp_poi <- round(doublet.rate*nrow(Seurat.object@meta.data))
  if(is.null(annotation)){
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi)
  }else{
    annotations <- Seurat.object@meta.data[, annotation]
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    pANN_value <- paste0("pANN_", pN_value, "_", pK_value, '_', nExp_poi)
    
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
    Seurat.object <- doubletFinder_v3(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = sct)
    label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi.adj)
  }
  Seurat.object@meta.data$Doublet <- Seurat.object@meta.data[, label]
  p <- DimPlot(Seurat.object,reduction = 'umap',group.by = 'Doublet') + 
    geom_text(aes(x = 0,y = 5),label = length(Seurat.object@meta.data$Doublet[Seurat.object@meta.data$Doublet=='Doublet'])) + 
    geom_text(aes(x = 0,y = 6),label = doublet.rate)
  print(p)
  return(Seurat.object)
}
