do_Elbow_quantitative<-function(seurat_integrated,harmony=T){
  #seurat_integrated<-tiss_subset_real_tumor 
  # Determine percent of variation associated with each PC
  if(harmony){pct <- seurat_integrated[["harmony"]]@stdev / sum(seurat_integrated[["harmony"]]@stdev) * 100
  }
  else{pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
  }

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  
  # Elbow plot to visualize 
  library(ggplot2)
  elbow=ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  
  
  print(elbow)
  res=list(fig=elbow,pcs=pcs)
  return(res)
  
}
