ejk_lpa <- function( dep_var_loc,
                     pos_col,
                     out_put_loc ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive_kmeans/Kmeans/dep_var.mat'
  pos_col_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive_kmeans/Kmeans/pos_col.mat'
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive_kmeans/Kmeans/'  
  
  library( 'R.matlab' )
  library( 'dplyr' )
  library( 'clValid' )
  library( 'NbClust' )
  
  ## Load Data ###################################################################################################################################################
  # Dependent Variable
  dep_var = readMat(dep_var_loc) 
  dep_var_nme = names( dep_var$dep.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var_nme = dep_var_nme[-1]
  
  num_dep_var = ncol(dep_var)-1
  
  # Possible Colors
  pos_col = readMat(pos_col_loc) 
  pos_col = as.matrix(pos_col$pos.col)

  ## K-Means ###################################################################################################################################################
  cls_fit = NbClust( dep_var[, -1], method='kmeans' )
  out_fit = t(cls_clc$Best.nc)
  write.csv( out_fit, paste( out_put_loc, '/', 'model_fit.csv', sep='') )
  
  mdl_num = max(cls_fit$Best.partition)
  
  out_cls = matrix( , nrow(dep_var), mdl_num+1)
  for( iF in (mdl_num-1):(mdl_num+1) ){
  fit <- kmeans( dep_var[, -1], iF)
  out_cls[,iF] = t(fit$cluster)
  }
  write.csv( out_cls, paste( out_put_loc, '/', 'model_cluster.csv', sep='') )
  
  dun_out = matrix( , 15, 1)
  for( iM in 2:15 ){
    fit <- kmeans( dep_var[, -1], iM)
    dun_out[iM,1] = dunn( clusters=fit$cluster, Data=dep_var[, -1])
  }
  write.csv( dun_out, paste( out_put_loc, '/', 'model_dunn.csv', sep='') )

  ## Colored Plot ###################################################################################################################################################
  for( iF in (mdl_num-1):(mdl_num+1) ){
    use_col = pos_col[ 1:iF, ]
    if (iF>1){use_col = apply( use_col, 1, function(x) rgb(x[1], x[2], x[3]) )}else { use_col = rgb( use_col[1], use_col[2], use_col[3]) }
    jpeg( paste(out_put_loc, '/', 'scatter_plot_colored_', as.character(iF), '.jpg', sep=''), width=1080, height=1080 )
    pairs( dep_var[,2:ncol(dep_var)] , pch=2, lower.panel=NULL, col=use_col[out_cls[,iF]])
    dev.off()
  }
}
