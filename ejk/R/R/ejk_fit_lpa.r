ejk_fit_lpa <- function( dep_var_loc,
                             out_put_loc ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/LPA/dep_var.mat'
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/'  
  
  library( 'R.matlab' )
  library( 'tidyLPA' )
  library( 'dplyr' )
  library( 'pracma' )
  library( 'psych' )
  
  ## Load Data ###################################################################################################################################################
  # Dependent Variable
  dep_var = readMat(dep_var_loc) 
  dep_var_nme = names( dep_var$dep.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var_nme = dep_var_nme[-1]
  
  ## SCATTER PLOT #################################################
  jpeg( paste(out_put_loc, '/', 'scatter_plot.jpg', sep=''), width=1080, height=1080 )
  pairs( dep_var[,2:ncol(dep_var)] , pch=2, lower.panel=NULL)
  dev.off()
  
  jpeg( paste(out_put_loc, '/', 'scatter_plot_fancy.jpg', sep=''), width=1080, height=1080 )
  pairs.panels( dep_var[,2:ncol(dep_var)] , method='pearson', density=FALSE, ellipses=FALSE)
  dev.off()
  
  ## LPA ###################################################################################################################################################
  
  tic()
  dta_mdl_hld = dep_var[,2:ncol(dep_var)] %>%
                  single_imputation() %>%
                  estimate_profiles(1:15, 
                      variances = c("equal" ),  # varying
                      covariances = c("zero" )) #%>% # varying
                  #compare_solutions( statistics=c('AIC', 'LogLik', 'SABIC', 'Entropy' ))
  toc()
  
  ## LPA ###################################################################################################################################################
  dta_mdl_fit = get_fit(dta_mdl_hld)
  
  ahp_hld = AHP(dta_mdl_fit)
  
  #library(fpc)
  #cluster.stats( dist(dep_var[,2:ncol(dep_var)]),  get_data(dta_mdl_hld$model_1_class_1)$Class, get_data(dta_mdl_hld$model_1_class_2)$Class )
  
  lrt_out = matrix(,14,4)
  for( iM in 2:15 ){
  lrt_hld = calc_lrt( nrow(dep_var), dta_mdl_fit$LogLik[iM-1], iM-1, 1, dta_mdl_fit$LogLik[iM], 7, iM)
  lrt_out[ iM-1, 1] = lrt_hld[1]
  lrt_out[ iM-1, 2] = lrt_hld[2]
  lrt_out[ iM-1, 3] = lrt_hld[3]
  lrt_out[ iM-1, 4] = lrt_hld[4]
  }
  
  dta_mdl_out = dta_mdl_fit[ , c('AIC', 'BIC', 'SABIC', 'AWE' ) ]
  
  dta_fit_out = compare_solutions( dta_mdl_hld, statistics=c('AIC', 'LogLik', 'SABIC', 'Entropy' ) )
  
  ## Save ###################################################################################################################################################
  write.csv( lrt_out, paste( out_put_loc, '/', 'model_fit_lrt.csv',sep='') )
  write.csv( ahp_hld, paste( out_put_loc, '/', 'model_fit_ahp.csv',sep='') )
  write.csv( dta_mdl_out, paste( out_put_loc, '/', 'model_fit_measures.csv',sep='') )
  write.csv( dta_fit_out$best, paste( out_put_loc, '/', 'model_fit_compare.csv',sep='') )
  
  

  
}