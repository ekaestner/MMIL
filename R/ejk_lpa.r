ejk_lpa <- function( dep_var_loc,
                     mdl_num,
                     out_put_loc ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/LPA_5/dep_var.mat'
  mdl_num     = 4
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/'  
  
  library( 'R.matlab' )
  library( 'tidyLPA' )
  library( 'dplyr' )
  library( 'pracma' )
  
  ## Load Data ###################################################################################################################################################
  # Dependent Variable
  dep_var = readMat(dep_var_loc) 
  dep_var_nme = names( dep_var$dep.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var_nme = dep_var_nme[-1]
  
  num_dep_var = ncol(dep_var)-1
  
  ## LPA ###################################################################################################################################################
  for( iF in (mdl_num-1):(mdl_num+1) ){
    
  dta_mdl_hld = dep_var[,2:ncol(dep_var)] %>%
    single_imputation() %>%
    estimate_profiles( iF, 
                       variances = c("equal" ),  # varying
                       covariances = c("zero" )) #%>% # varying
  #compare_solutions( statistics=c('AIC', 'LogLik', 'SABIC', 'Entropy' ))
  sbj_dta = get_data(dta_mdl_hld)
  
  sbj_dta_out = sbj_dta[ , 'Class' ]
  for ( iS in 1:nrow(sbj_dta_out) ){
        sbj_dta_out[ iS, 2 ] = sbj_dta[ iS,  paste( 'CPROB', sbj_dta_out[iS,1],sep='') ]
        }
  
  ## SCATTER PLOT #################################################
  use_col = pos_col[ 1:iF, ]
  if (iF>1){use_col = apply( use_col, 1, function(x) rgb(x[1], x[2], x[3]) )}else { use_col = rgb( use_col[1], use_col[2], use_col[3]) }
  jpeg( paste(out_put_loc, '/', 'model_', as.character(iF), '/', 'scatter_plot_colored_', as.character(iF), '.jpg', sep=''), width=1080, height=1080 )
  pairs( dep_var[,2:ncol(dep_var)] , pch=2, lower.panel=NULL, col=use_col[as.matrix(sbj_dta_out[,1])])
  dev.off()
  
  ## Save ###################################################################################################################################################
  write.csv( sbj_dta_out, paste( out_put_loc, '/', 'model_', as.character(iF), '/', 'model_', iF, '.csv', sep='') )
  }
  
}