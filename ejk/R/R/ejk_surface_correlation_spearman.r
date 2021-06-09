# 2 (N=18, N=20) x 3 (No covariate, Baseline Score) x 4 (LM2, LDFR, BVMT, VPAII) = 24 brains

ejk_surface_correlation <- function( lhs_srf_var_loc,
                                     rhs_srf_var_loc,
                                     cor_loc,
                                     grp_loc, 
                                     cov_loc,
                                     iG,
                                     fst_nme,
                                     out_put_loc ) {
  
  lhs_srf_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data/surf_wmparc_fa_lhs_sm313.mat'
  rhs_srf_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data/surf_wmparc_fa_rhs_sm313.mat'
  cor_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst/cor_var.mat'
  grp_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst/grp_var.mat'
  cov_loc = ''
  iG = 1
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst_spearman/'
  fst_nme = 'tle_post_3T_ATLonly_left' # 3_left # left # 'left_dominant'
  
  library( R.matlab )
  library( rstatix )
  library( ggpubr )
  library( dplyr )
  library( emmeans )
  library( parallel )
  library( doParallel )
  library( foreach )
  library( stringr )
  library( broom )
  library( tidyr )
  library(jtools)
  
  ## Load Data ###################################################################################################################################################
  # Correlation Variable
  cor_var = readMat(cor_loc) 
  cor_var_nme = names( cor_var$cor.var[,,1] )
  cor_var = as.data.frame( lapply(cor_var$cor.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=cor_var_nme )
  cor_var_nme = cor_var_nme[-1]
  
  # Group Variable
  grp_var = readMat(grp_loc) 
  grp_var_nme = names( grp_var$grp.var[,,1] )
  grp_var = as.data.frame( lapply(grp_var$grp.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=grp_var_nme )
  for (iGR in 2:ncol(grp_var))
  {
    grp_var[,iGR] = factor( grp_var[,iGR], exclude='N/A')
  }
  grp_var_nme = grp_var_nme[-1]
  
  ## ANCOVA #################################################################################################################################################
  srf_dta_nme = c('lhs_dep_var', 'rhs_dep_var')
  
  for( iH in 1:2){
    
    # Dependent Variable
    if (iH==1 ){ 
      srf_dta = readMat(lhs_srf_var_loc)
      dep_var = t(as.data.frame( lapply(srf_dta$srf.dta.sbj, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE ))
      srf_dta = srf_dta$srf.dta 
    } else if (iH==2){ 
      srf_dta = readMat(rhs_srf_var_loc) 
      dep_var = t(as.data.frame( lapply(srf_dta$srf.dta.sbj, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE ))
      srf_dta = srf_dta$srf.dta }
    
    colnames(dep_var) = c('sbj.nme')
    
    dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
    
    # Set up backend for parallelizing
    ncores <- detectCores() - 2
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    # Data merge
    use_dta = merge( dep_var, cor_var, by="sbj.nme" );
    use_dta = merge( use_dta, grp_var, by="sbj.nme" );

    # Run test for each vertex individually
    output <- foreach( iV = 1:163842, .combine = 'combine', .packages = c("tidyr", "broom", "emmeans", "pracma", "stringr","jtools"), .multicombine=TRUE, #
                       .init=list(list(), list()) ) %dopar% {
                         
                         # Put together data
                         use_dta_col = as.data.frame(srf_dta[, iV], fix.empty.names=FALSE )
                         colnames(use_dta_col) = c(srf_dta_nme[iH])
                         use_dta_col = cbind( dep_var, use_dta_col )
                         use_dta_par = merge( use_dta, use_dta_col, by="sbj.nme" )
                         use_dta_par = subset( use_dta_par, match(use_dta_par[[grp_var_nme[1]]],fst_nme)==1 )
                         
                         # Running cor.test()
                         analysis = cor.test( use_dta_par[[srf_dta_nme[iH]]], use_dta_par[[cor_var_nme[iG]]], method='spearman')
                         
                         # Putting result into one list
                         rvl <- analysis$estimate
                         pvl <- analysis$p.value
                         
                         # Current progress
                         list( rvl, pvl )
                         
                       }
    
    # Free Up CPU and Memory
    stopCluster(cl)
    if (iH==1 ){ 
      iV = 5483
    } else if (iH==2){ 
      iV = 39922
    }
    use_dta_col = as.data.frame(srf_dta[, iV], fix.empty.names=FALSE )
    colnames(use_dta_col) = c(srf_dta_nme[iH])
    use_dta_col = cbind( dep_var, use_dta_col )
    use_dta_par = merge( use_dta, use_dta_col, by="sbj.nme" )
    use_dta_par = subset( use_dta_par, match(use_dta_par[[grp_var_nme[1]]],fst_nme)==1 )
    write.csv(use_dta_par, paste(out_put_loc,"/","exact_data_",iH,".csv",sep=''))
    
    rm('srf_dta','use_dta', 'dep_var')
    
    # Create the file names of the variables 
    rvl_fle_nme = paste( out_put_loc,"/", "rvalues_" , srf_dta_nme[iH], "_", fst_nme, ".mat", sep = "")
    pvl_fle_nme = paste( out_put_loc,"/", "pvalues_" , srf_dta_nme[iH], "_", fst_nme, ".mat", sep = "")
    
    # Write result into a matrix for matlab
    output = unlist(output)
    
    writeMat(rvl_fle_nme, rvalues = as.numeric( output[seq(1,length(output),2)] ))
    writeMat(pvl_fle_nme, pvalues = as.numeric( output[seq(2,length(output),2)] ))
    
  }
  
}

#sgn_out = summary(analysis)$coefficients[nrow(summary(analysis)$coefficients),3]
#if (is.nan(sgn_out)){
#   sgn_out = 0;
# } else if (sgn_out>0){
#   sgn_out = 1
# } else if (sgn_out<0){
#   sgn_out = -1
# }

