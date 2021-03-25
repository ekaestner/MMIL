ejk_surface_1way_ancova <- function( lhs_srf_var_loc,
                                     rhs_srf_var_loc,
                                     grp_loc,
                                     cov_loc,
                                     iG,
                                     fst_nme,
                                     scd_nme,
                                     out_put_loc ) {
  
#  lhs_srf_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat'
#  rhs_srf_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/surf_aMRI_thickness_rhs_sm176_no_nan_no_mass.mat'
#  grp_loc     = '/home/ekaestner/Downloads/surface_ancova/grp_var.mat'
#  cov_loc     = '/home/ekaestner/Downloads/surface_ancova/cov_var.mat'
#  iG          = 1
#  fst_nme     = 'TLE'
#  scd_nme     = 'HC'
#  out_put_loc = '/home/ekaestner/Downloads/surface_ancova/'  
  
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

  ## Load Data ###################################################################################################################################################
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

  # Covariates
  cov_var = readMat(cov_loc) 
  cov_var_nme = names( cov_var$cov.var[,,1] )
  cov_var = as.data.frame( lapply(cov_var$cov.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=cov_var_nme )
  for (iCV in 2:ncol(cov_var))
  {
    if ( is.character(cov_var[1,iCV]) ) {
    cov_var[,iCV] = factor( cov_var[,iCV], exclude='N/A')
    }
  }
  cov_var_nme = cov_var_nme[-1]
  num_cov = length(cov_var_nme)
  
## ANCOVA #################################################################################################################################################
  srf_dta_nme = c('lhs_dep_var', 'rhs_dep_var')
  
  for( iH in 1:2){
  
    # Dependent Variable
    if (iH==1 ){ 
      srf_dta = readMat(lhs_srf_var_loc)
      dep_var = t(as.data.frame( lapply(srf_dta$srf.dta.sbj, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE ))
      srf_dta = srf_dta$srf.dta }
    else if (iH==2){ 
      srf_dta = readMat(rhs_srf_var_loc) 
      dep_var = t(as.data.frame( lapply(srf_dta$srf.dta.sbj, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE ))
      srf_dta = srf_dta$srf.dta }

    colnames(dep_var) = c('sbj.nme')
    
    # Merge
    plt_dta = merge( dep_var, grp_var, by="sbj.nme" );
    plt_dta = merge( plt_dta, cov_var, by="sbj.nme" );

    dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
    
    # Set up backend for parallelizing
    ncores <- detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    # Setup 
    equ_txt = paste( srf_dta_nme[iH], ' ~ ', sep='')
    for (iCV in 1:length(cov_var_nme)){ equ_txt = paste( equ_txt, cov_var_nme[iCV], ' + ', sep='') }
    equ_txt = paste( equ_txt, grp_var_nme[iG], sep='')
    equ_txt = as.formula( equ_txt )
    
    # Data merge
    use_dta = merge( dep_var, grp_var, by="sbj.nme" );
    use_dta = merge( use_dta, cov_var, by="sbj.nme" );
    
    # Run test for each vertex individually
    output <- foreach( iV = 1:163842, .combine = 'combine', .packages = c("tidyr", "broom", "emmeans", "pracma", "stringr"), .multicombine=TRUE, #
                       .init=list(list(), list(), list()) ) %dopar% {
                         
                         # Put together data
                         use_dta_col = as.data.frame(srf_dta[, iV], fix.empty.names=FALSE )
                         colnames(use_dta_col) = c(srf_dta_nme[iH])
                         use_dta_col = cbind( dep_var, use_dta_col )
                         use_dta_par = merge( use_dta, use_dta_col, by="sbj.nme" )
                         
                         # Running ANCOVA (Model: factor(sex) + factor(strength) + age + factor(group))
                         analysis <- aov( equ_txt, data=use_dta_par )
                         
                         ## Estimated Means
                         emm_hld = tidy(emmeans(analysis, grp_var_nme[iG], data=use_dta_par) )
                         F_mean <- emm_hld[which(str_detect(as.character(emm_hld[[grp_var_nme[iG]]]), fst_nme)), 'estimate']
                         S_mean <- emm_hld[which(str_detect(as.character(emm_hld[[grp_var_nme[iG]]]), scd_nme)), 'estimate']
                         
                         # Putting result into one list
                         emmean <- F_mean - S_mean
                         fValue <- summary(analysis)[[1]][[grp_var_nme[iG], "F value"]]
                         pValue <- summary(analysis)[[1]][[grp_var_nme[iG], "Pr(>F)"]]
                         
                         # Current progress
                         list(emmean, fValue, pValue)
                         
                       }
    
    # Free Up CPU and Memory
    stopCluster(cl)
    
    rm('srf_dta')
    
    # Create the file names of the variables 
    eFile <- paste( out_put_loc, "emmean_" , srf_dta_nme[iH], "_", fst_nme, "_", scd_nme, ".mat", sep = "")
    fFile <- paste( out_put_loc, "fValue_" , srf_dta_nme[iH], "_", fst_nme, "_", scd_nme, ".mat", sep = "")
    pFile <- paste( out_put_loc, "pValue_" , srf_dta_nme[iH], "_", fst_nme, "_", scd_nme, ".mat", sep = "")
    
    # Write result into a matrix for matlab
    output = unlist(output)

    writeMat(eFile, emmean = as.numeric( output[seq(1,length(output),3)] ))
    writeMat(fFile, fValue = as.numeric( output[seq(2,length(output),3)] ))
    writeMat(pFile, pValue = as.numeric( output[seq(3,length(output),3)] ))
    
  }

}
