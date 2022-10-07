ejk_surface_correlation_ <- function( grp_one_loc,
                                      grp_two_loc,
                                      out_put_loc ) {
  
  grp_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/surf/lm1_chg_ltle_atl_256'
  grp_one_nme = 'ltle_atl'
  grp_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/surf/lm1_chg_ltle_slah_256'
  grp_two_nme = 'ltle_slah'
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/surf/fishersZ_256_lm1'
  
  library( R.matlab )
  library( diffcor)
  library( emmeans )
  library( parallel )
  library( doParallel )
  library( foreach )
  library( stringr )
  library( broom )
  library( tidyr )
  library(jtools)
  
  ## Load Data ###################################################################################################################################################
  dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
  
  grp_one_num = read.csv(paste(grp_one_loc,"/","exact_data_","1",".csv",sep=''))
  grp_one_num = cor.test( grp_one_num[,3], grp_one_num[,5], method='pearson')
  grp_one_num = grp_one_num$parameter + 2
  
  grp_two_num = read.csv(paste(grp_two_loc,"/","exact_data_","1",".csv",sep=''))
  grp_two_num = cor.test( grp_two_num[,3], grp_two_num[,5], method='pearson')
  grp_two_num = grp_two_num$parameter + 2
  
  ## Load Data ###################################################################################################################################################
  srf_dta_nme = c('lhs_dep_var', 'rhs_dep_var')
  
  for( iH in 1:2){
  
    grp_one_rvl = readMat(paste(grp_one_loc,'/','rvalues_', srf_dta_nme[iH], '_',grp_one_nme,'.mat',sep='')) 
    grp_one_rvl = as.data.frame( lapply(grp_one_rvl$rvalues, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names='' )
    
    grp_two_rvl = readMat(paste(grp_two_loc,'/','rvalues_', srf_dta_nme[iH],'_',grp_two_nme,'.mat',sep='')) 
    grp_two_rvl = as.data.frame( lapply(grp_two_rvl$rvalues, unlist, use.names=FALSE), 
                                 fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names='' )  
    
    # Set up backend for parallelizing
    ncores <- detectCores() - 2
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    output <- foreach( iV = 1:163842, .combine = 'combine', .packages = c("diffcor","tidyr", "broom", "emmeans", "pracma", "stringr","jtools"), .multicombine=TRUE, #
                       .init=list(list(), list()) ) %dopar% {
                         
                         # Running diffcor.two()
                         tst_hld = diffcor.two(grp_one_rvl[1,iV],grp_two_rvl[1,iV],grp_one_num,grp_two_num,alternative = c("two.sided"))

                         # Putting result into one list
                         zvl <- tst_hld$z
                         pvl <- as.numeric(as.character(tst_hld$p))
                         
                         # Current progress
                         list( zvl, pvl )  
                         
                       }
    
    # Free Up CPU and Memory
    stopCluster(cl)

    # Create the file names of the variables 
    zvl_fle_nme = paste( out_put_loc,"/", "zvalues_" , grp_one_nme,'_',grp_two_nme ,'_', srf_dta_nme[iH], ".mat", sep = "")
    pvl_fle_nme = paste( out_put_loc,"/", "pvalues_" , grp_one_nme,'_',grp_two_nme ,'_', srf_dta_nme[iH], ".mat", sep = "")
    
    # Write result into a matrix for matlab
    output = unlist(output)
    
    writeMat(zvl_fle_nme, zvalues = as.numeric( output[seq(1,length(output),2)] ))
    writeMat(pvl_fle_nme, pvalues = as.numeric( output[seq(2,length(output),2)] ))
    
    
    
    
  }
}