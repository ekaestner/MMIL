ejk_cross_cor_robust <- function( dta_one_loc,
                           dta_two_loc,
                           out_put_loc,
                           pvl_cut, 
                           pvl_lib ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  dta_one_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Imaging/slh_atl_mem_v2/out/Correlations/statistics/robustcorrelation__sde_non_qc_QC/lft/dta_one.mat'
  dta_two_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Imaging/slh_atl_mem_v2/out/Correlations/statistics/robustcorrelation__sde_non_qc_QC/lft/dta_two.mat'
  out_put_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Imaging/slh_atl_mem_v2/out/Correlations/statistics/robustcorrelation__sde_non_qc_QC/lft/'
  #pvl_cut = .05
  #pvl_lib = .10 

  library( 'R.matlab' )  
  library('Hmisc')
  library('dplyr')
  library('pracma')
  library('robustbase')
  
  # LOAD  ####################################################################################
  dta_one = readMat(dta_one_loc) 
  dta_one_nme = names( dta_one$dta.one[,,1] )
  dta_one = as.data.frame( lapply(dta_one$dta.one, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dta_one_nme )
  dta_one_nme = dta_one_nme[-1]
  
  dta_two = readMat(dta_two_loc) 
  dta_two_nme = names( dta_two$dta.two[,,1] )
  dta_two = as.data.frame( lapply(dta_two$dta.two, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dta_two_nme )
  dta_two_nme = dta_two_nme[-1]

  row_nme = colnames(dta_one)
  col_nme = colnames(dta_two)
  
  
  if (!(sum(!is.na(match( dta_one_nme, dta_two_nme)))==length(dta_one_nme))){
    dta_use = merge( dta_one, dta_two, by="sbj.nme" )
  } else {
    dta_use = dta_one
  }
  dta_use = dta_use[,-1]  
  
  # Cross-Correlation  ####################################################################################
  row_nme = colnames(dta_one)
  row_nme = row_nme[-1]
  col_nme = colnames(dta_two)
  col_nme = col_nme[-1]
  
  crs_cor_rvl <- matrix(nrow=length(row_nme),ncol=length(col_nme),dimnames=list(row_nme,col_nme))
  crs_cor_pvl <- matrix(nrow=length(row_nme),ncol=length(col_nme),dimnames=list(row_nme,col_nme))
  crs_cor_num <- matrix(nrow=length(row_nme),ncol=length(col_nme),dimnames=list(row_nme,col_nme))
  
  for (iC in 1:length(col_nme)){
    for (iR in 1:length(row_nme)){
      
      if (iC != iR){
      
      equ_txt = as.formula( paste( col_nme[iC], ' ~ ', row_nme[iR], sep='') )
      hld = lmrob( equ_txt, data=dta_use, method = "MM", setting="KS2014")

      hld_sum = summary(hld)
      
      crs_cor_rvl[iR,iC] = hld$coefficients[2]
      if (nrow(hld_sum$coefficients)==2){crs_cor_pvl[iR,iC] = hld_sum$coefficients[2,4]}
      crs_cor_num[iR,iC] = nrow(dta_use)
      }
    }
  }
  
  # TABLE OUTPUT  ####################################################################################
  write.csv( crs_cor_rvl, paste( out_put_loc, '/', 'cross_correlation_rvalues.csv', sep='') )
  write.csv( crs_cor_pvl, paste( out_put_loc, '/', 'cross_correlation_pvalues.csv', sep='') )
  write.csv( crs_cor_num, paste( out_put_loc, '/', 'cross_correlation_n.csv', sep='') )

  cor_rvl_nrm <- as.data.frame(as.table(crs_cor_rvl))
  
  # Find Significant r-values  ####################################################################################
  #turn into a 3-column table
  cor_rvl_nrm <- as.data.frame(as.table(crs_cor_rvl))
  cor_pvl_nrm <- as.data.frame(as.table(crs_cor_pvl))
  cor_num_nrm <- as.data.frame(as.table(crs_cor_num))
  cor_rvl_nrm$P = cor_pvl_nrm$Freq
  cor_rvl_nrm$n = cor_num_nrm$Freq
  
  #select significant values  
  cor_dta_nrm <- subset(cor_rvl_nrm, abs(P) < pvl_cut )
  cor_dta_nrm = cor_dta_nrm[order(-abs(cor_dta_nrm$Freq)),]
  
  cor_dta_lib <- subset(cor_rvl_nrm, abs(P) < pvl_lib ) 
  cor_dta_lib = cor_dta_lib[order(-abs(cor_dta_lib$Freq)),]
  
  #turn corr back into matrix in order to plot with corrplot
  if (dim(cor_dta_nrm)[1]!=0){
  mtx_cor_rvl <- reshape2::acast(cor_dta_nrm, Var1~Var2, value.var="Freq")
  mtx_cor_pvl <- reshape2::acast(cor_dta_nrm, Var1~Var2, value.var="P")}
  
  if (dim(cor_dta_lib)[1]!=0){
  mtx_cor_lib_rvl <- reshape2::acast(cor_dta_lib, Var1~Var2, value.var="Freq")
  mtx_cor_lib_pvl <- reshape2::acast(cor_dta_lib, Var1~Var2, value.var="P")
  }
  
  # Output  ####################################################################################
  # write out tables
  if (dim(cor_dta_nrm)[1]!=0){
    write.csv( cor_dta_nrm, paste( out_put_loc, '/', 'cross_correlation_significances.csv', sep='') )
    write.csv( mtx_cor_rvl, paste( out_put_loc, '/', 'cross_correlation_subset_rvalues.csv', sep='') )  
    write.csv( mtx_cor_pvl, paste( out_put_loc, '/', 'cross_correlation_subset_pvalues.csv', sep='') )
    
  }
  
  if (dim(cor_dta_lib)[1]!=0){
    write.csv( cor_dta_lib, paste( out_put_loc, '/', 'cross_correlation_significances_liberal.csv', sep='') )
    write.csv( mtx_cor_lib_rvl, paste( out_put_loc, '/', 'cross_correlation_subset_liberal_rvalues.csv', sep='') )
    write.csv( mtx_cor_lib_pvl, paste( out_put_loc, '/', 'cross_correlation_subset_liberal_pvalues.csv', sep='') )
    
  }
  
}



  







