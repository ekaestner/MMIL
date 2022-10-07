ejk_cross_cor <- function( dta_one_loc,
                           dta_two_loc,
                           out_put_loc,
                           pvl_cut, 
                           pvl_lib ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  dta_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_robust/ltle_slah/dta_one.mat'
  dta_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_robust/ltle_slah//dta_two.mat'
  out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_robust/ltle_slah/'
  pvl_cut = .05
  pvl_lib = .10 

  library( 'R.matlab' )  
  library('corrplot')
  library('Hmisc')
  library('dplyr')
  library('pracma')
  library(robustbase)
  
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
  
  equ_txt = as.formula( paste( col_nme[3], ' ~ ', row_nme[1], ' + ', row_nme[10], sep='') )
  hld = lmrob( equ_txt, data=dta_use, method = "MM", setting="KS2014")
  summary(hld)
  
  equ_txt = as.formula( paste( col_nme[3], ' ~ ', row_nme[1], ' + ', row_nme[3], ' + ', row_nme[11], sep='') )
  hld = lmrob( equ_txt, data=dta_use, method = "MM", setting="KS2014")
  summary(hld)
  
  hld$coefficients[2,4]
  hld$r.squared
  hld$adj.r.squared
  hld$df[2]
  
  crs_cor_dta_frm <- rcorr( data.matrix(dta_use), type=cor_typ)
  
  row_ind = match( row_nme, rownames(crs_cor_dta_frm$r))
  col_ind = match( col_nme, colnames(crs_cor_dta_frm$r))
  
  crs_cor_rvl = crs_cor_dta_frm$r[ row_ind, col_ind ]
  crs_cor_pvl = crs_cor_dta_frm$P[ row_ind, col_ind ]
  crs_cor_num = crs_cor_dta_frm$n[ row_ind, col_ind ]
  
  # Shuffle control  ####################################################################################
  shf_cor_pvl = crs_cor_pvl
  
  if (shf_pvl==1) {
    for (iC in 1:ncol(crs_cor_pvl)){
      for (iR in 1:nrow(crs_cor_pvl)){
        if (!is.na(crs_cor_pvl[iR,iC]) & (crs_cor_pvl[iR,iC]<.05)){
          
          dta_one_tmp_cmp = dta_use[[dta_one_nme[iR]]]
          dta_two_tmp_cmp = dta_use[[dta_two_nme[iC]]]
          
          ind_use = !(is.nan(dta_one_tmp_cmp)) & !(is.nan(dta_two_tmp_cmp))
          
          dta_one_tmp_cmp = dta_one_tmp_cmp[ind_use]
          dta_two_tmp_cmp = dta_two_tmp_cmp[ind_use]
          
          rnd_hld <- numeric(shf_num)
          for (iSH in 1:shf_num){
            
            dta_one_tmp = dta_one_tmp_cmp
            dta_two_tmp = sample(dta_two_tmp_cmp, length(dta_one_tmp), replace = FALSE)
            
            rnd_hld[iSH] = cor( dta_one_tmp, dta_two_tmp, method=cor_typ, use="complete.obs" )
            
          }
          
          shf_cor_pvl[iR,iC] <- length( rnd_hld[rnd_hld >= crs_cor_rvl[iR,iC]])/shf_num
          
        } else{shf_cor_pvl[iR,iC]=NaN}
      }
    }
  }
  
  # TABLE OUTPUT  ####################################################################################
  write.csv( crs_cor_rvl, paste( out_put_loc, '/', 'cross_correlation_rvalues.csv', sep='') )
  write.csv( crs_cor_pvl, paste( out_put_loc, '/', 'cross_correlation_pvalues.csv', sep='') )
  write.csv( crs_cor_num, paste( out_put_loc, '/', 'cross_correlation_n.csv', sep='') )
  if (shf_pvl==1) { write.csv( shf_cor_pvl, paste( out_put_loc, '/', 'cross_correlation_pvalues_shuffle.csv', sep=''), na = "NaN" )}
  
  cor_rvl_nrm <- as.data.frame(as.table(crs_cor_rvl))
  
  jpeg(paste( out_put_loc, '/', 'cross_correlation_visualization.jpg', sep=''), width=1080, height=1080)
  crs_cor_rvl_plt = crs_cor_rvl
  crs_cor_rvl_plt[is.nan(crs_cor_rvl)] = NA
  crs_cor_rvl_plt[is.na(crs_cor_rvl)]  = 0
  corrplot(crs_cor_rvl_plt, is.corr=FALSE, tl.col="black", na.label=" ", cl.lim = c(-1, 1))
  dev.off()  
  
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
    
    if (!(nrow(mtx_cor_lib_rvl)==1 & ncol(mtx_cor_rvl)==1)){
    jpeg(paste( out_put_loc, '/', 'significant_correlations.jpg', sep=''), width=1080, height=1080)
    corrplot(mtx_cor_rvl, is.corr=FALSE, tl.col="black", na.label=" ")
    dev.off()}
    
  }
  
  if (dim(cor_dta_lib)[1]!=0){
    write.csv( cor_dta_lib, paste( out_put_loc, '/', 'cross_correlation_significances_liberal.csv', sep='') )
    write.csv( mtx_cor_lib_rvl, paste( out_put_loc, '/', 'cross_correlation_subset_liberal_rvalues.csv', sep='') )
    write.csv( mtx_cor_lib_pvl, paste( out_put_loc, '/', 'cross_correlation_subset_liberal_pvalues.csv', sep='') )
    
    if (!(nrow(mtx_cor_lib_rvl)==1 & ncol(mtx_cor_lib_rvl)==1)){
    jpeg(paste( out_put_loc, '/', 'significant_correlations_liberal.jpg', sep=''), width=1080, height=1080)
    corrplot(mtx_cor_lib_rvl, is.corr=FALSE, tl.col="black", na.label=" ")
    dev.off()}
  }
  
}











