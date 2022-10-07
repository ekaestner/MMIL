ejk_partial_correlation <- function( dta_one_loc,
                                     dta_two_loc,
                                     dta_thr_loc,
                                     out_put_loc,
                                     cor_typ ) {
  
  #dta_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/partial_correlation__surgery_pst_cog_dti/ltle_atl//dta_one.mat'
  #dta_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/partial_correlation__surgery_pst_cog_dti/ltle_atl//dta_two.mat'
  #dta_thr_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/partial_correlation__surgery_pst_cog_dti/ltle_atl/dta_thr.mat'
  #out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/partial_correlation__surgery_pst_cog_dti/ltle_atl/'
  #cor_typ     = 'spearman'
  
  library( 'R.matlab' )
  library( 'ppcor' )
  
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
  
  dta_thr = readMat(dta_thr_loc) 
  dta_thr_nme = names( dta_thr$dta.thr[,,1] )
  dta_thr = as.data.frame( lapply(dta_thr$dta.thr, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dta_thr_nme )
  dta_thr_nme = dta_thr_nme[-1]
  
  dta_use = merge( merge( dta_one, dta_two, by="sbj.nme" ), dta_thr, by="sbj.nme" )
  dta_use = dta_use[,-1]  
  write.csv(dta_use, paste( out_put_loc,'/', 'data_use.csv', sep='') )
  
  # Partial-Correlation  ####################################################################################
  for (iO in 1:length(dta_one_nme)){
    
    dir.create( paste(out_put_loc,'/',dta_one_nme[iO],sep=''), showWarnings = FALSE)
    dir.create( paste(out_put_loc,'/',dta_one_nme[iO],'/','semipartial',sep=''), showWarnings = FALSE)
    dir.create( paste(out_put_loc,'/',dta_one_nme[iO],'/','partial',sep=''), showWarnings = FALSE)
    
    rvl_out = data.frame(matrix(nrow=length(dta_two_nme), ncol=length(dta_thr_nme)+1))
      rownames(rvl_out) = dta_two_nme
      colnames(rvl_out) = c('bivariate',dta_thr_nme)
    pvl_out = data.frame(matrix(nrow=length(dta_two_nme), ncol=length(dta_thr_nme)+1))
      rownames(pvl_out) = dta_two_nme
      colnames(pvl_out) = c('bivariate',dta_thr_nme)
    num_out = data.frame(matrix(nrow=length(dta_two_nme), ncol=length(dta_thr_nme)+1))
      rownames(num_out) = dta_two_nme
      colnames(num_out) = c('bivariate',dta_thr_nme)
    
    for (iC in 1:length(dta_two_nme)){

      cor_hld = cor.test( dta_use[[dta_one_nme[iO]]], dta_use[[dta_two_nme[iC]]], method=cor_typ)
      rvl_out[iC,1] = cor_hld$estimate
      pvl_out[iC,1] = cor_hld$p.value
      num_out[iC,1] = cor_hld$n
      
      spc_dta = cbind( as.data.frame(dta_use[,dta_one_nme[iO]]), as.data.frame(dta_use[,dta_two_nme[iC]]), dta_use[,dta_thr_nme])
      colnames(spc_dta)[1:2] = c(dta_one_nme[iO],dta_two_nme[iC])
      
      pcr_hld = pcor(spc_dta[complete.cases(spc_dta),], method=cor_typ)
      write.csv(pcr_hld$estimate, paste( out_put_loc,'/',dta_one_nme[iO],'/','partial', '/', 'partial_correlation_', dta_two_nme[iC],'_rvalues.csv', sep='') )
      write.csv(pcr_hld$p.value, paste( out_put_loc,'/',dta_one_nme[iO],'/','partial', '/', 'partial_correlation_', dta_two_nme[iC],'_pvalues.csv', sep='') )
      write.csv(pcr_hld$n, paste( out_put_loc,'/',dta_one_nme[iO],'/','partial', '/', 'partial_correlation_', dta_two_nme[iC],'_num.csv', sep='') )
      
      spc_hld = spcor(spc_dta[complete.cases(spc_dta),], method=cor_typ)
      write.csv(spc_hld$estimate, paste( out_put_loc,'/',dta_one_nme[iO],'/','semipartial', '/', 'semi_partial_correlation_', dta_two_nme[iC],'_rvalues.csv', sep='') )
      write.csv(spc_hld$p.value, paste( out_put_loc,'/',dta_one_nme[iO],'/','semipartial', '/', 'semi_partial_correlation_', dta_two_nme[iC],'_pvalues.csv', sep='') )
      write.csv(spc_hld$n, paste( out_put_loc,'/',dta_one_nme[iO],'/','semipartial', '/', 'semi_partial_correlation_', dta_two_nme[iC],'_num.csv', sep='') )
      
      for (iR in 1:length(dta_thr_nme)){
        
            tst_dta = cbind(dta_use[[dta_one_nme[iO]]], dta_use[[dta_two_nme[iC]]], dta_use[[dta_thr_nme[iR]]])
            colnames(tst_dta) = c(dta_one_nme[iO],dta_two_nme[iC],dta_thr_nme[iR])
            tst_dta_cmp = data.frame(tst_dta[ complete.cases(tst_dta), ])
            
            prt_hld = pcor.test( tst_dta_cmp[[dta_one_nme[iO]]], tst_dta_cmp[[dta_two_nme[iC]]], tst_dta_cmp[[dta_thr_nme[iR]]], method=cor_typ)
            
            rvl_out[iC,iR+1] = prt_hld$estimate
            pvl_out[iC,iR+1] = prt_hld$p.value
            num_out[iC,iR+1] = prt_hld$n
      }
    }
      
      write.csv(rvl_out, paste( out_put_loc,'/',dta_one_nme[iO], '/', 'partial_correlation_rvalues.csv', sep='') )
      write.csv(pvl_out, paste( out_put_loc,'/',dta_one_nme[iO], '/', 'partial_correlation_pvalues.csv', sep='') )
      write.csv(num_out, paste( out_put_loc,'/',dta_one_nme[iO], '/', 'partial_correlation_num.csv', sep='') )
  }
  
}



