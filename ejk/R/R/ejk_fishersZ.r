ejk_fishersZ <- function( rvl_one_loc,
                          rvl_two_loc,
                          num_one_loc,
                          num_two_loc,
                          out_put_loc,
                          out_put_pre) {
  
  #rvl_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_spearman/FishersZ//rvl_one.mat'
  #rvl_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_spearman/FishersZ//rvl_two.mat'
  #num_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_spearman/FishersZ//num_one.mat'
  #num_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_spearman/FishersZ//num_two.mat'
  #out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/correlation__surgery_pst_cog_spearman/FishersZ/'
  #out_put_pre = 'ltle_surgery'
  
  library( 'R.matlab' ) 
  library('diffcor')
  library('cocor')
  
  # LOAD  ####################################################################################
  rvl_one = readMat(rvl_one_loc) 
  rvl_one_nme = names( rvl_one$rvl.one[,,1] )
  rvl_one = as.data.frame( lapply(rvl_one$rvl.one, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=rvl_one_nme )

  rvl_two = readMat(rvl_two_loc) 
  rvl_two_nme = names( rvl_two$rvl.two[,,1] )
  rvl_two = as.data.frame( lapply(rvl_two$rvl.two, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=rvl_two_nme )
  
  num_one = readMat(num_one_loc) 
  num_one_nme = names( num_one$num.one[,,1] )
  num_one = as.data.frame( lapply(num_one$num.one, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=num_one_nme )
  
  num_two = readMat(num_two_loc) 
  num_two_nme = names( num_two$num.two[,,1] )
  num_two = as.data.frame( lapply(num_two$num.two, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=num_two_nme )
  
  # Fisher's Z  ####################################################################################
  fsh_pvl_out = matrix(data = NaN, nrow = nrow(rvl_one), ncol = ncol(rvl_one));
  zou_int_out = data.frame(matrix(data = "", nrow = nrow(rvl_one), ncol = ncol(rvl_one)),stringsAsFactors=FALSE)
  for (iC in 1:ncol(rvl_one)){
    for (iR in 1:nrow(rvl_one)){
      tst_hld = diffcor.two(rvl_one[iR,iC],rvl_two[iR,iC],num_one[iR,iC],num_two[iR,iC],alternative = c("two.sided"))
      fsh_pvl_out[iR,iC] = as.numeric(as.character(tst_hld$p))
      
      tst_hld = cocor.indep.groups(rvl_one[iR,iC],rvl_two[iR,iC],num_one[iR,iC],num_two[iR,iC],alternative = c("two.sided"))
      zou_int_out[iR,iC] = paste('[',signif(tst_hld@zou2007$conf.int[1],digits=3),' ',signif(tst_hld@zou2007$conf.int[2],digits=3),']',sep='')
    }
  }
  
  # TABLE OUTPUT  ####################################################################################
  write.csv( fsh_pvl_out, paste( out_put_loc, '/', out_put_pre, '_fishersZ.csv', sep='') )
  write.csv( zou_int_out, paste( out_put_loc, '/', out_put_pre, '_Zou2007.csv', sep='') )
  
}