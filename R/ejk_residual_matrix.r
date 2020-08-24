ejk_residual_matrix <- function( dta_var_loc,
                             cov_loc,
                             out_put_loc ) {
  
  #dta_var_loc = '/home/ekaestner/Downloads/Braintest/residuals/dta_var.mat'
  #cov_loc     = '/home/ekaestner/Downloads/Braintest/residuals/cov_var.mat'
  #out_put_loc = '/home/ekaestner/Downloads/Braintest/residuals/'  
  
  library( R.matlab )
  
  ## Load Data ###########################
  # Dependent Variable
  dta_var = readMat(dta_var_loc) 
  dta_var_nme = names( dta_var$dta.var[,,1] )
  dta_var = as.data.frame( lapply(dta_var$dta.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dta_var_nme )
  dta_var_nme = dta_var_nme[-1]
  
  # Covariates
  cov_var = readMat(cov_loc) 
  cov_var_nme = names( cov_var$cov.var[,,1] )
  cov_var = as.data.frame( lapply(cov_var$cov.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=cov_var_nme )
  for (iG in 2:ncol(cov_var))
  {
    if ( is.character(cov_var[1,iG]) ) {
      cov_var[,iG] = factor( cov_var[,iG], exclude='N/A')
    }
  }
  cov_var_nme = cov_var_nme[-1]
  num_cov = length(cov_var_nme)

  plt_dta = merge( dta_var, cov_var, by="sbj.nme" );
  
  dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
  
  ## Initiatlize ###########################
  rsd_tbl_out = matrix( vector(), nrow(dta_var), ncol(dta_var)-1 )
  
  ##  Residuals #################################################################
  for ( iV in 1:length(dta_var_nme) )
  {
  
  equ_txt = paste( dta_var_nme[iV], ' ~ ', sep='')
  for (iCV in 1:length(cov_var_nme)){ equ_txt = paste( equ_txt, cov_var_nme[iCV], ' + ', sep='') }
  equ_txt = substr(equ_txt,1,nchar(equ_txt)-3)
  equ_txt = as.formula( equ_txt )
  
  rsd_mdl = lm( equ_txt, data=plt_dta)
  
  rsd_hld = resid(rsd_mdl)
  #rsd_hld = as.data.frame(rsd_hld)
  rsd_tbl_out[,iV] = as.matrix(rsd_hld)
  }
  
  rsd_dta_frm_out = as.data.frame( rsd_tbl_out )
  colnames(rsd_dta_frm_out) = dta_var_nme
  rownames(rsd_dta_frm_out) = dta_var[,1]
  
  write.csv( rsd_dta_frm_out, paste( out_put_loc,'/','/','residual_matrix.csv',sep='') )
  
}

#iV = 5
#
#rsd_plt_var = data.frame( org_dta=double(nrow(dta_var)),
#                          rsd_dta=double(nrow(dta_var)),
#                          age_dta=double(nrow(dta_var)) )
#
#rsd_plt_var$org_dta = dta_var[[ dta_var_nme[iV] ]]
#rsd_plt_var$rsd_dta = rsd_dta_frm_out[[ dta_var_nme[iV] ]]
#rsd_plt_var$age_dta = cov_var$Age
#
#dta_one_plt = ggscatter( rsd_plt_var, 'rsd_dta', 'org_dta'  )
#dta_two_plt = ggscatter( rsd_plt_var, 'rsd_dta', 'age_dta'  )
#dta_thr_plt = ggscatter( rsd_plt_var, 'org_dta', 'age_dta'  )
#
#ggarrange( dta_one_plt, dta_two_plt, dta_thr_plt, ncol = 2, nrow = 2)


