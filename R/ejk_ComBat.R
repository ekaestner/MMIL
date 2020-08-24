ejk_ComBat <- function( dta_var_loc,
                        btc_var_loc,
                        cov_var_loc,
                        out_put_loc ) {

  #dta_var_loc = '/home/ekaestner/Downloads/ComBat/dta_var.mat'
  #btc_var_loc = '/home/ekaestner/Downloads/ComBat/btc_var.mat'
  #cov_var_loc = '/home/ekaestner/Downloads/ComBat/cov_var.mat'
  #out_put_loc = '/home/ekaestner/Downloads/ComBat/'
  
  library( R.matlab )  
  library( neuroCombat )

  ## Load Data ################################
  # Dependent Variable
  dep_var = readMat(dta_var_loc) 
  dep_var_nme = names( dep_var$dta.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dta.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var     = dep_var[-1]
  dep_var_nme = dep_var_nme[-1]
  
  # Group Variable
  btc_var = readMat(btc_var_loc) 
  btc_var_nme = names( btc_var$btc.var[,,1] )
  btc_var = as.data.frame( lapply(btc_var$btc.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=btc_var_nme )
  btc_var = btc_var[-1]
  btc_var_nme = btc_var_nme[-1]
  
  # Covariates (optional)
  if ( cov_var_loc != '' ){
  cov_var = readMat(cov_var_loc) 
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
  }
  
  ## ComBat Data ################################
  if ( cov_var_loc != '' ){
    
    equ_txt = ' ~ '
    for (iCV in 1:length(cov_var_nme)){ 
      equ_txt = paste( equ_txt, paste(cov_var_nme[iCV],'_hld',sep=''), ' + ', sep='')
      assign( paste(cov_var_nme[iCV],'_hld',sep=''), cov_var[[cov_var_nme[iCV]]])
      }
    equ_txt = substr(equ_txt,1,nchar(equ_txt)-3)
    equ_txt = as.formula( equ_txt )
    
    cov_mdl = model.matrix( equ_txt )
    
    com_bat_dta = neuroCombat( as.matrix(t(dep_var)), as.character(btc_var[[btc_var_nme[1]]]), cov_mdl )
    out_dta = com_bat_dta$dat.combat
    
  } else {
    com_bat_dta = neuroCombat( as.matrix(t(dep_var)), as.character(btc_var[[btc_var_nme[1]]]) )
    out_dta = com_bat_dta$dat.combat
  }
  
  ## Save Data
  write.csv( out_dta, paste( out_put_loc, '/', 'out_dta.csv', sep='') )
}
