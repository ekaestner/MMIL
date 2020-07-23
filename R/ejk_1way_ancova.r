ejk_1way_ancova <- function( dep_var_loc,
                            grp_loc,
                            cov_loc,
                            out_put_loc,
                            alternative='two.sided') {
  
  # https://www.datanovia.com/en/lessons/ancova-in-r/
  
  #  dep_var_loc = '/home/ekaestner/Downloads/ancova_tot/dep_var.mat'
  #  grp_loc     = '/home/ekaestner/Downloads/ancova_tot/grp_var.mat'
  #  cov_loc     = '/home/ekaestner/Downloads/ancova_tot/cov_var.mat'
  #  out_put_loc = '/home/ekaestner/Downloads/ancova_tot/'  
  
  library( R.matlab )
  library( rstatix )
  library( ggpubr )
  library( dplyr )
  library( emmeans )
  
  ## Load Data
  # Dependent Variable
  dep_var = readMat(dep_var_loc) 
  dep_var_nme = names( dep_var$dep.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var_nme = dep_var_nme[-1]
  
  # Group Variable
  grp_var = readMat(grp_loc) 
  grp_var_nme = names( grp_var$grp.var[,,1] )
  grp_var = as.data.frame( lapply(grp_var$grp.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=grp_var_nme )
  for (iG in 2:ncol(grp_var))
  {
    grp_var[,iG] = factor( grp_var[,iG], exclude='N/A')
  }
  grp_var_nme = grp_var_nme[-1]

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
  
  plt_dta = merge( dep_var, grp_var, by="sbj.nme" );
  plt_dta = merge( plt_dta, cov_var, by="sbj.nme" );
  
  dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
  
  for ( iG in 1:length(grp_var_nme) )
  {
    
    dir.create( paste(out_put_loc,'/',grp_var_nme[iG],sep=''), showWarnings = FALSE)
    dir.create( paste(out_put_loc,'/',grp_var_nme[iG],'/','tests','/',sep=''), showWarnings = FALSE)
    
    # Put together table
    pvl_tbl = data.frame( DV=character(length(dep_var_nme)),
                          statistic=double(length(dep_var_nme)),
                          pvalue=double(length(dep_var_nme)),
                          df_num=double(length(dep_var_nme)),
                          df_den=double(length(dep_var_nme)),
                          report=character(length(dep_var_nme)),
                          par_eta_sq=double(length(dep_var_nme)),
                          stringsAsFactors=FALSE)
    
    # Pairwise Setup
    # Check ejk_1way_anova
    
    # Iterate through Dependent Variables
    for ( iV in 1:length(dep_var_nme) )
    {
      
      ## Put together 
      equ_txt = paste( dep_var_nme[iV], ' ~ ', sep='')
      for (iCV in 1:length(cov_var_nme)){ equ_txt = paste( equ_txt, cov_var_nme[iCV], ' + ', sep='') }
      equ_txt = paste( equ_txt, grp_var_nme[iG], sep='')
      equ_txt = as.formula( equ_txt )
      
      use_dta = as.data.frame(list(plt_dta[['sbj.nme']], plt_dta[[grp_var_nme[iG]]], plt_dta[[dep_var_nme[iV]]]), col.names=c('sbj.nme', grp_var_nme[iG], dep_var_nme[iV]))
      use_dta = merge( use_dta, cov_var, by="sbj.nme" );
      
      ## Leven Test
      # Check ejk_1way_anova
      
      ## Shapiro Test
      # Check ejk_1way_anova
      
      ## Run 1-way anova
      ttt = use_dta %>% anova_test( equ_txt, effect.size='pes' )
      
      ## Put together table
      pvl_tbl$DV[iV]          = dep_var_nme[iV]
      pvl_tbl$pvalue[iV]      = signif(ttt$p[length(ttt$p)],digits=3)
      pvl_tbl$statistic[iV]   = ttt$F[length(ttt$p)]
      pvl_tbl$df_num[iV]      = ttt$DFn[length(ttt$p)]
      pvl_tbl$df_den[iV]      = ttt$DFd[length(ttt$p)]
     
      pvl_tbl$report[iV]      = paste('F(',signif(ttt$DFn,digits=3),',',signif(ttt$DFd,digits=3),')=',signif(ttt$F,digits=3),', p=',signif(ttt$p,digits=2),sep='')
      
      pvl_tbl$par_eta_sq[iV]  = signif(ttt$pes[length(ttt$p)],digits=3)
    
      ## Estimated Means
      equ_est = as.formula( paste( dep_var_nme[iV], ' ~ ', grp_var_nme[iG], sep='') )
      est_mdl = lm( equ_txt, data=use_dta )
      emmeans( est_mdl, 'Site' )
      #with(use_dta, tapply( rhs.rostral.STG, Site, mean) )
      
      ## Make Visual Reports
      box_plt <- ggboxplot( use_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
                            ylab = dep_var_nme[iV], xlab = "Groups", add = "jitter" )
      ggsave(paste(out_put_loc,'/',grp_var_nme[iG],'/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)
      
      ## Put together pairwise
      # Check ejk_1way_anova
      
      ## Make Reports
      # Check ejk_1way_anova
      
    }
    
    write.csv( pvl_tbl, paste(out_put_loc,'/',grp_var_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
    cov_var_nme
  }
  
}

# sum(p.adjust( pvl_tbl$pvalue, method = 'fdr' )<.0
