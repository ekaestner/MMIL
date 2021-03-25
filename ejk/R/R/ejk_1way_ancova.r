ejk_1way_ancova <- function( dep_var_loc,
                            grp_loc,
                            cov_loc,
                            out_put_loc,
                            alternative='two.sided') {
  
  # https://www.datanovia.com/en/lessons/ancova-in-r/
  
  #dep_var_loc = '/home/ekaestner/Downloads/Braintest/desikan/orig/dep_var.mat'
  #grp_loc     = '/home/ekaestner/Downloads/Braintest/desikan/orig/grp_var.mat'
  #cov_loc     = '/home/ekaestner/Downloads/Braintest/desikan/orig/cov_var.mat'
  #out_put_loc = '/home/ekaestner/Downloads/Braintest/desikan/'  
  
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
  num_cov = length(cov_var_nme)
  
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
    equ_txt = as.formula( paste( dep_var_nme[1], ' ~ ', grp_var_nme[iG], sep='') )
    pwc = plt_dta %>% tukey_hsd( equ_txt )
    
    par_wse_nme = levels(plt_dta[[grp_var_nme[1]]])
    
    par_wse_num = length(par_wse_nme)
    par_wse_num = (par_wse_num*(par_wse_num-1)) / 2
    
    tuk_org_tbl = data.frame(matrix(vector(), length(dep_var_nme), par_wse_num))
    tuk_org_nme = data.frame(name=character(par_wse_num),stringsAsFactors=FALSE)
   
    for ( iC in 1:par_wse_num )
    {
      tuk_org_nme$name[iC] = paste( pwc$group1[iC], '_', pwc$group2[iC], '_tuk', sep='')
    }
    
    colnames(tuk_org_tbl) = tuk_org_nme$name

    # Estimated Means Table Setup
    est_men_hld = data.frame(matrix(vector(), length(dep_var_nme), length(levels(grp_var[[grp_var_nme[iG]]]))+1))
    colnames(est_men_hld) = c('DV',levels(grp_var[[grp_var_nme[iG]]]))
    
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
     
      pvl_tbl$report[iV]      = paste('F(',signif(ttt$DFn[num_cov+1],digits=3),',',signif(ttt$DFd[num_cov+1],digits=3),')=',signif(ttt$F[num_cov+1],digits=3),', p=',signif(ttt$p[num_cov+1],digits=2),sep='')
      
      pvl_tbl$par_eta_sq[iV]  = signif(ttt$pes[length(ttt$p)],digits=3)
    
      ## Estimated Means
      est_mdl = lm( equ_txt, data=use_dta )
      est_men_dta = emmeans( est_mdl, grp_var_nme[iG] )
      est_men_dta_frm = as.data.frame(est_men_dta)
      for ( iC in 1:nrow(est_men_dta_frm) )
      {
        est_men_hld[[ as.character(est_men_dta_frm[iC,1]) ]][iV] = paste(signif(est_men_dta_frm[iC,2],digits=3),' (',signif(est_men_dta_frm[iC,3],digits=3),')',sep='')
      }
      est_men_hld$DV[iV]          = dep_var_nme[iV]
      
      ## Put together pairwise
      pwc = as.data.frame(pairs(est_men_dta))
      for ( iC in 1:par_wse_num )
      {
        tuk_org_tbl[[tuk_org_nme$name[iC]]][iV] = signif(pwc$p.value[iC],digits=3)
      }
      
      ## Make Visual Reports
      box_plt <- ggboxplot( use_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
                            ylab = dep_var_nme[iV], xlab = "Groups", add = "jitter" )
      ggsave(paste(out_put_loc,'/',grp_var_nme[iG],'/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)
      
      ## Put together pairwise
      # Check ejk_1way_anova
      
      ## Make Reports
      # Check ejk_1way_anova
      
    }
    
    out_tbl = cbind( pvl_tbl, tuk_org_tbl) 
    write.csv( out_tbl, paste(out_put_loc,'/',grp_var_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
    write.csv( est_men_hld, paste(out_put_loc,'/',grp_var_nme[iG],'/','output_estimated_means_table.csv',sep=''), row.names=FALSE)
    
  }
  
}

# sum(p.adjust( pvl_tbl$pvalue, method = 'fdr' )<.0
