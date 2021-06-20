ejk_ttest1 <- function( dep_var_loc,
                        out_put_loc,
                        mean_val = 0,
                        alternative='two.sided') {
  
  #dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Cognitive/TLE_post_ttest_left/dep_var.mat'
  #out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Cognitive/TLE_post_ttest_left/'  
  
  library( R.matlab )
  library( rstatix )
  library( ggpubr )
  library( dplyr )
  
  dep_var = readMat(dep_var_loc) 
  dep_var_nme = names( dep_var$dep.var[,,1] )
  dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                           fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
  dep_var_nme = dep_var_nme[-1]
  
  dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)

  pvl_tbl = data.frame( DV=character(length(dep_var_nme)),
                        statistic=double(length(dep_var_nme)),
                        pvalue=double(length(dep_var_nme)),
                        df=double(length(dep_var_nme)),
                        method=character(length(dep_var_nme)),
                        alternative=character(length(dep_var_nme)), 
                        report=character(length(dep_var_nme)),
                        cohensd=double(length(dep_var_nme)),
                        shapiro1=double(length(dep_var_nme)),
                        stringsAsFactors=FALSE)
  
  for ( iV in 1:length(dep_var_nme) )
  {
    
    ## Test assumptions
    one_shp = shapiro.test( dep_var[[dep_var_nme[iV]]] )

    ## Run t-test
    ttt = t.test( dep_var[[dep_var_nme[iV]]], 
                  mu=mean_val,
                  alternative = alternative )
    
    ## Run effect size
    equ_txt = as.formula( paste( dep_var_nme[iV], ' ~ ', '1', sep='') )
    
    chd_val = cohens_d( dep_var, equ_txt, mu=mean_val)
    
    ## Update Output Table
    pvl_tbl$DV[iV]          = dep_var_nme[iV]
    pvl_tbl$pvalue[iV]      = ttt$p.value
    pvl_tbl$statistic[iV]   = ttt$statistic
    pvl_tbl$df[iV]          = ttt$parameter
    pvl_tbl$alternative[iV] = ttt$alternative
    pvl_tbl$method[iV]      = ttt$method
    
    pvl_tbl$report[iV]      = paste('t(',signif(ttt$parameter,digits=3),')=',signif(ttt$statistic,digits=3),'; p=',signif(ttt$p.value,digits=2),sep='')
    
    pvl_tbl$cohensd[iV]     = chd_val$effsize
    
    pvl_tbl$shapiro1[iV]    = one_shp$p.value

    ## Make reports
    #bxp <- ggboxplot( plt_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
    #                  ylab = dep_var_nme[iV], xlab = "Groups", add = "jitter" )
    box_plt = ggboxplot( dep_var, y = dep_var_nme[iV], 
                         ylab = dep_var_nme[iV], xlab = "Groups", add = "dotplot" )

    ggsave(paste(out_put_loc,'/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)
      
    write.csv( pvl_tbl, paste(out_put_loc,'/','output_table.csv',sep=''), row.names=FALSE)
    
  }
  
}