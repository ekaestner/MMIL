ejk_ttest2 <- function( dep_var_loc,
                        grp_loc,
                        out_put_loc,
                        alternative='two.sided') {

#dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/TotalCor/Cognitive/Post_TLE_ttest/dep_var.mat'
#grp_loc     = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/TotalCor/Cognitive/Post_TLE_ttest/grp_var.mat'
#out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/TotalCor/Cognitive/Post_TLE_ttest'  
    
library( R.matlab )
library( rstatix )
library( ggpubr )
library( dplyr )
  
dep_var = readMat(dep_var_loc) 
dep_var_nme = names( dep_var$dep.var[,,1] )
dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
dep_var_nme = dep_var_nme[-1]

grp_var = readMat(grp_loc) 
grp_var_nme = names( grp_var$grp.var[,,1] )
grp_var = as.data.frame( lapply(grp_var$grp.var, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=grp_var_nme )
for (iG in 2:ncol(grp_var))
{
  grp_var[,iG] = factor( grp_var[,iG], exclude='N/A')
}
grp_var_nme = grp_var_nme[2:ncol(grp_var)]

plt_dta = merge( dep_var, grp_var, by="sbj.nme" );

dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)

for ( iG in 1:length(grp_var_nme) )
{
  
  dir.create( paste(out_put_loc,'/',grp_var_nme[iG],sep=''), showWarnings = FALSE)
  
  pvl_tbl = data.frame( DV=character(length(dep_var_nme)),
                        statistic=double(length(dep_var_nme)),
                        pvalue=double(length(dep_var_nme)),
                        df=double(length(dep_var_nme)),
                        method=character(length(dep_var_nme)),
                        alternative=character(length(dep_var_nme)), 
                        report=character(length(dep_var_nme)),
                        cohensd=double(length(dep_var_nme)),
                        levene=double(length(dep_var_nme)),
                        shapiro1=double(length(dep_var_nme)),
                        shapiro2=double(length(dep_var_nme)),
                        stringsAsFactors=FALSE)
  
  for ( iV in 1:length(dep_var_nme) )
  {
    
    equ_txt = as.formula( paste( dep_var_nme[iV], ' ~ ', grp_var_nme[iG], sep='') )
    
    ## Test assumptions
    fac_lvl = levels(plt_dta[[grp_var_nme[iG]]])
    one_shp = shapiro.test( plt_dta[[dep_var_nme[iV]]][which(plt_dta[[grp_var_nme[iG]]]==fac_lvl[1])] )
    two_shp = shapiro.test( plt_dta[[dep_var_nme[iV]]][which(plt_dta[[grp_var_nme[iG]]]==fac_lvl[2])] )
    
    lev_pvl = levene_test( plt_dta, equ_txt )
    if ( lev_pvl$p <.05 ){ var_typ = FALSE } else{var_typ=TRUE}
   
    
    ## Run t-test
    ttt = t.test( plt_dta[[dep_var_nme[iV]]] ~ plt_dta[[grp_var_nme[iG]]], 
                  alternative = 'two.sided', 
                  var.equal = var_typ )
    
    ## Run effect size
    chd_val = cohens_d( plt_dta, equ_txt, var.equal = var_typ)
    
    ## Update Output Table
    pvl_tbl$DV[iV]          = dep_var_nme[iV]
    pvl_tbl$pvalue[iV]      = ttt$p.value
    pvl_tbl$statistic[iV]   = ttt$statistic
    pvl_tbl$df[iV]          = ttt$parameter
    pvl_tbl$alternative[iV] = ttt$alternative
    pvl_tbl$method[iV]      = ttt$method
    
    pvl_tbl$report[iV]      = paste('t(',signif(ttt$parameter,digits=3),')=',signif(ttt$statistic,digits=3),', p=',signif(ttt$p.value,digits=2),sep='')
    
    pvl_tbl$cohensd[iV]     = chd_val$effsize
    
    pvl_tbl$levene[iV]      = lev_pvl$p
    pvl_tbl$shapiro1[iV]    = one_shp$p.value
    pvl_tbl$shapiro2[iV]    = two_shp$p.value
    
    ## Make reports
    #bxp <- ggboxplot( plt_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
    #                  ylab = dep_var_nme[iV], xlab = "Groups", add = "jitter" )
    box_plt = ggboxplot( plt_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
                      ylab = dep_var_nme[iV], xlab = "Groups", add = "dotplot" )
    qqq_plt = ggqqplot( plt_dta, x = dep_var_nme[iV], facet.by = grp_var_nme[iG])
    out_plt = ggarrange( box_plt, qqq_plt, ncol = 1, nrow = 2)
    
    ggsave(paste(out_put_loc,'/',grp_var_nme[iG],'/',dep_var_nme[iV],'.jpg',sep=''), plot=out_plt)

  }
  
  write.csv( pvl_tbl, paste(out_put_loc,'/',grp_var_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
  
}

}