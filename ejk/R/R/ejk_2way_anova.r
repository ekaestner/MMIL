ejk_2way_anova <- function( dep_var_loc,
                            grp_loc_one,
                            grp_loc_two,
                            out_put_loc,
                            alternative='two.sided') {
  
#dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/Cognitive_2wayANOVA_surgery_pst_cog/dep_var.mat'
#grp_loc_one = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/Cognitive_2wayANOVA_surgery_pst_cog/grp_var_one.mat'
#grp_loc_two = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/Cognitive_2wayANOVA_surgery_pst_cog/grp_var_two.mat'
#out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis//stats/Cognitive_2wayANOVA_surgery_pst_cog/'  

library( R.matlab )
library( rstatix )
library( ggpubr )
library( dplyr )

dep_var = readMat(dep_var_loc) 
dep_var_nme = names( dep_var$dep.var[,,1] )
dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
dep_var_nme = dep_var_nme[-1]

grp_var_one = readMat(grp_loc_one) 
grp_var_one_nme = names( grp_var_one$grp.var.one[,,1] )
grp_var_one = as.data.frame( lapply(grp_var_one$grp.var.one, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=grp_var_one_nme )
for (iG in 2:ncol(grp_var_one))
{
  grp_var_one[,iG] = factor( grp_var_one[,iG], exclude='N/A')
}
grp_var_one_nme = grp_var_one_nme[-1]

grp_var_two = readMat(grp_loc_two) 
grp_var_two_nme = names( grp_var_two$grp.var.two[,,1] )
grp_var_two = as.data.frame( lapply(grp_var_two$grp.var.two, unlist, use.names=FALSE), 
                             fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=grp_var_two_nme )
for (iG in 2:ncol(grp_var_two))
{
  grp_var_two[,iG] = factor( grp_var_two[,iG], exclude='N/A')
}
grp_var_two_nme = grp_var_two_nme[-1]

plt_dta = merge(merge( dep_var, grp_var_one, by="sbj.nme" ),grp_var_two, by="sbj.nme");

dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)

for ( iG in 1:length(grp_var_one_nme) )
{
  
  dir.create( paste(out_put_loc,'/',grp_var_one_nme[iG],sep=''), showWarnings = FALSE)
  dir.create( paste(out_put_loc,'/',grp_var_one_nme[iG],'/','tests','/',sep=''), showWarnings = FALSE)
  
  # Put together table
  pvl_tbl = data.frame( DV=character(length(dep_var_nme)),
                        statistic_one=double(length(dep_var_nme)),
                        statistic_two=double(length(dep_var_nme)),
                        statistic_int=double(length(dep_var_nme)),
                        pvalue_one=double(length(dep_var_nme)),
                        pvalue_two=double(length(dep_var_nme)),
                        pvalue_int=double(length(dep_var_nme)),
                        df_num_one=double(length(dep_var_nme)),
                        df_num_two=double(length(dep_var_nme)),
                        df_num_int=double(length(dep_var_nme)),
                        df_den_one=double(length(dep_var_nme)),
                        df_den_two=double(length(dep_var_nme)),
                        df_den_int=double(length(dep_var_nme)),
                        #method=character(length(dep_var_nme)),
                        report_one=character(length(dep_var_nme)),
                        report_two=character(length(dep_var_nme)),
                        report_int=character(length(dep_var_nme)),
                        par_eta_sq_one=double(length(dep_var_nme)),
                        par_eta_sq_two=double(length(dep_var_nme)),
                        par_eta_sq_int=double(length(dep_var_nme)),
                        levene=double(length(dep_var_nme)),
                        shapiro=double(length(dep_var_nme)),
                        stringsAsFactors=FALSE)
  
  # Put together table for ttests
  equ_txt = as.formula( paste( dep_var_nme[1], ' ~ ', grp_var_one_nme[iG],'*',grp_var_two_nme[iG], sep='') )
  pwc = plt_dta %>% tukey_hsd( equ_txt )
  
  par_wse_nme = paste(pwc$group1,'_VS_',pwc$group2,sep='')
  
  tuk_org_tbl = data.frame(matrix(vector(), length(dep_var_nme), length(par_wse_nme)))
    tuk_org_nme =par_wse_nme
    colnames(tuk_org_tbl) = par_wse_nme
  
  # Iterate through Dependent Variables
  for ( iV in 1:length(dep_var_nme) )
  {
   
    equ_txt = as.formula( paste( dep_var_nme[iV], ' ~ ', grp_var_one_nme[iG],'*',grp_var_two_nme[iG], sep='') )
    
    use_dta = as.data.frame(list(plt_dta[[grp_var_one_nme[iG]]],plt_dta[[grp_var_two_nme[iG]]], plt_dta[[dep_var_nme[iV]]]), col.names=c(grp_var_one_nme[iG],grp_var_two_nme[iG],dep_var_nme[iV]))
    
    ## Leven Test
    lev_tst = use_dta %>% levene_test(equ_txt)
    
    ## Shapiro Test
    nor_mdl = lm(equ_txt, data = use_dta)
    shp_tst = shapiro_test(residuals(nor_mdl))
    
    ## Run 1-way anova
    ttt = use_dta %>% anova_test( equ_txt, effect.size='pes',  type=2 )
    
    pwc = use_dta %>% tukey_hsd( equ_txt )
    
    ## Put together table
    pvl_tbl$DV[iV]          = dep_var_nme[iV]
    
    pvl_tbl$pvalue_one[iV]      = signif(ttt$p[1],digits=3)
    pvl_tbl$pvalue_two[iV]      = signif(ttt$p[2],digits=3)
    pvl_tbl$pvalue_int[iV]      = signif(ttt$p[3],digits=3)
    
    pvl_tbl$statistic_one[iV]   = ttt$F[1]
    pvl_tbl$statistic_two[iV]   = ttt$F[2]
    pvl_tbl$statistic_int[iV]   = ttt$F[3]
    
    pvl_tbl$df_num_one[iV]      = ttt$DFn[1]
    pvl_tbl$df_num_two[iV]      = ttt$DFn[2]
    pvl_tbl$df_num_int[iV]      = ttt$DFn[3]
    
    pvl_tbl$df_den_one[iV]      = ttt$DFd[1]
    pvl_tbl$df_den_two[iV]      = ttt$DFd[2]
    pvl_tbl$df_den_int[iV]      = ttt$DFd[3]
    
    pvl_tbl$report_one[iV]      = paste('F(',signif(ttt$DFn[1],digits=3),';',signif(ttt$DFd[1],digits=3),')=',signif(ttt$F[1],digits=3),'; p=',signif(ttt$p[1],digits=2),sep='')
    pvl_tbl$report_two[iV]      = paste('F(',signif(ttt$DFn[2],digits=3),';',signif(ttt$DFd[2],digits=3),')=',signif(ttt$F[2],digits=3),'; p=',signif(ttt$p[2],digits=2),sep='')
    pvl_tbl$report_int[iV]      = paste('F(',signif(ttt$DFn[3],digits=3),';',signif(ttt$DFd[3],digits=3),')=',signif(ttt$F[3],digits=3),'; p=',signif(ttt$p[3],digits=2),sep='')
    
    pvl_tbl$par_eta_sq_one[iV]  = signif(ttt$pes[1],digits=3)
    pvl_tbl$par_eta_sq_two[iV]  = signif(ttt$pes[2],digits=3)
    pvl_tbl$par_eta_sq_int[iV]  = signif(ttt$pes[3],digits=3)
    
    pvl_tbl$levene[iV]    = signif(lev_tst$p,digits=3)
    pvl_tbl$shapiro[iV]    = signif(shp_tst$p.value,digits=3)

    ## Put together pairwise
    for ( iC in 1:length(par_wse_nme) )
    {
      tuk_org_tbl[[tuk_org_nme[iC]]][iV] = pwc$p.adj[iC]
    }
    
    ## Make Reports
    box_plt <- ggboxplot( use_dta, x = grp_var_one_nme[iG], y = dep_var_nme[iV], color = grp_var_two_nme[iG],
                          ylab = dep_var_nme[iV], xlab = "Groups", add = "dotplot" )
    ggsave(paste(out_put_loc,'/',grp_var_one_nme[iG],'/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)
    
    nor_plt     = ggqqplot(residuals(nor_mdl))
    nor_grp_plt = ggqqplot(use_dta, dep_var_nme[iV], facet.by = grp_var_one_nme[iG])
    out_plt     = ggarrange( nor_plt, nor_grp_plt, ncol = 2, nrow = 1)
    ggsave(paste(out_put_loc,'/',grp_var_one_nme[iG],'/','tests','/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)

  }
 
  out_tbl = cbind( pvl_tbl, tuk_org_tbl) 
  write.csv( out_tbl, paste(out_put_loc,'/',grp_var_one_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
  
}
    
}

