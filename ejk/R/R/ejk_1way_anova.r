ejk_1way_anova <- function( dep_var_loc,
                            grp_loc,
                            out_put_loc,
                            alternative='two.sided') {
  
#dep_var_loc =  '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA/ANOVA/dep_var.mat'
#grp_loc     = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA/ANOVA/grp_var.mat'
#out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/epd_emo/tractFA/ANOVA'  

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
grp_var_nme = grp_var_nme[-1]

plt_dta = merge( dep_var, grp_var, by="sbj.nme" );

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
                        #method=character(length(dep_var_nme)),
                        report=character(length(dep_var_nme)),
                        par_eta_sq=double(length(dep_var_nme)),
                        levene=double(length(dep_var_nme)),
                        shapiro=double(length(dep_var_nme)),
                        stringsAsFactors=FALSE)
  
  # Put together table for ttests
  equ_txt = as.formula( paste( dep_var_nme[1], ' ~ ', grp_var_nme[iG], sep='') )
  pwc = plt_dta %>% tukey_hsd( equ_txt )
  
  par_wse_nme = levels(plt_dta[[grp_var_nme[1]]])
  
  par_wse_num = length(par_wse_nme)
  par_wse_num = (par_wse_num*(par_wse_num-1)) / 2
  
  tuk_org_tbl = data.frame(matrix(vector(), length(dep_var_nme), par_wse_num))
    tuk_org_nme = data.frame(name=character(par_wse_num),stringsAsFactors=FALSE)
  tts_org_tbl = data.frame(matrix(vector(), length(dep_var_nme), par_wse_num))
    tts_org_nme = data.frame(name=character(par_wse_num),stringsAsFactors=FALSE)
  tts_hlm_tbl = data.frame(matrix(vector(), length(dep_var_nme), par_wse_num))
    tts_hlm_nme = data.frame(name=character(par_wse_num),stringsAsFactors=FALSE)
  tts_fdr_tbl = data.frame(matrix(vector(), length(dep_var_nme), par_wse_num))
    tts_fdr_nme = data.frame(name=character(par_wse_num),stringsAsFactors=FALSE)
  for ( iC in 1:par_wse_num )
  {
    tuk_org_nme$name[iC] = paste( pwc$group1[iC], '_', pwc$group2[iC], '_tuk', sep='')
    tts_org_nme$name[iC] = paste( pwc$group1[iC], '_', pwc$group2[iC], '_tts_org', sep='')
    tts_hlm_nme$name[iC] = paste( pwc$group1[iC], '_', pwc$group2[iC], '_tts_hlm', sep='')
    tts_fdr_nme$name[iC] = paste( pwc$group1[iC], '_', pwc$group2[iC], '_tts_fdr', sep='')
  }
  
  colnames(tuk_org_tbl) = tuk_org_nme$name
  colnames(tts_org_tbl) = tts_org_nme$name
  colnames(tts_hlm_tbl) = tts_hlm_nme$name
  colnames(tts_fdr_tbl) = tts_fdr_nme$name
  
  # Iterate through Dependent Variables
  for ( iV in 1:length(dep_var_nme) )
  {
   
    equ_txt = as.formula( paste( dep_var_nme[iV], ' ~ ', grp_var_nme[iG], sep='') )
    
    use_dta = as.data.frame(list(plt_dta[[grp_var_nme[iG]]], plt_dta[[dep_var_nme[iV]]]), col.names=c(grp_var_nme[iG],dep_var_nme[iV]))
    
    ## Leven Test
    lev_tst = use_dta %>% levene_test(equ_txt)
    
    ## Shapiro Test
    nor_mdl = lm(equ_txt, data = use_dta)
    shp_tst = shapiro_test(residuals(nor_mdl))
    
    ## Run 1-way anova
    ttt = use_dta %>% anova_test( equ_txt, effect.size='pes' )
    
    pwc = use_dta %>% tukey_hsd( equ_txt )
    
    pwt_hlm = use_dta %>% pairwise_t_test( equ_txt, p.adjust.method='holm' )
    
    pwt_fdr = use_dta %>% pairwise_t_test( equ_txt, p.adjust.method='fdr' )
    
    ## Put together table
    pvl_tbl$DV[iV]          = dep_var_nme[iV]
    pvl_tbl$pvalue[iV]      = signif(ttt$p,digits=3)
    pvl_tbl$statistic[iV]   = ttt$F
    pvl_tbl$df_num[iV]      = ttt$DFn
    pvl_tbl$df_den[iV]      = ttt$DFd
    #pvl_tbl$method[iV]      = ttt$method
    
    pvl_tbl$report[iV]      = paste('F(',signif(ttt$DFn,digits=3),';',signif(ttt$DFd,digits=3),')=',signif(ttt$F,digits=3),'; p=',signif(ttt$p,digits=2),sep='')
    
    pvl_tbl$par_eta_sq[iV]  = signif(ttt$pes,digits=3)
    
    pvl_tbl$levene[iV]    = signif(lev_tst$p,digits=3)
    pvl_tbl$shapiro[iV]    = signif(shp_tst$p.value,digits=3)

    ## Put together pairwise
    for ( iC in 1:par_wse_num )
    {
      tuk_org_tbl[[tuk_org_nme$name[iC]]][iV] = pwc$p.adj[iC]
      tts_org_tbl[[tts_org_nme$name[iC]]][iV] = pwt_fdr$p[iC]
      tts_hlm_tbl[[tts_hlm_nme$name[iC]]][iV] = pwt_hlm$p.adj[iC]
      tts_fdr_tbl[[tts_fdr_nme$name[iC]]][iV] = pwt_fdr$p.adj[iC]
    }
    
    ## Make Reports
    box_plt <- ggboxplot( use_dta, x = grp_var_nme[iG], y = dep_var_nme[iV], 
                          ylab = dep_var_nme[iV], xlab = "Groups", add = "dotplot" )
    ggsave(paste(out_put_loc,'/',grp_var_nme[iG],'/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)
    
    nor_plt     = ggqqplot(residuals(nor_mdl))
    nor_grp_plt = ggqqplot(use_dta, dep_var_nme[iV], facet.by = grp_var_nme[iG])
    out_plt     = ggarrange( nor_plt, nor_grp_plt, ncol = 2, nrow = 1)
    ggsave(paste(out_put_loc,'/',grp_var_nme[iG],'/','tests','/',dep_var_nme[iV],'.jpg',sep=''), plot=box_plt)

  }
 
  out_tbl = cbind( pvl_tbl, tuk_org_tbl, tts_org_tbl, tts_hlm_tbl, tts_fdr_tbl) 
  write.csv( out_tbl, paste(out_put_loc,'/',grp_var_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
  
}
    
}

