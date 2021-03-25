ejk_pearson_correlation <- function( dta_one_loc,
                                     dta_two_loc,
                                     out_put_loc ) {
  
  #dta_one_loc = '/home/ekaestner/Downloads/pearson/dta_one.mat'
  #dta_two_loc = '/home/ekaestner/Downloads/pearson/dta_two.mat'
  #out_put_loc = '/home/ekaestner/Downloads/pearson/'  
  
  library( R.matlab )
  library( rstatix )
  library( ggpubr )
  library( dplyr )
  
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
  
  plt_dta = merge( dta_one, dta_two, by="sbj.nme" );
  
  dir.create( paste(out_put_loc,sep=''), showWarnings = FALSE)
  
  for ( iG in 1:length(dta_two_nme) )
  {
    
    dir.create( paste(out_put_loc,'/',dta_two_nme[iG],sep=''), showWarnings = FALSE)
    dir.create( paste(out_put_loc,'/',dta_two_nme[iG],'/','tests','/',sep=''), showWarnings = FALSE)
    
    # Put together table
    pvl_tbl = data.frame( DV=character(length(dta_one_nme)),
                          statistic=double(length(dta_one_nme)),
                          pvalue=double(length(dta_one_nme)),
                          df_num=double(length(dta_one_nme)),
                          report=character(length(dta_one_nme)),
                          r_sq=double(length(dta_one_nme)),
                          shapiro1=double(length(dta_one_nme)),
                          shapiro2=double(length(dta_one_nme)),
                          stringsAsFactors=FALSE)
    
    # Iterate through Dependent Variables
    for ( iV in 1:length(dta_one_nme) )
    {
      
      use_dta = as.data.frame(list(plt_dta[[dta_two_nme[iG]]], plt_dta[[dta_one_nme[iV]]]), col.names=c(dta_two_nme[iG],dta_one_nme[iV]))
      
      ## Run tests
      one_shp = shapiro.test( use_dta[[dta_one_nme[iV]]] )
      two_shp = shapiro.test( use_dta[[dta_two_nme[iG]]] )
      
      ## Run Pearson Correlation
      ttt = cor.test(use_dta[[dta_two_nme[iG]]], use_dta[[dta_one_nme[iV]]], method = "pearson")
      
      ## Put together table
      pvl_tbl$DV[iV]          = dta_one_nme[iV]
      pvl_tbl$pvalue[iV]      = signif(ttt$p.value,digits=3)
      pvl_tbl$statistic[iV]   = ttt$estimate
      pvl_tbl$df_num[iV]      = ttt$parameter
      pvl_tbl$method[iV]      = ttt$method
      pvl_tbl$report[iV]      = paste('r=',signif(ttt$estimate,digits=2),', p=',signif(ttt$p.value,digits=2),sep='')
      pvl_tbl$r_sq[iV]        = ttt$estimate^2
      pvl_tbl$shapiro1[iV]    = one_shp$p.value
      pvl_tbl$shapiro2[iV]    = two_shp$p.value
      
      ## Make Reports
      sct_plt = ggscatter(use_dta, x = dta_two_nme[iG], y = dta_one_nme[iV], 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = dta_two_nme[iG], ylab = dta_one_nme[iV])
      ggsave(paste(out_put_loc,'/',dta_two_nme[iG],'/',dta_one_nme[iV],'.jpg',sep=''), plot=sct_plt)
      
      dta_one_plt = ggqqplot(use_dta[[dta_one_nme[iV]]], ylab = dta_one_nme[iG])
      dta_two_plt = ggqqplot(use_dta[[dta_two_nme[iG]]], ylab = dta_two_nme[iG] )
      out_plt     = ggarrange( dta_one_plt, dta_two_plt, ncol = 2, nrow = 1)
      ggsave(paste(out_put_loc,'/',dta_two_nme[iG],'/','tests','/',dta_one_nme[iV],'.jpg',sep=''), plot=out_plt)
      
    }

    write.csv( pvl_tbl, paste(out_put_loc,'/',dta_two_nme[iG],'/','output_table.csv',sep=''), row.names=FALSE)
    
  }
  
}

