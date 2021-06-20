ejk_fisher_test <- function( dta_one_loc,
                             dta_two_loc,
                             out_put_loc ) {
  
  # https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
  
  #dta_one_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/dta_one.mat'
  #dta_two_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/dta_two.mat'
  #out_put_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/'

  library( 'R.matlab' )  
  library('corrplot')
  library('Hmisc')
  library('dplyr')
  library('pracma')
  library('gridExtra')
  
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
  
  if (!(sum(!is.na(match( dta_one_nme, dta_two_nme)))==length(dta_one_nme))){
    dta_use = merge( dta_one, dta_two, by="sbj.nme" )
  } else {
    dta_use = dta_one
  }
  dta_use = dta_use[,-1]
  
  pvl_tbl = data.frame( DV=character(length(dta_one_nme)),
                        statistic=double(length(dta_one_nme)),
                        pvalue=double(length(dta_one_nme)),
                        method=character(length(dta_one_nme)),
                        alternative=character(length(dta_one_nme)), 
                        report=character(length(dta_one_nme)),
                        stringsAsFactors=FALSE)
  
  for (iC in 1:length(dta_one_nme)){
    
    ttt = fisher.test( table( dta_use[[ dta_one_nme[iC]]], dta_use[[ dta_two_nme[iC]]] ) )
    
    pvl_tbl$DV[iC]          = dta_one_nme[iC]
    pvl_tbl$pvalue[iC]      = ttt$p.value
    if (!is.null(ttt$estimate)){ pvl_tbl$statistic[iC]   = ttt$estimate }
    pvl_tbl$alternative[iC] = ttt$alternative
    pvl_tbl$method[iC]      = ttt$method
    if (!is.null(ttt$estimate)){ pvl_tbl$report[iC]      = paste('FET=',signif(ttt$estimate,digits=3),'; p=',signif(ttt$p.value,digits=2),sep='')
    } else {pvl_tbl$report[iC]      = paste('FET; p=',signif(ttt$p.value,digits=2),sep='')}
    
    png( paste(out_put_loc,'/',dta_one_nme[iC],'.png',sep='') )
    p<-tableGrob( table( dta_use[[ dta_one_nme[iC]]], dta_use[[ dta_two_nme[iC]]] ) )
    grid.arrange(p)
    dev.off()
    
  }
  
  write.csv( pvl_tbl, paste(out_put_loc,'/','output_table.csv',sep=''), row.names=FALSE)
  
  }