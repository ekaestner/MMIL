library( R.matlab )
library(robustbase)
library(ppcor)

####### LOAD #######
wmp_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_left/dta_one.mat'
fib_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_left/dta_one.mat'
cvr_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/Correlation/LTLE_post_cln/dta_one.mat'
prd_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_left/dta_two.mat'

wmp_var = readMat(wmp_var_loc) 
wmp_var_nme = names( wmp_var$dta.one[,,1] )
wmp_var = as.data.frame( lapply(wmp_var$dta.one, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=wmp_var_nme )
wmp_var_nme = wmp_var_nme[-1]

fib_var = readMat(fib_var_loc) 
fib_var_nme = names( fib_var$dta.one[,,1] )
fib_var = as.data.frame( lapply(fib_var$dta.one, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=fib_var_nme )
fib_var_nme = fib_var_nme[-1]

cvr_var = readMat(cvr_var_loc) 
cvr_var_nme = names( cvr_var$dta.one[,,1] )
cvr_var = as.data.frame( lapply(cvr_var$dta.one, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=cvr_var_nme )
cvr_var_nme = cvr_var_nme[-1]

prd_var = readMat(prd_var_loc) 
prd_var_nme = names( prd_var$dta.two[,,1] )
prd_var = as.data.frame( lapply(prd_var$dta.two, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=prd_var_nme )
prd_var_nme = prd_var_nme[-1]

####### BNT / R-FUSIFORM #######
dta_hld = cbind( prd_var$bnt.raw.scr.pst, wmp_var$xrh.fusiform, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )

####### ANT / R-FUSIFORM #######
dta_hld = cbind( prd_var$ant.mem.raw.scr.pst, wmp_var$xrh.fusiform, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )

####### BNT / L-ILF #######
dta_hld = cbind( prd_var$bnt.raw.scr.pst, fib_var$xL.ILF, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )

####### ANT / L-ILF #######
dta_hld = cbind( prd_var$ant.mem.raw.scr.pst, fib_var$xL.ILF, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )

####### BNT / L-IFOF #######
dta_hld = cbind( prd_var$bnt.raw.scr.pst, fib_var$xL.IFO, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )

####### ANT / L-IFOF #######
dta_hld = cbind( prd_var$ant.mem.raw.scr.pst, fib_var$xL.IFO, cvr_var$Educ)
dta_hld = as.data.frame(na.omit(dta_hld))

cor.test( dta_hld[,1], dta_hld[,2], method=c('spearman') )
cor.test( dta_hld[,1], dta_hld[,2], method=c('pearson') )
pcor.test( dta_hld[,1], dta_hld[,2], dta_hld[,3], method=c('spearman') )










