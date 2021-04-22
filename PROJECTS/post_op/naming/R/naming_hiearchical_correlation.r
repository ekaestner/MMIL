library( R.matlab )
library(robustbase)

dep_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SpecificCor/HiearchicalRegression/dep_var.mat'
prd_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SpecificCor/HiearchicalRegression/prd_var.mat'
prd_nme_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SpecificCor/HiearchicalRegression/prd_nme.mat'

blk_nme = c( 'blk.1', 'blk.2', 'blk.3', 'blk.4', 'blk.5', 'blk.6', 'blk.7', 'blk.8', 'blk.9' )

#############

dep_var = readMat(dep_var_loc) 
dep_var_nme = names( dep_var$dep.var[,,1] )
dep_var = as.data.frame( lapply(dep_var$dep.var, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dep_var_nme )
dep_var_nme = dep_var_nme[-1]

prd_var = readMat(prd_var_loc) 
prd_var_nme = names( prd_var$prd.var[,,1] )
prd_var = as.data.frame( lapply(prd_var$prd.var, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=prd_var_nme )

prd_nme = readMat(prd_nme_loc) 
prd_nme_nme = names( prd_nme$prd.nme[,,1] )
prd_nme = as.data.frame( lapply(prd_nme$prd.nme, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=prd_nme_nme )

use_dta = merge( dep_var, prd_var, by="sbj.nme" )

# Check - BNT ###############################
cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xrh.fusiform, method=c('pearson') )
cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xrh.fusiform, method=c('spearman') )
summary(lm(xbnt.raw.scr.pst ~ xrh.fusiform, data=use_dta))
summary(lmrob(xbnt.raw.scr.pst ~ xrh.fusiform, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xL.ILF, method=c('pearson') )
cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xL.ILF, method=c('spearman') )
summary(lm(xbnt.raw.scr.pst ~ xL.ILF, data=use_dta))
summary(lmrob(xbnt.raw.scr.pst ~ xL.ILF, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xLeft.Hippocampus, method=c('pearson') )
cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xLeft.Hippocampus, method=c('spearman') )
summary(lm(xbnt.raw.scr.pst ~ xLeft.Hippocampus, data=use_dta))
summary(lmrob(xbnt.raw.scr.pst ~ xLeft.Hippocampus, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xbnt.raw.scr, method=c('pearson') )
cor.test( use_dta$xbnt.raw.scr.pst, use_dta$xbnt.raw.scr, method=c('spearman') )
summary(lmrob(xbnt.raw.scr.pst ~ xbnt.raw.scr, data=use_dta, method = "MM", setting="KS2014"))

# Hiearchical Regression - BNT ###############################
blk_one = lm( xbnt.raw.scr.pst ~  xbnt.raw.scr + xLeft.Hippocampus, use_dta)
summary(blk_one)

blk_two = lm( xbnt.raw.scr.pst ~  xbnt.raw.scr + xLeft.Hippocampus + xrh.fusiform + xL.ILF, use_dta)
summary(blk_two)

anova( blk_one, blk_two)

# Robust Regression BNT ###############################
blk_one = lmrob(xbnt.raw.scr.pst ~ xbnt.raw.scr + xLeft.Hippocampus, data=use_dta, method = "MM", setting="KS2014")
summary(blk_one)

blk_two = lmrob(xbnt.raw.scr.pst ~ xbnt.raw.scr + xLeft.Hippocampus + xrh.fusiform + xL.ILF, data=use_dta, method = "MM", setting="KS2014")
summary(blk_two)

# Check - ANT ###############################
cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xrh.fusiform, method=c('pearson') )
cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xrh.fusiform, method=c('spearman') )
summary(lmrob(xant.mem.raw.scr.pst ~ xrh.fusiform, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xL.ILF, method=c('pearson') )
cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xL.ILF, method=c('spearman') )
summary(lmrob(xant.mem.raw.scr.pst ~ xL.ILF, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xLeft.Hippocampus, method=c('pearson') )
cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xLeft.Hippocampus, method=c('spearman') )
summary(lmrob(xant.mem.raw.scr.pst ~ xLeft.Hippocampus, data=use_dta, method = "MM", setting="KS2014"))

cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xant.mem.raw.scr, method=c('pearson') )
cor.test( use_dta$xant.mem.raw.scr.pst, use_dta$xant.mem.raw.scr, method=c('spearman') )
summary(lmrob(xant.mem.raw.scr.pst ~ xant.mem.raw.scr, data=use_dta, method = "MM", setting="KS2014"))

# Hierchical Regression ANT ###############################
blk_one = lm( xant.mem.raw.scr.pst ~  xant.mem.raw.scr + xLeft.Hippocampus, use_dta)
summary(blk_one)

blk_two = lm( xant.mem.raw.scr.pst ~ xant.mem.raw.scr + xLeft.Hippocampus + xL.ILF + xrh.fusiform, use_dta)
summary(blk_two)

blk_one_neu = lm( xant.mem.raw.scr.pst ~ xL.ILF + xrh.fusiform, use_dta)
summary(blk_one_neu)

anova( blk_one, blk_two)

# Robust Regression ANT ###############################
blk_one = lmrob(xant.mem.raw.scr.pst ~ xant.mem.raw.scr + xLeft.Hippocampus, data=use_dta, method = "MM", setting="KS2014")
summary(blk_one)

blk_one_neu = lmrob(xant.mem.raw.scr.pst ~ xrh.fusiform + xL.ILF, data=use_dta, method = "MM", setting="KS2014")
summary(blk_one_neu)

blk_two = lmrob(xant.mem.raw.scr.pst ~ xant.mem.raw.scr + xLeft.Hippocampus + xrh.fusiform + xL.ILF, data=use_dta, method = "MM", setting="KS2014")
summary(blk_two)
