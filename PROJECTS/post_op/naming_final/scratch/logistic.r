library( R.matlab )
library(pscl)
library(ROCR)
library(lmtest)
library(rcompanion)

####### LOAD #######
wmp_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/wmparc_FA_wm/Raw/tle_post_3T_ATLonly_left/dta_one.mat'
fib_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_left/dta_one.mat'
cvr_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/Correlation/LTLE_post_cln/dta_one.mat'
cog_var_loc = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Cognitive/Correlation/LTLE_pre_post/dta_one.mat'
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

cog_var = readMat(cog_var_loc) 
cog_var_nme = names( cog_var$dta.one[,,1] )
cog_var = as.data.frame( lapply(cog_var$dta.one, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=cog_var_nme )
cog_var_nme = cog_var_nme[-1]

prd_var = readMat(prd_var_loc) 
prd_var_nme = names( prd_var$dta.two[,,1] )
prd_var = as.data.frame( lapply(prd_var$dta.two, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=prd_var_nme )
prd_var_nme = prd_var_nme[-1]

### BNT - Clinical #######################################################
dta_hld = cbind( prd_var$bnt.raw.scr.pst, cvr_var$Educ, wmp_var$xrh.fusiform, fib_var$xL.ILF, cog_var$bnt.raw.scr)
dta_hld = as.data.frame(na.omit(dta_hld))
dta_hld = cbind( dta_hld, ifelse(dta_hld[,1]<=-1.5,'Impaired','Non-Impaired') )
names(dta_hld)[6]='V6'

mdl_hld_cln = glm( V6 ~ V2 + V5, family=binomial(link='logit'), data=dta_hld )

#summary(mdl_hld_cln)
#anova(mdl_hld_cln, test="Chisq")
pR2(mdl_hld_cln)

#fit_rsl = predict( mdl_hld_cln, dta_hld, type='response')
#fit_rsl = ifelse(fit_rsl<.5, 'Impaired','Non-Impaired')
#err_hld = mean( fit_rsl != dta_hld[,5])
#1-err_hld

fit_auc     = predict( mdl_hld_cln, dta_hld, type='response')
fit_auc_prd = prediction( fit_auc, dta_hld$V6 )
#fit_auc_prf = performance(fit_auc_prd, measure='tpr', x.measure='fpr')
#plot(fit_auc_prf)

auc_out = performance( fit_auc_prd, measure='auc' )
auc_out@y.values[[1]]

### BNT - WHITE #######################################################
dta_hld = cbind( prd_var$bnt.raw.scr.pst, cvr_var$Educ, wmp_var$xrh.fusiform, fib_var$xL.ILF, cog_var$bnt.raw.scr)
dta_hld = as.data.frame(na.omit(dta_hld))
dta_hld = cbind( dta_hld, ifelse(dta_hld[,1]<=-1.5,'Impaired','Non-Impaired') )
names(dta_hld)[6]='V6'

mdl_hld = glm( V6 ~ V3 + V4, family=binomial(link='logit'), data=dta_hld )

#summary(mdl_hld)
anova(mdl_hld, test="Chisq")
lrtest(mdl_hld)
anova(mdl_hld,update(mdl_hld, ~1), test='Chisq')
pR2(mdl_hld)

#fit_rsl = predict( mdl_hld, dta_hld, type='response')
#fit_rsl = ifelse(fit_rsl<.5, 'Impaired','Non-Impaired')
#err_hld = mean( fit_rsl != dta_hld[,5])
#1-err_hld

fit_auc     = predict( mdl_hld, dta_hld, type='response')
fit_auc_prd = prediction( fit_auc, dta_hld$V6 )
#fit_auc_prf = performance(fit_auc_prd, measure='tpr', x.measure='fpr')
#plot(fit_auc_prf)

auc_out = performance( fit_auc_prd, measure='auc' )
auc_out@y.values[[1]]

### COMPARE #######################################################
anova(mdl_hld_cln,mdl_hld,test ="Chisq")

lrtest(mdl_hld_cln, mdl_hld)

### ANT - Clinical #######################################################
dta_hld = cbind( prd_var$ant.mem.raw.scr.pst, cvr_var$Educ, fib_var$xL.IFO, fib_var$xL.ILF, cog_var$ant.mem.raw.scr)
dta_hld = as.data.frame(na.omit(dta_hld))
dta_hld = cbind( dta_hld, ifelse(dta_hld[,1]<=-1.5,'Impaired','Non-Impaired') )
names(dta_hld)[6]='V6'

mdl_hld_cln = glm( V6 ~ V2 + V5, family=binomial(link='logit'), data=dta_hld )

#summary(mdl_hld_cln)
anova(mdl_hld_cln, test="Chisq")
pR2(mdl_hld_cln)

#fit_rsl = predict( mdl_hld_cln, dta_hld, type='response')
#fit_rsl = ifelse(fit_rsl<.5, 'Impaired','Non-Impaired')
#err_hld = mean( fit_rsl != dta_hld[,5])
#1-err_hld

fit_auc     = predict( mdl_hld_cln, dta_hld, type='response')
fit_auc_prd = prediction( fit_auc, dta_hld$V6 )
#fit_auc_prf = performance(fit_auc_prd, measure='tpr', x.measure='fpr')
#plot(fit_auc_prf)

auc_out = performance( fit_auc_prd, measure='auc' )
auc_out@y.values[[1]]

### ANT #######################################################
dta_hld = cbind( prd_var$ant.mem.raw.scr.pst, cvr_var$Educ, fib_var$xL.IFO, fib_var$xL.ILF, cog_var$ant.mem.raw.scr)
dta_hld = as.data.frame(na.omit(dta_hld))
dta_hld = cbind( dta_hld, ifelse(dta_hld[,1]<=-1.5,'Impaired','Non-Impaired') )
names(dta_hld)[6]='V6'

mdl_hld = glm( V6 ~ V3 + V4, family=binomial(link='logit'), data=dta_hld )

#summary(mdl_hld)
anova(mdl_hld, test="Chisq")
lrtest(mdl_hld)
anova(mdl_hld,update(mdl_hld, ~1), test='Chisq')
pR2(mdl_hld)

#fit_rsl = predict( mdl_hld, dta_hld, type='response')
#fit_rsl = ifelse(fit_rsl<.5, 'Impaired','Non-Impaired')
#err_hld = mean( fit_rsl != dta_hld[,5])
#1-err_hld

fit_auc     = predict( mdl_hld, dta_hld, type='response')
fit_auc_prd = prediction( fit_auc, dta_hld$V6 )
#fit_auc_prf = performance(fit_auc_prd, measure='tpr', x.measure='fpr')
#plot(fit_auc_prf)

auc_out = performance( fit_auc_prd, measure='auc' )
auc_out@y.values[[1]]

### COMPARE #######################################################
anova(mdl_hld_cln,mdl_hld,test ="Chisq")

lrtest(mdl_hld_cln, mdl_hld)





