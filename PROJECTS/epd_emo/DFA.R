library(MASS)
library( R.matlab )
library(klaR)
library( ggpubr )

dta_frs = readMat('/home/ekaestner/Downloads/out_dta_frs.mat') 

dta_frs_nme = names( dta_frs$out.dta.frs[,,1] )
dta_frs = as.data.frame( lapply(dta_frs$out.dta.frs, unlist, use.names=FALSE), 
                         fix.empty.names=FALSE, stringsAsFactors=FALSE, col.names=dta_frs_nme )
dta_frs$sbj.grp = factor(dta_frs$sbj.grp)

ttt_dfa = lda( sbj.grp ~ lhs.ist.cng + rhs.ist.cng + lhs.cng.pra + rhs.cng.pra, na.action="na.omit", data=dta_frs, CV=TRUE )
ttt_cta = table(dta_frs$sbj.grp[c(-24,-25,-29,-40,-54,-95)], ttt_dfa$class) # -24,-25,-29,-40,-54,-95 /// -59,-70,-71,-81,-84,-99,-104,-105,-106
diag(prop.table(ttt_cta,1))
sum(diag(prop.table(ttt_cta,1)))/2

ttt_plt = partimat(sbj.grp ~ lhs.ist.cng + rhs.ist.cng + lhs.cng.pra + rhs.cng.pra, data=dta_frs, method='lda')
ggsave('/home/ekaestner/Downloads/Dep_Non.jpg', plot=ttt_plt)

ttt_dta =        c(1,1,1,1,1,2,2,2,2,12,12,12,12,14,14,14,14,14,2,2)
ttt_grp = factor(c(1,1,1,1,1,1,1,1,1,1, 1, 2, 2, 2, 2, 2, 2, 2, 2,2))  

ttt_dfa = lda( ttt_grp ~ ttt_dta, CV=TRUE )
ttt_cta = table(ttt_grp, ttt_dfa$class)
diag(prop.table(ttt_cta,1))
sum(diag(prop.table(ttt_cta,1)))
