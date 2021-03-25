
ttt_add = as.data.frame(sample(32), row.names=row.names(mtcars) )
colnames(ttt_add) = c('rando')

ttt_new = merge( mtcars, ttt_add, by=0 )
rownames(ttt_new) = ttt_new$Row.names

Model = lm( mpg ~ rando, data=ttt_new)

resid_hld = resid(Model)
resid_dta = as.data.frame(resid_hld)
colnames(resid_dta) = c('resid_hld')

ttt_new = merge( ttt_new, resid_dta, by=0 )
ttt_new = ttt_new[-1]

dta_one_plt = ggscatter( ttt_new, 'resid_hld', 'mpg'  )
dta_two_plt = ggscatter( ttt_new, 'resid_hld', 'rando'  )
dta_thr_plt = ggscatter( ttt_new, 'mpg', 'rando'  )
ggarrange( dta_one_plt, dta_two_plt, dta_thr_plt, ncol = 2, nrow = 2)

cbind( ttt_new$mpg, ttt_new$resid_hld, ttt_add$rando )