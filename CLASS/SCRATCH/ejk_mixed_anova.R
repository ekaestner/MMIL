ejk_mixed_anova = function( dta_tbl, sbj_nme_name, within_name, between_name ){

  ## GATHER DATA ################################################
  col_nme = colnames(dta_tbl)
  
  wth_col_nme = setdiff( col_nme, c( sbj_nme_name, between_name ))

  dta_use = gather( dta_tbl, key = within, value="score", 
                    all_of(wth_col_nme))

  ## Plot Data ################################################  
  bxp <- ggboxplot(
    dta_use, x = "within", y = "score",
    color = between_name, palette = "jco"
  )
  bxp
  
  
}