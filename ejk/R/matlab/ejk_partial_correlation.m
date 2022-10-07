function ejk_partial_correlation( cfg )

%%
dta_one.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.lbl_one)
    try dta_one.(cfg.lbl_one{iD}) = cfg.dta_one(:,iD); catch fprintf('Failure Field Name: %s \n',cfg.lbl_one{iD}); end
end

dta_two.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.lbl_two)
    try dta_two.(cfg.lbl_two{iG}) = cfg.dta_two(:,iG); catch fprintf('Failure Field Name: %s \n',cfg.lbl_two{iG}); end
end

dta_thr.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.lbl_cov)
    try dta_thr.(cfg.lbl_cov{iG}) = cfg.dta_cov(:,iG); catch fprintf('Failure Field Name: %s \n',cfg.lbl_cov{iG}); end
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_one.mat' ], 'dta_one' )
save( [ cfg.out_dir '/' 'dta_two.mat' ], 'dta_two' )
save( [ cfg.out_dir '/' 'dta_thr.mat' ], 'dta_thr' )

%%
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_partial_correlation.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_partial_correlation( ''%s/dta_one.mat'', ''%s/dta_two.mat'', ''%s/dta_thr.mat'', ''%s'', ''%s'')', sve_cmd, ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), cfg.cor_typ );
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' ejk_fix_path(cfg.out_dir) '/example_R_script.r' ] );

%%

end