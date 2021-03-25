function dta_out = ejk_residual_matrix( cfg )

%% 
dta_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    dta_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD);
end

cov_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.cov_nme)
    cov_var.(cfg.cov_nme{iG}) = cfg.cov(:,iG);
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_var.mat' ], 'dta_var' )
save( [ cfg.out_dir '/' 'cov_var.mat' ], 'cov_var' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_residual_matrix.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_residual_matrix( ''%s/dta_var.mat'', ''%s/cov_var.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' cfg.out_dir '/example_R_script.r' ] );

%% Load and return
dta_out = mmil_readtext( [ cfg.out_dir '/' 'residual_matrix.csv' ] );
dta_out = cell2mat( dta_out(2:end, 2:end) );

end