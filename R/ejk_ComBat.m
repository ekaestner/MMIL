function dta_out = ejk_ComBat( cfg )


%% 
dta_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    dta_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD);
end

btc_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.btc_nme)
    btc_var.(cfg.btc_nme{iG}) = cfg.btc(:,iG);
end

if isfield( cfg, 'cov')
    cov_var.sbj_nme   = cfg.sbj_nme;
    for iG = 1:numel(cfg.cov_nme)
        cov_var.(cfg.cov_nme{iG}) = cfg.cov(:,iG);
    end
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_var.mat' ], 'dta_var' )
save( [ cfg.out_dir '/' 'btc_var.mat' ], 'btc_var' )
if isfield( cfg, 'cov'); save( [ cfg.out_dir '/' 'cov_var.mat' ], 'cov_var' ); end

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_ComBat.R'')\n',sve_cmd);

if ~isfield( cfg, 'cov')
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', '''', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir);
else
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', ''%s/cov_var.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir, cfg.out_dir);
end

cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' cfg.out_dir '/example_R_script.r' ] );

%% Load and return
dta_out = mmil_readtext( [ cfg.out_dir '/' 'out_dta.csv' ] );
dta_out = cell2mat( dta_out(2:end, 2:end) )';

end