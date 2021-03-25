function dta_out = ejk_ComBat( cfg )

tmp_out_dir = [ cfg.out_dir '/' 'ComBat' '/' ];

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

ejk_chk_dir(tmp_out_dir)
save( [ tmp_out_dir '/' 'dta_var.mat' ], 'dta_var' )
save( [ tmp_out_dir '/' 'btc_var.mat' ], 'btc_var' )
if isfield( cfg, 'cov'); save( [ tmp_out_dir '/' 'cov_var.mat' ], 'cov_var' ); end

%% Initial Figure
btc_nme = unique( cfg.btc );
btc_num = numel(btc_nme);
col_hld = distinguishable_colors(btc_num);
avg_tfa_dta = mean(cfg.dta,2);

bin_wdt = (cfg.ylm(2) - cfg.ylm(1)) / 100;
figure('Visible','off'); xlim( cfg.ylm ); hold on;
for iPL = 1:btc_num
    fig_hld = histogram( avg_tfa_dta(strcmpi( cfg.btc, btc_nme{iPL}),1) );
    fig_hld.Normalization='probability';
    fig_hld.BinWidth = bin_wdt;
    fig_hld.FaceColor = col_hld(iPL,:);
    fig_hld.EdgeColor = 'none';
    fig_hld.FaceAlpha = 0.50;
end
print( gcf, [ tmp_out_dir '/' 'OriginalData' '_' cfg.btc_nme{1} '.png'], '-dpng' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_ComBat.R'')\n',sve_cmd);

if ~isfield( cfg, 'cov')
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', '''', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir, tmp_out_dir);
else
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', ''%s/cov_var.mat'', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir, tmp_out_dir, tmp_out_dir);
end

cell2csv( [tmp_out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' tmp_out_dir '/example_R_script.r' ] );

%% Load and return
dta_out = mmil_readtext( [ tmp_out_dir '/' 'out_dta.csv' ] );
dta_out = cell2mat( dta_out(2:end, 2:end) )';

%% Second Figure
btc_nme = unique( cfg.btc );
btc_num = numel(btc_nme);
col_hld = distinguishable_colors(btc_num);
avg_tfa_dta = mean(dta_out,2);

bin_wdt = (cfg.ylm(2) - cfg.ylm(1)) / 100;
figure('Visible','off'); xlim( cfg.ylm ); hold on;
for iPL = 1:btc_num
    fig_hld = histogram( avg_tfa_dta(strcmpi( cfg.btc, btc_nme{iPL}),1) );
    fig_hld.Normalization='probability';
    fig_hld.BinWidth = bin_wdt;
    fig_hld.FaceColor = col_hld(iPL,:);
    fig_hld.EdgeColor = 'none';
    fig_hld.FaceAlpha = 0.50;
end
print( gcf, [ tmp_out_dir '/' 'ComBatData' '_' cfg.btc_nme{1} '.png'], '-dpng' )

end