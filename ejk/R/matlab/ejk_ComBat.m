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
ejk_chk_dir([ tmp_out_dir '/' 'Original' '/' ])

btc_nme = unique( cfg.btc );
btc_num = numel(btc_nme);
col_hld = distinguishable_colors(btc_num);

for iP = 1:size(cfg.dta,2)
    
    figure('Visible','off'); hold on; 
    if isfield(cfg,'ylm'); xlim( cfg.ylm ); bin_wdt = (cfg.ylm(2) - cfg.ylm(1)) / 50;
    else ylm(1) = min(cfg.dta(:,iP));
    	 ylm(2) = max(cfg.dta(:,iP)); 
         ylm(1) = ylm(1) - range(ylm)*.1;
         ylm(2) = ylm(2) + range(ylm)*.1;
    	 xlim( ylm ); 
         bin_wdt = (ylm(2) - ylm(1)) / 50; end
    
    for iPL = 1:btc_num
        fig_hld = histogram( cfg.dta(strcmpi( cfg.btc, btc_nme{iPL}),iP) );
        fig_hld.Normalization='probability';
        fig_hld.BinWidth = bin_wdt;
        fig_hld.FaceColor = col_hld(iPL,:);
        fig_hld.EdgeColor = 'none';
        fig_hld.FaceAlpha = 0.50;
    end
    legend(btc_nme)
    print( gcf, [ tmp_out_dir '/' 'Original' '/' 'OriginalData' '_' cfg.btc_nme{1} '_' cfg.dta_nme{iP} '.png'], '-dpng' )
    close all
    
end
    
%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_ComBat.r'')\n',sve_cmd);

if ~isfield( cfg, 'cov')
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', '''', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir, tmp_out_dir);
else
    sve_cmd = sprintf('%sejk_ComBat( ''%s/dta_var.mat'', ''%s/btc_var.mat'', ''%s/cov_var.mat'', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir, tmp_out_dir, tmp_out_dir);
end

cell2csv( [tmp_out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' tmp_out_dir '/example_R_script.r' ] );

%% Load and return
dta_out = mmil_readtext( [ tmp_out_dir '/' 'out_dta.csv' ] );
dta_out(cellfun(@(x) strcmpi(x,'NA'),dta_out)) = {NaN};
dta_out = cell2mat( dta_out(2:end, 2:end) )';

%% Second Figure
if cfg.plt
    ejk_chk_dir([ tmp_out_dir '/' 'ComBat' '/' ])
    
    btc_nme = unique( cfg.btc );
    btc_num = numel(btc_nme);
    col_hld = distinguishable_colors(btc_num);
    
    for iP = 1:size(dta_out,2)       
        figure('Visible','off'); hold on;
        if isfield(cfg,'ylm'); xlim( cfg.ylm ); bin_wdt = (cfg.ylm(2) - cfg.ylm(1)) / 50;
        else ylm(1) = min(dta_out(:,iP));
            ylm(2) = max(dta_out(:,iP));
            ylm(1) = ylm(1) - range(ylm)*.1;
            ylm(2) = ylm(2) + range(ylm)*.1;
             xlim( ylm ); 
             bin_wdt = (ylm(2) - ylm(1)) / 50; end
         
        for iPL = 1:btc_num
            fig_hld = histogram( dta_out(strcmpi( cfg.btc, btc_nme{iPL}),iP) );
            fig_hld.Normalization='probability';
            fig_hld.BinWidth = bin_wdt;
            fig_hld.FaceColor = col_hld(iPL,:);
            fig_hld.EdgeColor = 'none';
            fig_hld.FaceAlpha = 0.50;
        end
        legend(btc_nme)
        print( gcf, [ tmp_out_dir '/' 'ComBat'  '/' 'ComBatData' '_' cfg.btc_nme{1} '_' cfg.dta_nme{iP} '.png'], '-dpng' )
        close all
    end
    
end

end