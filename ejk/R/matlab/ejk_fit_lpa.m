function ejk_fit_lpa( cfg )

tmp_out_dir = [ cfg.out_dir '/'  ]; % 'LPA' '/'

%% Save Data
dep_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    dep_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD);
end

ejk_chk_dir(tmp_out_dir)
save( [ tmp_out_dir '/' 'dep_var.mat' ], 'dep_var' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_fit_lpa.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_fit_lpa( ''%s/dep_var.mat'', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir);
cell2csv( [tmp_out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' tmp_out_dir '/example_R_script.r' ] );

%% Make plot
dta_plt = mmil_readtext( [ tmp_out_dir '/' 'model_fit_measures.csv' ] );
dta_plt_nme = dta_plt(1,2:end);
dta_plt = cell2mat( dta_plt(2:end, 2:end) );

num_plt = ceil( sqrt(numel(dta_plt_nme)) );
figure('Visible','off','Position',[0 0 1080 1080])
for iP = 1:numel(dta_plt_nme)
    subplot( num_plt, num_plt, iP)
    plot( 1:size( dta_plt, 1), [nan abs(diff(dta_plt(:,iP)))'], 'k', 'Marker', 'o' )
    xlim([0 size( dta_plt, 1)])
    title(dta_plt_nme{iP})    
end
tightfig();
print( gcf, [ tmp_out_dir '/' 'Model_Fit_Diff.png'], '-dpng' );
close all