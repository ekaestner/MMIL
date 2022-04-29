
function ejk_cross_cor( cfg )

if ~isfield(cfg,'force_plot'); cfg.force_plot = 0; end

%%
dta_one.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.lbl_one)
    try dta_one.(cfg.lbl_one{iD}) = cfg.dta_one(:,iD); catch fprintf('Failure Field Name: %s \n',cfg.lbl_one{iD}); end
end

dta_two.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.lbl_two)
    try dta_two.(cfg.lbl_two{iG}) = cfg.dta_two(:,iG); catch fprintf('Failure Field Name: %s \n',cfg.lbl_two{iG}); end
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_one.mat' ], 'dta_one' )
save( [ cfg.out_dir '/' 'dta_two.mat' ], 'dta_two' )

%%
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_cross_cor.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_cross_cor( ''%s/dta_one.mat'', ''%s/dta_two.mat'', ''%s'', ''%s'', %f, %f)', sve_cmd, ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), ejk_fix_path(cfg.out_dir), cfg.cor_typ, cfg.pvl_cut, cfg.pvl_lib);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' ejk_fix_path(cfg.out_dir) '/example_R_script.r' ] );

%% Make Scatter plots of significant correlations
ejk_chk_dir([  cfg.out_dir '/' 'plots' '/'])

if exist([ cfg.out_dir '/' 'cross_correlation_significances_liberal.csv' ]) && ~cfg.force_plot
    
    plt_dta = mmil_readtext([ cfg.out_dir '/' 'cross_correlation_significances_liberal.csv' ]);
    plt_dta = plt_dta(2:end, :);
    
    for iP = 1:size( plt_dta, 1)
        
        fcfg = [];
        
        fcfg.xdt     = { dta_two.(mmil_spec_char(plt_dta{iP,3},{'.'})) };
        fcfg.ydt     = { dta_one.(mmil_spec_char(plt_dta{iP,2},{'.'})) };
        
        fcfg.trd_lne = [ 1 0 ];
        
        fcfg.fce_col = { rgb('blue')  };
        fcfg.edg_col = { rgb('black') };
        
        fcfg.xlb = { mmil_spec_char(plt_dta{iP,3},{'.'}) };
        fcfg.ylb = { mmil_spec_char(plt_dta{iP,2},{'.'}) };
        
        fcfg.ttl = ['r = ' num2str(plt_dta{iP,4}) '  p = ' num2str(roundsd(plt_dta{iP,5},2))];
        
        fcfg.out_dir = [  cfg.out_dir '/' 'plots' '/'];
        fcfg.out_nme = [ mmil_spec_char(plt_dta{iP,3},{'.'}) '__BY__' mmil_spec_char(plt_dta{iP,2},{'.'}) ];
        
        ejk_scatter(fcfg)
        
    end
    
elseif cfg.force_plot
    
    rvl_dta = mmil_readtext([ cfg.out_dir '/' 'cross_correlation_rvalues.csv' ]);
    pvl_dta = mmil_readtext([ cfg.out_dir '/' 'cross_correlation_pvalues.csv' ]);
    
    for iR = 2:size( rvl_dta, 1)
        for iC = 2:size( rvl_dta, 2)
            try
            fcfg = [];
            
            fcfg.xdt     = { dta_two.(mmil_spec_char(rvl_dta{1,iC},{'.'})) };
            fcfg.ydt     = { dta_one.(mmil_spec_char(rvl_dta{iR,1},{'.'})) };
            
            fcfg.trd_lne = [ 1 0 ];
            
            fcfg.fce_col = { rgb('blue')  };
            fcfg.edg_col = { rgb('black') };
            
            fcfg.xlb = { mmil_spec_char(rvl_dta{1,iC},{'.'}) };
            fcfg.ylb = { mmil_spec_char(rvl_dta{iR,1},{'.'}) };
            
            if ~ischar(class(pvl_dta{iR,iC}))
            fcfg.ttl = ['r = ' num2str(rvl_dta{iR,iC}) '  p = ' num2str(roundsd(pvl_dta{iR,iC},2))];
            end
            
            fcfg.out_dir = [ cfg.out_dir '/' 'plots' '/'];
            fcfg.out_nme = [ mmil_spec_char(rvl_dta{1,iC},{'.'}) '__BY__' mmil_spec_char(rvl_dta{iR,1},{'.'}) ];
            
            ejk_scatter(fcfg)
            catch fprintf('WARNING: failed scatter on %s \n',[ mmil_spec_char(rvl_dta{1,iC},{'.'}) '__BY__' mmil_spec_char(rvl_dta{iR,1},{'.'}) ]); end
        end
    end
    
end

end