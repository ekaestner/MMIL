function ejk_fisher_test( cfg )

%%
dta_one.sbj_nme = cfg.sbj;
for iD = 1:numel(cfg.lbl_one)
    dta_one.(cfg.lbl_one{iD}) = cfg.dta_one(:,iD);
    if isfield(cfg,'lvl')
        use_lvl = regexp(cfg.lvl{iD},'/','split');
        dta_one.(cfg.lbl_one{iD})(~ismember(dta_one.(cfg.lbl_one{iD}),use_lvl)) = {'NA'};
    end
end

dta_two.sbj_nme   = cfg.sbj;
for iG = 1:numel(cfg.lbl_two)
    dta_two.(cfg.lbl_two{iG}) = cfg.dta_two(:,iG);
end

cfg.out_dir = [ cfg.out_dir '/' cfg.grp_nme];

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_one.mat' ], 'dta_one' )
save( [ cfg.out_dir '/' 'dta_two.mat' ], 'dta_two' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_fisher_test.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_fisher_test( ''%s/dta_one.mat'', ''%s/dta_two.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' mmil_spec_char(cfg.out_dir,{' '},{'\ '}) '/example_R_script.r' ] );


end