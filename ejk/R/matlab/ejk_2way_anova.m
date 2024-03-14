function ejk_2way_anova( cfg )

%% Save Data
dep_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    try dep_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD); catch fprintf('Failure Field Name: %s \n',cfg.dta_nme{iD}); end
end

grp_var_one.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.grp_nme_one)
    grp_var_one.(cfg.grp_nme_one{iG}) = cfg.grp_one(:,iG);
end

grp_var_two.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.grp_nme_two)
    grp_var_two.(cfg.grp_nme_two{iG}) = cfg.grp_two(:,iG);
end


ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dep_var.mat' ], 'dep_var' )
save( [ cfg.out_dir '/' 'grp_var_one.mat' ], 'grp_var_one' )
save( [ cfg.out_dir '/' 'grp_var_two.mat' ], 'grp_var_two' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_2way_anova.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_2way_anova( ''%s/dep_var.mat'', ''%s/grp_var_one.mat'', ''%s/grp_var_two.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir, cfg.out_dir);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' mmil_spec_char(cfg.out_dir,{' '},{'\ '}) '/example_R_script.r' ] );

end