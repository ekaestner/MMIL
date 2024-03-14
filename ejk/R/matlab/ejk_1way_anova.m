function ejk_1way_anova( cfg )

% cfg.sbj_nme = sbj-by-1 cell of subject names
%
% cfg.dta     = sbj-by-dv matrix of numbers
% cfg.dta_nme = 1-by-dv cell of names of the dependent variables
%
% cfg.grp     = sbj-by-iv cell and/or matrix of group memberships
%
% cfg.out_dir = string of location to put dta
% 
% made by erik kaestner
% created: 5/11/2020 (hello from quarantine!)

%% Save Data
dep_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    try dep_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD); catch fprintf('Failure Field Name: %s \n',cfg.dta_nme{iD}); end
end

grp_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.grp_nme)
    grp_var.(cfg.grp_nme{iG}) = cfg.grp(:,iG);
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dep_var.mat' ], 'dep_var' )
save( [ cfg.out_dir '/' 'grp_var.mat' ], 'grp_var' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_1way_anova.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_1way_anova( ''%s/dep_var.mat'', ''%s/grp_var.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' mmil_spec_char(cfg.out_dir,{' '},{'\ '}) '/example_R_script.r' ] );

end