function ejk_pearson_correlation( cfg )

% cfg.sbj_nme = sbj-by-1 cell of subject names
%
% cfg.dta_one     = sbj-by-dv matrix of numbers
% cfg.dta_one_nme = 1-by-dv cell of names of the dependent variables
%
% cfg.dta_two     = sbj-by-dv matrix of numbers
% cfg.dta_two_nme = 1-by-dv cell of names of the dependent variables
%
% cfg.grp     = sbj-by-iv cell and/or matrix of group memberships
%
% cfg.out_dir = string of location to put dta
% 
% made by erik kaestner
% created: 5/11/2020 (hello from quarantine!)

%% Save Data
dta_one.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_one_nme)
    dta_one.(cfg.dta_one_nme{iD}) = cfg.dta_one(:,iD);
end

dta_two.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.dta_two_nme)
    dta_two.(cfg.dta_two_nme{iG}) = cfg.dta_two(:,iG);
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'dta_one.mat' ], 'dta_one' )
save( [ cfg.out_dir '/' 'dta_two.mat' ], 'dta_two' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_pearson_correlation.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_pearson_correlation( ''%s/dta_one.mat'', ''%s/dta_two.mat'', ''%s'' )', sve_cmd, cfg.out_dir, cfg.out_dir, cfg.out_dir);
cell2csv( [cfg.out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' cfg.out_dir '/example_R_script.r' ] );

end