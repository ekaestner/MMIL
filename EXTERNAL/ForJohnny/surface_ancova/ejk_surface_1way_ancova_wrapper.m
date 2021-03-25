%% ejk_surface_1way_ancova %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_chs = .05;
pvl_cls = .05;
smt_stp = 176;

%
dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA';

new_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/epd_age_covariates_no_nan_no_mass.csv');

grp_dta     = new_dta(:,[1 6]);
grp_sbj_nme = grp_dta(2:end,1);
grp_roi_nme = grp_dta(1,2:end);
grp_dta     = grp_dta(2:end,2:end);

cov_dta     = new_dta(:,[1 2 3 4 8]);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%
fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = grp_sbj_nme;

fcfg.dta_lhs = [ dta_dir '/' 'surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_aMRI_thickness_rhs_sm176_no_nan_no_mass.mat'];

fcfg.grp     = grp_dta;
fcfg.grp_cmp = { { { 'TLE' 'HC' } { 'MCI' 'HC' } {'TLE' 'MCI'} } };
fcfg.grp_nme = grp_roi_nme;

fcfg.cov     = cov_dta;
fcfg.cov_nme = cov_roi_nme;

fcfg.out_dir = '/home/ekaestner/Downloads/surface_ancova';

ejk_surface_1way_ancova( fcfg );