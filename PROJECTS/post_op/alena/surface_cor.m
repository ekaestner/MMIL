clear; clc; 

%% Surface Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lod_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData/FA_MD.csv');
        lod_dta( cellfun(@(x) strcmpi(x,'NA'),lod_dta) ) = {NaN}; 

thk_dta =  lod_dta(:,[1 56 67 109 47 48]);
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = thk_dta(2:end,2:end); 

dta_dir = '/home/ekaestner/Downloads/surface_correlations/JohnnyData';

% % %

pvl_chs = .05;
pvl_cls = .05;
smt_stp = 176;

%
tst_use = 1;
thk_use = 4;
grp_use = 4;
cov_use = 1;
sbj_use = 2:58;

% cor_dta     = mem_scr_dta(sbj_use,tst_use);
% cor_sbj_nme = mem_scr_sbj_nme(sbj_use);
% cor_roi_nme = mem_scr_roi_nme(tst_use);

cor_dta     = thk_dta(sbj_use,thk_use);
cor_sbj_nme = thk_sbj_nme(sbj_use);
cor_roi_nme = thk_roi_nme(thk_use);

grp_dta     = dem_sbj_dta(sbj_use,grp_use);
grp_sbj_nme = dem_sbj_sbj_nme(sbj_use);
grp_roi_nme = dem_sbj_roi_nme(grp_use);

cov_dta     = dem_sbj_dta(sbj_use,cov_use);
cov_sbj_nme = dem_sbj_sbj_nme(sbj_use);
cov_roi_nme = dem_sbj_roi_nme(cov_use);
    cov_roi_nme{1} = 'FieldStrength';

%
fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = cor_sbj_nme;

fcfg.dta_lhs = [ dta_dir '/' 'surf_wmparc_fa_lhs_sm313.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_wmparc_fa_rhs_sm313.mat'];

fcfg.cor     = cor_dta;
fcfg.cor_nme = cor_roi_nme;

fcfg.grp     = grp_dta;
fcfg.grp_nme = grp_roi_nme;
fcfg.grp_cmp = { '3_left' };

fcfg.cov     = cov_dta;
fcfg.cov_nme = cov_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006';

ejk_surface_correlations( fcfg );

%% Memory

tst_use = 8;
cor_dta     = mem_scr_dta(sbj_use,tst_use);
cor_sbj_nme = mem_scr_sbj_nme(sbj_use);
cor_roi_nme = mem_scr_roi_nme(tst_use);

fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = cor_sbj_nme;

fcfg.dta_lhs = [ dta_dir '/' 'surf_wmparc_fa_lhs_sm313.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_wmparc_fa_rhs_sm313.mat'];

fcfg.cor     = cor_dta;
fcfg.cor_nme = cor_roi_nme;

fcfg.grp     = grp_dta;
fcfg.grp_nme = grp_roi_nme;
fcfg.grp_cmp = { 'left' };

fcfg.cov     = []; %cellfun(@num2str,cov_dta,'UniformOutput',false);
fcfg.cov_nme = ''; cov_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_January';
fcfg.out_pre = 'no_cov_n18';

ejk_surface_correlations( fcfg );


