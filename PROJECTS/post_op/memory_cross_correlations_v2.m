clear; clc;

dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData';
out_put = '';

%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri_dta = mmil_readtext('');



%% ORGANIZE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_col  = find(strcmpi(dem_sbj_roi_nme, 'Field Strength'));
    TS3_ind = find(cell2mat( dem_sbj_dta(:,str_col))==3);
    
sde_sze_col = find(strcmpi(dem_sbj_roi_nme, 'SideOfSeizureFocus'));
    lft_sze_ind = intersect(find(strcmpi( dem_sbj_dta(:,sde_sze_col), 'left')), TS3_ind);
    rgh_sze_ind = intersect(find(strcmpi( dem_sbj_dta(:,sde_sze_col), 'right')), TS3_ind);
    
%% RUN CROSS CORRRELATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEVHARIOBAL COVARIANCE %%%%%%%%%%%%%%%

% FA_TRACTS %%%%%%%%%%%%%%%
% LEFT SIDE
fcfg = [];

fcfg.sbj_nme = dem_sbj_sbj_nme(lft_sze_ind,1);

fcfg.dta_one = mem_scr_dta(lft_sze_ind,:);
fcfg.lbl_one = mem_scr_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = tfa_trc_dta(lft_sze_ind,:);
fcfg.lbl_two = tfa_trc_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/CrossCor/FA_Tract_by_Memory/';

ejk_cross_cor( fcfg ); 

% FA_SWM %%%%%%%%%%%%%%%
% LEFT SIDE
fcfg = [];

fcfg.sbj_nme = dem_sbj_sbj_nme(lft_sze_ind,1);

fcfg.dta_one = mem_scr_dta(lft_sze_ind,:);
fcfg.lbl_one = mem_scr_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = tfa_srf_dta(lft_sze_ind,:);
fcfg.lbl_two = tfa_srf_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/CrossCor/FA_SWM_by_Memory/';

ejk_cross_cor( fcfg );

% VOLUME %%%%%%%%%%%%%%%
% LEFT SIDE
fcfg = [];

fcfg.sbj_nme = dem_sbj_sbj_nme(lft_sze_ind,1);

fcfg.dta_one = mem_scr_dta(lft_sze_ind,:);
fcfg.lbl_one = mem_scr_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = tfa_vol_dta(lft_sze_ind,:);
fcfg.lbl_two = tfa_vol_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/CrossCor/VOLUME_by_Memory/';

ejk_cross_cor( fcfg );

%% Surface Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lod_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData/T1.csv');
        lod_dta( cellfun(@(x) strcmpi(x,'NA'),lod_dta) ) = {NaN}; 

thk_dta =  lod_dta(:,[1 28 29 31 32]);
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = thk_dta(2:end,2:end); 


dta_dir = '/home/ekaestner/Downloads/surface_correlations/JohnnyData';
% thk_lhs_srf = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_aMRI_thickness_lhs_sm176.mat');

% [ ~, ovr_lap_sbj ] = intersect( dem_sbj_sbj_nme, thk_lhs_srf.srf_dta_sbj );
% [ dem_sbj_sbj_nme(sbj_use) thk_lhs_srf.srf_dta_sbj'  ]

% % %

pvl_chs = .20;
pvl_cls = .20;
smt_stp = 176;

%
tst_use = 1;
thk_use = 1;
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

fcfg.dta_lhs = [ dta_dir '/' 'surf_aMRI_thickness_lhs_sm176.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_aMRI_thickness_rhs_sm176.mat'];

fcfg.cor     = cor_dta;
fcfg.cor_nme = cor_roi_nme;

fcfg.grp     = grp_dta;
fcfg.grp_nme = grp_roi_nme;
fcfg.grp_cmp = { '3_left' };

fcfg.cov     = cov_dta;
fcfg.cov_nme = cov_roi_nme;

fcfg.out_dir = dta_dir;

ejk_surface_correlations( fcfg );









