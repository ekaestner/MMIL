clear; clc;

%% Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_dir = '/home/ekaestner/Downloads/surface_correlations/JohnnyData';

lod_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData' '/' 'Memory_FA_Updated_Raw.csv' ] );
    lod_dta( cellfun(@isempty,lod_dta) ) = {NaN}; 
    
% DEM_SBJ
dem_sbj_dta =  lod_dta(:,[1 29 3 4:16 30 ]);
dem_sbj_sbj_nme = dem_sbj_dta(2:end,1);
dem_sbj_roi_nme = dem_sbj_dta(1,2:end);
dem_sbj_dta     = dem_sbj_dta(2:end,2:end); 

% pre_mem
pre_mem_dta =  lod_dta(:,[1 17:20]);
pre_mem_sbj_nme = pre_mem_dta(2:end,1);
pre_mem_roi_nme = pre_mem_dta(1,2:end);
pre_mem_dta     = pre_mem_dta(2:end,2:end); 

% pst_mem
pst_mem_dta =  lod_dta(:,[1 31 21:24 32]);
pst_mem_sbj_nme = pst_mem_dta(2:end,1);
pst_mem_roi_nme = pst_mem_dta(1,2:end);
pst_mem_dta     = pst_mem_dta(2:end,2:end); 

% roi_nme
roi_nme_dta =  lod_dta(:,[1 25:28]);
roi_nme_sbj_nme = roi_nme_dta(2:end,1);
roi_nme_roi_nme = roi_nme_dta(1,2:end);
roi_nme_dta     = roi_nme_dta(2:end,2:end); 

%% Include dominance
dom_col = 10;
rep_col = 1:2;

for iC = 1:numel(rep_col)
    dem_sbj_dta(:,rep_col(iC)) = strcat(dem_sbj_dta(:,rep_col(iC)),'_',dem_sbj_dta(:,dom_col));    
end

% for iC = 1:numel(rep_col)
%     dem_sbj_dta(:,rep_col(iC)) = dem_sbj_dta(:,dom_col);    
% end


%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_nme     = { 'Entorhinal' 'LDFR' 'LM2' 'VPA2' 'BVMT' 'EntorhinalRight' };
cov_nme     = { 'none' 'scannerstrength' 'presurgical' };
    cov_tst_num = [2 1 2 3 4 2];
sbj_num_nme = { 'n18'    'n21'  };
       grp_nme = { '3_left_dominant' 'left_dominant' }; % { 'dominant' 'dominant' }; % { '3_left_dominant' 'left_dominant' }; % grp_nme = { '3_left' 'left' };

%% Run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst_num_use = 4; % [ 1 2 3 4 5 ]
cov_typ_use = 1; % [ 1 2 3 ]
sbj_num_use = 2; % [ 1 2 ]

pvl_chs = .05;
pvl_cls = .05;
smt_stp = 176; % [ 176 313 ]

% Subject Names
sbj_nme_use = dem_sbj_sbj_nme;

% Data
cor_dta     = pst_mem_dta(:, tst_num_use);
cor_roi_nme = dta_nme(tst_num_use);

% Covariates
if cov_typ_use==1
    cov         = [];
    cov_nme_hld = {''};
elseif cov_typ_use==2
    cov         = cellfun( @num2str, dem_sbj_dta(:,3), 'uni', 0);
    cov_nme_hld = {'Strength'};
elseif cov_typ_use==3
    cov         = pre_mem_dta(:,cov_tst_num(tst_num_use));
    cov_nme_hld = pre_mem_roi_nme(:,cov_tst_num(tst_num_use));
end

% Group Names
grp_dta     = dem_sbj_dta(:, sbj_num_use);
grp_roi_nme = dem_sbj_roi_nme(sbj_num_use);
grp_cmp     = grp_nme(sbj_num_use);

%
out_pre_fix = [ cor_roi_nme{1} '_' sbj_num_nme{sbj_num_use} '_' 'Covariate' '_' cov_nme{cov_typ_use} '_' grp_nme{sbj_num_use} ];

%% Plot
fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = sbj_nme_use;

fcfg.dta_lhs = [ dta_dir '/' 'surf_wmparc_fa_lhs_sm' num2str(smt_stp) '_epd006switch.mat']; % [ dta_dir '/' 'surf_wmparc_fa_lhs_sm' num2str(smt_stp) '_epd006switch_epd096add.mat']
fcfg.dta_rhs = [ dta_dir '/' 'surf_wmparc_fa_rhs_sm' num2str(smt_stp) '_epd006switch.mat']; % [ dta_dir '/' 'surf_wmparc_fa_rhs_sm' num2str(smt_stp) '_epd006switch_epd096add.mat']

fcfg.cor     = cor_dta;
fcfg.cor_nme = cor_roi_nme;

fcfg.grp     = grp_dta;
fcfg.grp_nme = grp_roi_nme;
fcfg.grp_cmp = grp_cmp;

fcfg.cov     = cov; 
fcfg.cov_nme = cov_nme_hld; 

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February';
fcfg.out_pre = [out_pre_fix '_' 'ATL_only_Raw'];

ejk_surface_correlations( fcfg );
