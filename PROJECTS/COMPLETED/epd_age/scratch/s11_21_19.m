clear; clc;

dta_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/11_21_19/';

pvl_chs = .01;
emn_chs = .05;

%% Test ANCOVA_TLE_vs_HC
plt_nme = 'ANCOVA_TLE_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;
 
% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 

%% Test ANCOVA_MCI_vs_HC
plt_nme = 'ANCOVA_MCI_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_MCI_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_MCI_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_MCI_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean*-1;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_MCI_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean*-1;
    
% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 

%% Test ANCOVA_TLE_vs_MCI
plt_nme = 'ANCOVA_TLE_vs_MCI';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_MCI.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_MCI.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_MCI.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_MCI.mat');
    men_dff_rhs = men_dff_rhs.emmean;
   
% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)    
    
%% Test ANCOVA_LATE_vs_HC
plt_nme = 'ANCOVA_LATE_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Late_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Late_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Late_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Late_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 

%% Test ANCOVA_EARLY_vs_HC
plt_nme = 'ANCOVA_EARLY_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Early_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Early_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Early_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Early_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 

%% Check Vertices
% Data
sbj_nme     = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/epd_age_covariates_no_nan.csv');
    sbj_nme(cellfun(@isempty,sbj_nme(:,10)),10) = {''};
ely_ind = find(strcmpi(sbj_nme(:,10),'Early'))-1;
lte_ind = find(strcmpi(sbj_nme(:,10),'Late'))-1;
con_ind = find(strcmpi(sbj_nme(:,6),'HC'))-1;
mci_ind = find(strcmpi(sbj_nme(:,6),'MCI'))-1;
epd_ind = find(strcmpi(sbj_nme(:,6),'EPD_Old'))-1;

org_thk_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_no_nan.mat');
org_thk_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_rhs_no_nan.mat');

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Late_HC.mat');
    pvl_lhs_lte = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Late_HC.mat');
    pvl_rhs_lte = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Late_HC.mat');
    men_dff_lhs_lte = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Late_HC.mat');
    men_dff_rhs_lte = men_dff_rhs.emmean;

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Early_HC.mat');
    pvl_lhs_ely = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Early_HC.mat');
    pvl_rhs_ely = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Early_HC.mat');
    men_dff_lhs_ely = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Early_HC.mat');
    men_dff_rhs_ely = men_dff_rhs.emmean;    

% LEFT / THIN
[ ~ , thn_vtx_lhs ] = sort(men_dff_lhs_lte); %%%
thn_vtx_lhs = thn_vtx_lhs(20);

figure() %%%
subplot(4,1,1)
hist( org_thk_lhs.srf_dta(con_ind,thn_vtx_lhs) , linspace(1,4,100) )
ylabel('CONTROL'); xlim([1 4]); ylim([0 10]);
subplot(4,1,2)
hist( org_thk_lhs.srf_dta(mci_ind,thn_vtx_lhs) , linspace(1,4,100) )
ylabel('MCI'); xlim([1 4]); ylim([0 10]);
subplot(4,1,3)
hist( org_thk_lhs.srf_dta(ely_ind,thn_vtx_lhs) , linspace(1,4,100) )
ylabel('EARLY'); xlim([1 4]); ylim([0 10]);
subplot(4,1,4);
hist( org_thk_lhs.srf_dta(lte_ind,thn_vtx_lhs) , linspace(1,4,100) )
ylabel('LATE'); xlim([1 4]); ylim([0 10]);

tightfig();
print(gcf,[dta_out '/' 'lhs_thn' '.png'],'-dpng','-r200')
close all

% LEFT / THICK
[ ~ , thk_vtx_lhs ] = max(men_dff_lhs_ely); %%%

figure() %%%
subplot(4,1,1)
hist( org_thk_lhs.srf_dta(con_ind,thk_vtx_lhs) , linspace(1,4,100) )
ylabel('CONTROL'); xlim([1 4]); ylim([0 10]);
subplot(4,1,2)
hist( org_thk_lhs.srf_dta(mci_ind,thk_vtx_lhs) , linspace(1,4,100) )
ylabel('MCI'); xlim([1 4]); ylim([0 10]);
subplot(4,1,3)
hist( org_thk_lhs.srf_dta(ely_ind,thk_vtx_lhs) , linspace(1,4,100) )
ylabel('EARLY'); xlim([1 4]); ylim([0 10]);
subplot(4,1,4);
hist( org_thk_lhs.srf_dta(lte_ind,thk_vtx_lhs) , linspace(1,4,100) )
ylabel('LATE'); xlim([1 4]); ylim([0 10]);

tightfig();
print(gcf,[dta_out '/' 'lhs_thk' '.png'],'-dpng','-r200')
close all

% RIGHT / THIN
[ ~ , thn_vtx_rhs ] = sort(men_dff_rhs_lte); %%%
thn_vtx_rhs = thn_vtx_rhs(20);

figure() %%%
subplot(4,1,1)
hist( org_thk_rhs.srf_dta(con_ind,thn_vtx_rhs) , linspace(1,4,100) )
ylabel('CONTROL'); xlim([1 4]); ylim([0 10]);
subplot(4,1,2)
hist( org_thk_rhs.srf_dta(mci_ind,thn_vtx_rhs) , linspace(1,4,100) )
ylabel('MCI'); xlim([1 4]); ylim([0 10]);
subplot(4,1,3)
hist( org_thk_rhs.srf_dta(ely_ind,thn_vtx_rhs) , linspace(1,4,100) )
ylabel('EARLY'); xlim([1 4]); ylim([0 10]);
subplot(4,1,4);
hist( org_thk_rhs.srf_dta(lte_ind,thn_vtx_rhs) , linspace(1,4,100) )
ylabel('LATE'); xlim([1 4]); ylim([0 10]);

tightfig();
print(gcf,[dta_out '/' 'rhs_thn' '.png'],'-dpng','-r200')
close all

% RIGHT / THICK
[ ~ , thk_vtx_rhs ] = max(men_dff_rhs_ely); %%%

figure() %%%
subplot(4,1,1)
hist( org_thk_rhs.srf_dta(con_ind,thk_vtx_rhs) , linspace(1,4,100) )
ylabel('CONTROL'); xlim([1 4]); ylim([0 10]);
subplot(4,1,2)
hist( org_thk_rhs.srf_dta(mci_ind,thk_vtx_rhs) , linspace(1,4,100) )
ylabel('MCI'); xlim([1 4]); ylim([0 10]);
subplot(4,1,3)
hist( org_thk_rhs.srf_dta(ely_ind,thk_vtx_rhs) , linspace(1,4,100) )
ylabel('EARLY'); xlim([1 4]); ylim([0 10]);
subplot(4,1,4);
hist( org_thk_rhs.srf_dta(lte_ind,thk_vtx_rhs) , linspace(1,4,100) )
ylabel('LATE'); xlim([1 4]); ylim([0 10]);

tightfig();
print(gcf,[dta_out '/' 'rhs_thk' '.png'],'-dpng','-r200')
close all

fcfg = []; %%%

fcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/11_21_19/';
fcfg.out_pre_fix = 'vertex_highlights';

fcfg.vtx_cor = { [ thn_vtx_lhs        thk_vtx_lhs ]      [ thn_vtx_rhs        thk_vtx_rhs ]};
fcfg.vtx_col = { { rgb('bright blue') rgb('bright red')} { rgb('bright blue') rgb('bright red')} };

ejk_highlight_vertex(fcfg)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test ANCOVA_LATE_vs_HC_3T
plt_nme = 'ANCOVA_LATE_vs_HC_3T';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Late_HC_3T.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Late_HC_3T.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Late_HC_3T.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Late_HC_3T.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 

%% Test ANCOVA_EARLY_vs_HC_3T
plt_nme = 'ANCOVA_EARLY_vs_HC_3T';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Early_HC_3T.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Early_HC_3T.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Early_HC_3T.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Early_HC_3T.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% EM-Thresh
lhs_bdd_ind = find(abs(men_dff_lhs)<emn_chs);
rhs_bdd_ind = find(abs(men_dff_rhs)<emn_chs);

men_dff_lhs_emn = men_dff_lhs;
men_dff_rhs_emn = men_dff_rhs;

men_dff_lhs_emn(lhs_bdd_ind) = 0;
men_dff_rhs_emn(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff'];

pcfg.plt_dta = { men_dff_lhs_emn' men_dff_rhs_emn' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% pvalue-Thresh
lhs_bdd_ind = find(pvl_lhs>pvl_chs);
rhs_bdd_ind = find(pvl_rhs>pvl_chs);

men_dff_lhs_pvl = men_dff_lhs;
men_dff_rhs_pvl = men_dff_rhs;

men_dff_lhs_pvl(lhs_bdd_ind) = 0;
men_dff_rhs_pvl(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = dta_out;
pcfg.out_pre_fix = [ plt_nme '_estmean_diff_p01'];

pcfg.plt_dta = { men_dff_lhs_pvl' men_dff_rhs_pvl' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.05 0.05 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg) 











