clear; clc;

dta_loc = '/home/mmilmcdRSI/MetaData/MCD_Phenotypes/SurfGroupAvgs/output_n116/';
dta_fle = 'CortMD_diff_tval_HC_Language & memory-sm313-lh.mgh';

dta_out = '/home/ekaestne/PROJECTS/OUTPUT/mri_srf_cls/';
dta_lhs = 'CortMD_diff_tval_HC_LanguageMemory-sm313-lh.mgh';
dta_rhs = 'CortMD_diff_tval_HC_LanguageMemory-sm313-rh.mgh';

pvl_dir = '';
pvl_fle = '';

pvl_chs = .005;

%% Check cluster threshold
srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

%3549 2150
fcfg = [];

fcfg.nverts = 163842-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(2819);
fcfg.df     = 144;
fcfg.alpha  = .05;
fcfg.pval   = pvl_chs;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ... 
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );

pvl_nme = num2str(pvl_chs); pvl_nme = pvl_nme(3:end);
cls_nme = num2str(round(cls_thr));

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
    
% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% 
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
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

% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% 
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
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

% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% 
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
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

% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% 
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
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

% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% 
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster_don/';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)
