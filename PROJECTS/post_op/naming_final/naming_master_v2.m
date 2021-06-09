clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Naming_final';

%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_06_03.csv';

aln_dom_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Alena_Dominance.csv';

cog_tst_nme = { 'bnt_raw_scr'     'ant_mem_raw_scr'     ...
                'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst' };
cog_tst_inc = [1 2 3 4];

dem_nme = { 'sbj_sex' 'sbj_hnd' 'sbj_edu' 'sbj_age' }; % FieldStrength
sze_nme = { 'sbj_sde_ons' 'sbj_age_ons' 'sbj_mts' 'sbj_sze_frq' 'sbj_aed_num' };
srg_nme = { 'srg_age' 'srg_sde' 'srg_typ' 'eng_out' };

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

rcn_fle = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv';

mri_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcdRSI_2020_12_17.csv';
mri_mse     = { 'subcort_vol' 'subcort_vol_ICV_cor'   };
mri_roi_use = { [ 0 ]         [ -1 ] };

dti_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/DTI_all_mmilmcdRSI_2020_12_17_v2.csv';
dti_mse     = { 'wmparc_FA_wm' 'fiber_FA'  };
dti_roi_use = { [ 2 ]          [ 0 ]       };

%% Load Data
naming_load_v2

naming_clean_v2

%% Group Comparison
naming_group_comparison_v2

%% Cross Correlations
naming_correlation_v3

naming_summarize_cc_v3
%%%%%%%%
%% Surface Correlations
naming_surface_correlations

%% Hiearchical Correlation
naming_hiearchical_correlation

%% Logistic Regression
naming_logistic_regression

%% Tables/Figures
% Table 1: Pre-operative table
Naming_Table_1

% Table 2: Post-operative table
Naming_Table_2

% Table 3: Cognitive scores
Naming_Table_3

% Table 4: Pre-operative correlations 
Naming_Table_4

% Table 5: Post-operative correlations
Naming_Table_5







