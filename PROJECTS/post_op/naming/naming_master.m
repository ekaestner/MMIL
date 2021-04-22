clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Naming';

%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_02_24.csv';

aln_dom_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Alena_Dominance.csv';

cog_tst_nme = { 'bnt_raw_scr'     'ant_mem_raw_scr'     'cat_flu_nor_scr'     ...
                'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst' 'cat_flu_nor_scr_pst' };
cog_tst_inc = [1 2 4 5];

dem_nme = { 'sbj_sex' 'sbj_hnd' 'sbj_edu' 'sbj_age' }; % FieldStrength
sze_nme = { 'sbj_sde_ons' 'sbj_age_ons' 'sbj_mts' 'sbj_sze_frq' 'sbj_aed_num' };
srg_nme = { 'srg_age' 'srg_sde' 'srg_typ' 'eng_out' };

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

rcn_fle = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv';

fmr_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/fMRI_N_FF_nzvoxels_2019_05_06.csv';
fmr_mse     = { 'destr' 'alicia' };
fmr_roi_use = { [0]     [0] };

mri_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcdRSI_2020_12_17.csv';
mri_mse     = { 'cort_thick_ctx' 'subcort_vol_ICV_cor' };
mri_roi_use = { [ 2 ]            [ 0 ] };

dti_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/DTI_all_mmilmcdRSI_2020_12_17_v2.csv';
dti_mse     = { 'wmparc_FA_wm' 'fiber_FA'  }; % 'wmparc_MD_wm' 'fiber_MD' 'aseg_FA' 'aseg_MD'
dti_roi_use = { [ 2 ]          [ 0 ]       }; % [ 2 ]          [ 0 ]      [0]       [0] 

rsf_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/rsBOLD_fparc_var_2020_12_17.csv';
rsf_mse     = { 'var_ctx' 'var_vol' 'network' };
rsf_roi_use = { [ 2 ]     [ 0 ]     [ 0 ]      };

%% Load Data
naming_load

naming_clean

%% Group Comparison
naming_group_comparison

%% Cross Correlations
naming_correlation_v2

naming_summarize_cc_v3

%% Surface Correlations
naming_surface_correlations

%% Hiearchical Correlation
naming_hiearchical_correlation

%% Logistic Regression
naming_logistic_regression

%% Tables/Figures


% Table 2: Cognitive
Naming_Table_2










