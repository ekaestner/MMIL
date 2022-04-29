clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'slh_atl_mem';

%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_06_22.csv';

rcn_fle = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv';

cog_tst_nme = { 'log_mem_nor_scr_two'     'vp2_nor_scr' ...
                'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' };
            
roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

mri_238_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcdRSI_2020_12_17.csv';
mri_238_pre     = 'MRI_238';
mri_238_suf     = { '238'                '238'  };
mri_238_mse     = { 'subcort_vol'        'cort_thick_ctx'  };
mri_238_roi     = { [ 0 ]                [2] };
mri_238_lat     = [ 1                    1 ];
mri_238_icv     = { 'IntracranialVolume' '' };

dti_238_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/DTI_all_mmilmcdRSI_2020_12_17_v2.csv';
dti_238_pre     = 'DTI_238';
dti_238_suf     = { '238'                '238'  };
dti_238_mse     = { 'wmparc_FA_wm' 'fiber_FA'  };
dti_238_roi     = { [ 2 ]          [ 0 ]       };
dti_238_lat     = [ 1                1 ];
dti_238_icv     = { ''             '' };

mri_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcd2_2021_08_03.csv';
mri_dev_pre     = 'MRI_dev';
mri_dev_suf     = { 'dev'                'dev'  };
mri_dev_mse     = { 'subcort_vol'        'cort_thick_ctx'  };
mri_dev_roi     = { [ 0 ]                [ 2 ] };
mri_dev_lat     = [ 1                    1 ];
mri_dev_icv     = { 'IntracranialVolume' '' };

dti_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/DTI_all_mmilmcd2_2021_08_03.csv';
dti_dev_pre     = 'DTI_dev';
dti_dev_suf     = { 'dev'                'dev'  };
dti_dev_mse     = { 'wmparc_FA_wm' 'fiber_FA'  };
dti_dev_roi     = { [ 2 ]          [ 0 ]       };
dti_dev_lat     = [ 1              1 ];
dti_dev_icv     = { ''             '' };

%% Load Data
slh_atl_mem_load_v2

slh_atl_mem_clean

%% Groups
slh_atl_mem_groups

%% Cross Correlations
% Clinical
slh_atl_mem_correlation_clinical

% Cognitive
slh_atl_mem_correlation_cognitive

% Neurobio between groups
slh_atl_mem_comparison_neurobio

%% Misc
Presentation_prep

AlenaCompare_Tables

rvalue_examine

slh_atl_phenotype_explore

