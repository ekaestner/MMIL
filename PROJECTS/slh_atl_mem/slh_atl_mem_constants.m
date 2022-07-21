% Project locations %%%%%%%%%%%%%%%%%%%%%%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'slh_atl_mem';
dta_dir = '/home/ekaestne/PROJECTS/DATA';

% File locations %%%%%%%%%%%%%%%%%%%%%%%
% Clinical & Cognitive data location
red_cap_fle = [ dta_dir '/' 'csv' '/' 'redcap' '/' 'Redcap_2022_01_10.csv'];

aln_dom_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Alena_Dominance.csv';

% Define measures of interest %%%%%%%%%%%%%%%%%%%%%%%
% Clinical
dem_nme = { 'sbj_sex' 'sbj_hnd' 'sbj_edu' 'sbj_age' }; % FieldStrength
sze_nme = { 'sbj_sde_ons' 'sbj_age_ons' 'sbj_mts' 'sbj_sze_frq' 'sbj_aed_num' };
srg_nme = { 'srg_age' 'srg_sde' 'srg_typ' 'eng_out' };

% Cognitive
cog_tst_nme = { 'log_mem_nor_scr_two'     'vp2_nor_scr'     'cvl_lfr_nor_scr' ... % LM = SS, VP = SS, CVLT = Tscore
                'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' 'cvl_lfr_nor_scr_pst' };

% Neurobiological
mri_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcd2_2021_08_03.csv';
dti_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/DTI_all_mmilmcd2_2021_08_03.csv';

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

rcn_fle = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv';

mri_dev_pre     = 'MRI_dev';
mri_dev_suf     = { 'dev'                'dev'  };
mri_dev_mse     = { 'subcort_vol'        'cort_thick_ctx'  };
mri_dev_roi     = { [ 0 ]                [ 2 ] };
mri_dev_lat     = [ 1                    1 ];
mri_dev_icv     = { 'IntracranialVolume' '' };

dti_dev_pre     = 'DTI_dev';
dti_dev_suf     = { 'dev'          'dev'      'dev' };
dti_dev_mse     = { 'wmparc_FA_wm' 'fiber_FA' 'gwcsurf_FA_wm' };
dti_dev_roi     = { [ 2 ]          [ 0 ]      [ 0 ] };
dti_dev_lat     = [ 1              1          1 ];
dti_dev_icv     = { ''             ''         '' };