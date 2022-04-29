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
