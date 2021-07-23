clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'slh_atl_mem';

%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_06_22.csv';

cog_tst_nme = { 'log_mem_nor_scr_two'     'vp2_nor_scr' ...
                'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' };
cog_tst_inc = [1 2 3 4];
            
roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

mri_238_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcdRSI_2020_12_17.csv';
mri_238_mse     = { 'subcort_vol' 'subcort_vol_ICV_cor'   };
mri_238_roi_use = { [ 0 ]         [ -1 ] };

dti_238_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/';
dti_238_mse     = { 'wmparc_FA_wm' 'fiber_FA'  };
dti_238_roi_use = { [ 2 ]          [ 0 ]       };

mri_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/MRI_all_mmilmcdRSI_2020_12_17.csv';
mri_dev_mse     = { 'subcort_vol' 'subcort_vol_ICV_cor'   };
mri_dev_roi_use = { [ 0 ]         [ -1 ] };

dti_dev_fle     = '/home/ekaestne/PROJECTS/DATA/csv/ROIHOLD/';
dti_dev_mse     = { 'wmparc_FA_wm' 'wmparc_FA_wm' 'fiber_FA'  };
dti_dev_roi_use = { [ 2 ]          [ 2 ]          [ 0 ]       };

%% Load Data
slh_atl_mem_load