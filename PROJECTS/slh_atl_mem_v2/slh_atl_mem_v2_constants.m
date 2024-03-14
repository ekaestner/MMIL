%% Project locations %%%%%%%%%%%%%%%%%%%%%%%
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Lead/slh_atl_mem_v2/';
fprintf('Project directory: %s\n',prj_dir)

dta_dir = [ prj_dir '/' 'Data' '/' ];
qal_dir = [ dta_dir '/' 'QC' '/' ];
dta_cph = 'demographic_variables.csv';

%% Demographics / Clinical
red_cap_fle = 'Redcap_2023_07_31.csv';

emy_cog_fle = 'LMforErik_2023_04_04_add_Dan.csv';
emy_dem_fle = 'demographic_variables_RR.csv';

%% Imaging
mri_fle     = 'MRI_all_2023_05_14.csv';
dti_all_fle = 'DTI_all_2023_05_14.csv';
% dti_flx_fle = 'DTI_inner_flex_2023_05_11.csv';

rcn_fle = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv';

%
mri_dev_fle     = [ dta_dir '/' mri_fle ];
dti_dev_fle     = [ dta_dir '/' dti_all_fle ];

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

mri_dev_pre     = 'MRI_dev';
mri_dev_suf     = { 'dev'                'dev'  };
mri_dev_mse     = { 'subcort_vol'        'cort_thick_ctx'  };
mri_dev_roi     = { [ 0 ]                [ 2 ] };
mri_dev_lat     = [ 1                    1 ];
mri_dev_icv     = { 'IntracranialVolume' '' };
mri_var_int     = { {'Hippocampus' 'Thalamus' 'Amygdala' } ...
                    {'fusiform' 'lateralorbitofrontal' 'entorhinal' } };

dti_dev_pre     = 'DTI_dev';
dti_dev_suf     = { 'dev'          'dev'      'dev' };
dti_dev_mse     = { 'wmparc_FA_wm' 'fiber_FA' 'gwcsurf_FA_wm_ctx' };
dti_dev_roi     = { [ 2 ]          [ 0 ]      [ 2 ] };
dti_dev_lat     = [ 1              1          1 ];
dti_dev_icv     = { ''             ''         '' };
dti_var_int     = { {'fusiform' 'lateralorbitofrontal' 'entorhinal' } ...
                    {'ILF' 'Unc' 'IFO' 'tSLF'} ...
                    {'fusiform' 'lateralorbitofrontal' 'entorhinal' } };

% 