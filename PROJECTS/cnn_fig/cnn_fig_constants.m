
%% Overall Project Directory
prj_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/';

dta_dir = '/space/mcdonald-syn01/1/data/ai/';

bid_dir = '/space/mcdonald-syn01/1/BIDS/';

atl_dir = '/space/mcdonald-syn01/1/projects/ekaestner/resample_test/';

%% Data locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%   shared directories %%%%%%
dte_str = '2024_11_22';

new_dta_dir = [ prj_dir '/' 'new_data' '/'];

ult_t1w_c12_dir = [ new_dta_dir '/' '_cat12_t1w_data_repo' '/' ];

cov_dir = [ new_dta_dir '/' 'covariates' '/' ];

dta_set_ifo_dir = [ new_dta_dir '/' 'dataset_info' '/' ];

cnt_loc = '/space/mcdonald-syn01/1/software/containers/';
mni_loc = '/opt/micapipe/MNI152Volumes/MNI152_T1_0.8mm.nii.gz';

bid_nme = { 'enigma_conglom' 'ECP_dataset_bids' 'CAPES' 'MTL_proj' 'BrACE' 'Emory_UW_comb' 'ADNI-3' 'ucla_bilinguals' };

atl_ovr_nme = { 'Schaefer2018_100Parcels_17Networks_order' 'Schaefer2018_200Parcels_17Networks_order' 'aal3' 'julichbrain' };

%%   ENIGMA-Epilepsy %%%%%%
bid_nme_use = bid_nme{1};

ses_exc.(bid_nme_use) = { '-02' '-FU1' '-Rescan' '-2' '-3' '-4' '-5' '-1.tar.gz' '-6' '-7' '-8' '-9' };

drv_fld.(bid_nme_use) = { 'cat12' 'cat12_2' 'cat12_3' };

red_cap_fle.(bid_nme_use) = 'ENIGMAEpilepsy-CNN_DATA_2024-10-23_1021.csv';

%%   ECP %%%%%%
bid_nme_use = bid_nme{2};

ses_exc.(bid_nme_use) = { '-1b' '-2a' '-2b' '-3a'};

drv_fld.(bid_nme_use) = { 'cat12' };

red_cap_fle.(bid_nme_use) = 'ECP_Covariates_edited.csv';

%% CAPES %%%%%%
bid_nme_use = bid_nme{3};

drv_fld.(bid_nme_use) = { 'cat12' };

red_cap_fle.(bid_nme_use) = 'CAPES_Covariates_ejk_stopgap.csv';

%% MTL %%%%%%
bid_nme_use = bid_nme{4};

ses_exc.(bid_nme_use) = { '-capes' '-1b' };

%% BrACE %%%%%%
bid_nme_use = bid_nme{5};

ses_exc.(bid_nme_use) = { '-2' };

%% UW %%%%%%
bid_nme_use = bid_nme{6};

ses_exc.(bid_nme_use) = { '-2' '-1b' '-3'};

%% SEE-GAAN %%%%%%%
mdl_fld = [ prj_dir '/' 'see_gaan' '/'];
mdl_nme = { 'ejk_oct_2024' 'Animesh_version'  };

%% HOLD / FUTURE DELETE %%%%%%
% eng_ses_hld = mmil_readtext([ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'enigma_session.csv']);
% red_cap_fle = 'ENIGMAEpilepsy-CNN_DATA_2024-05-24_1338_comma.csv';
% t1w_fsr_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'freesurfer' '/' 't1'];
% t1w_mca_dir = [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'micapipe' '/' 't1' '/'];
% mca_drv_dir = [ eng_dir '/' 'derivatives' '/' 'micapipe_v0.2.0' '/' ];
% t1w_mca_drv_txt = [  'ses-01' '/' 'anat' '/'];
% t1w_mca_drv_txt_no0 = [  'ses-1' '/' 'anat' '/'];
% t1w_mca_drv_txt_nse = [  'anat' '/'];
% t1w_mca_drv_fle     = '_ses-01_space-fsnative_T1w.nii.gz';
% t1w_mca_drv_fle_no0 = '_ses-1_space-fsnative_T1w.nii.gz';
% t1w_mca_drv_fle_nse = '_space-fsnative_T1w.nii.gz';
 
% eng_c12_drv_dir = [ bid_dir '/' 'enigma_conglom' '/' 'derivatives' '/' ];
% c12_drv_ver_one_dir = [ eng_c12_drv_dir '/' 'cat12' '/' ];
% c12_drv_ver_two_dir = [ eng_c12_drv_dir '/' 'cat12_2' '/' ];
% 
% t1w_c12_eng_dir = [ new_dta_dir '/' 'enigma_conglom' '/' 'cat12' '/' 't1' '/'];
%     
% ecp_c12_drv_dir = [ bid_dir '/' 'ECP_dataset_bids' '/' 'derivatives' '/' 'cat12' '/'];
% 
% t1w_c12_ecp_dir = [ prj_dir '/' 'new_data' '/' 'ecp' '/' 'cat12' '/' 't1' '/'];
% red_cap_ecp_dir = [ prj_dir '/' 'new_data' '/' 'ecp' '/' 'covariates' '/'];






