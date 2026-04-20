
%% Overall Project Directories
dte_str = '2025_03_26';

dta_dir = '/space/mcdonald-syn01/1/data/ai/';

bid_dir = '/space/mcdonald-syn01/1/BIDS/';

atl_dir = '/space/mcdonald-syn01/1/projects/ekaestner/resample_test/';

%% Defined directories
new_dta_dir = [ dta_dir '/' dte_str '/'];

cov_dir = [ new_dta_dir '/' 'covariates' '/' ];

dta_set_ifo_dir = [ new_dta_dir '/' 'dataset_info' '/' ];

%% Atlases
org_atl_dir = '/space/mcdonald-syn01/1/data/ai/original_atlases/';

src_cmd = 'source /usr/pubsw//bin/SetUpAFNI.csh 2021.01.02';

thr_dim_rsm_cmd{1} = '3dresample -rmode NN -dxyz ';
thr_dim_rsm_cmd{2} = '-input ';
thr_dim_rsm_cmd{3} = '-prefix ';

rsm_dim = [ 113 256 ];
pad_dim = { [7, 0, 7; 8, 0, 8] []};

%% Setup Run(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall Data
dta_typ_nme = { 'cat12_nii_preprocess' 'cat12_v9' };
fld_ptt = { 'anat' 'anat/mri' };
fle_ptt = { 'mwp1' 'mwp1' };

bid_nme = { 'enigma_conglom' 'ECP_dataset_bids' 'CAPES' 'MTL_proj' 'BrACE' 'Emory_UW_comb' 'ADNI-3' 'ucla_bilinguals' };

%%   ENIGMA-Epilepsy %%%%%% | bidsdir
bid_nme_use = bid_nme{1};

drv_fld{1}.(bid_nme_use) = { 'cat12_copy' 'cat12_2' 'cat12_3' 'cat12_4' };
ses_exc{1}.(bid_nme_use) = { '-02' '-FU1' '-Rescan' '-2' '-3' '-4' '-5' '-1.tar.gz' '-6' '-7' '-8' '-9' };

drv_fld{2}.(bid_nme_use) = { 'CAT12.9' };
ses_exc{2}.(bid_nme_use) = { '-02' '-FU1' '-Rescan' '-2' '-3' '-4' '-5' '-1.tar.gz' '-6' '-7' '-8' '-9' };

red_cap_fle.(bid_nme_use) = 'ENIGMAEpilepsy-CNN_DATA_2026-01-27_1529.csv';

%%   ECP %%%%%% | bidsdir
bid_nme_use = bid_nme{2};

drv_fld{1}.(bid_nme_use) = { 'cat12' };
ses_exc{1}.(bid_nme_use) = { '-1b' '-2a' '-2b' '-3a'};

drv_fld{2}.(bid_nme_use) = { 'CAT12.9' };
ses_exc{2}.(bid_nme_use) = { '-1b' '-2a' '-2b' '-3a'};

red_cap_fle.(bid_nme_use) = 'ECPFMRIData-BasicCovariates_DATA_2025-03-27_1005.csv';

%% CAPES %%%%%% | bidsdir
bid_nme_use = bid_nme{3};

drv_fld{1}.(bid_nme_use) = { 'cat12' };

drv_fld{2}.(bid_nme_use) = { 'CAT12.9' };

red_cap_fle.(bid_nme_use) = 'CAPES_Covariates_ejk_stopgap_FINAL.csv';

%% MTL %%%%%% | bidsdir
bid_nme_use = bid_nme{4};

ses_exc{1}.(bid_nme_use) = { '-capes' '-1b' };
ses_exc{2}.(bid_nme_use) = { '-capes' '-1b' };

%% BrACE %%%%%% | bidsdir
bid_nme_use = bid_nme{5};

ses_exc{1}.(bid_nme_use) = { '-2' };
ses_exc{2}.(bid_nme_use) = { '-2' };

%% UW %%%%%% | bidsdir
bid_nme_use = bid_nme{6};

ses_exc{1}.(bid_nme_use) = { '-2' '-1b' '-3'};
ses_exc{2}.(bid_nme_use) = { '-2' '-1b' '-3'};

%% HOLD / FUTURE DELETE %%%%%%
%dta_rep_nme = '_cat12_niipreprocess_t1w_data_repo'; - removal candidate

%ult_t1w_c12_dir = [ new_dta_dir '/' dta_rep_nme '/' ]; - removal candidate

%cnt_loc = '/space/mcdonald-syn01/1/software/containers/'; - removal candidate
%mni_loc = '/opt/micapipe/MNI152Volumes/MNI152_T1_0.8mm.nii.gz'; - removal candidate

%atl_ovr_nme = { 'Schaefer2018_100Parcels_17Networks_order' 'Schaefer2018_200Parcels_17Networks_order' 'aal3' 'julichbrain' }; - removal candidate

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






