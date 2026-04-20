
%% Overall Project Directories
dte_str = '2026_06_01';

dta_dir = '/space/mcdonald-syn01/1/data/ai/';

bid_dir = '/space/mcdonald-syn01/1/BIDS/';

atl_dir = '/space/mcdonald-syn01/1/projects/ekaestner/resample_test/';

%% Defined directories
new_dta_dir = [ dta_dir '/' dte_str '/'];

cov_dir = [ new_dta_dir '/' 'covariates' '/' ];

dta_set_ifo_dir = [ new_dta_dir '/' 'dataset_info' '/' ];

typ_nme = {'cmig_struct_func_session1_run1'};
raw_fle_nme = { {'t1w_raw' 'func_raw'} };
dta_typ_nme = { {{'mwp1_processed'} {'reho_postproc' 'alff_postproc' 'falff_postproc'}} };

%% Atlases
org_atl_dir = '/space/mcdonald-syn01/1/data/ai/original_atlases/';

src_cmd = 'source /usr/pubsw//bin/SetUpAFNI.csh 2021.01.02';

thr_dim_rsm_cmd{1} = '3dresample -rmode NN -dxyz ';
thr_dim_rsm_cmd{2} = '-input ';
thr_dim_rsm_cmd{3} = '-prefix ';

rsm_dim = [ 113 256 ];
pad_dim = { [7, 0, 7; 8, 0, 8] []};

%% Biovariates
drv_fld = { 'cat12-r2166' 'freesurfer_7.3.2_MCRv79_mcrLinked'  }; % 'mriqc_v24.0.2'
drv_mes = { { 'aggregate_biovariates_aparc_DK40_wide.csv' 'aggregate_report_QC_qualitymeasures.csv' } ...
            { 'crossSectional/euler_crossSectional.csv' ...
              'crossSectional/aparc.DKTatlas_crossSectional_wideTable.csv' ...
              'crossSectional/aseg_crossSectional_wideTable.csv' }  };

%% Setup Run(s) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall Data
bid_nme = { 'ENIGMA_dataset' 'ECP_dataset_bids' 'CAPES' 'IDEAS_dataset' 'HEP_dataset_bids' 'DLBS_dataset' 'MTL_proj' 'VU_dataset' }; % 'ADNI-3'

%%   ENIGMA-Epilepsy %%%%%% | bidsdir
bid_nme_use = bid_nme{1};

drv_dte.(bid_nme_use) = { '20260322' '' };

red_cap_fle.(bid_nme_use) = 'ENIGMAEpilepsy-CNN_DATA_2026-03-24_0848.csv';

%%   ECP %%%%%% | bidsdir
bid_nme_use = bid_nme{2};

drv_dte.(bid_nme_use) = { '20260318' '' };

red_cap_fle.(bid_nme_use) = 'ECPFMRIData-BasicCovariates_DATA_2025-03-27_1005.csv';

%% CAPES %%%%%% | bidsdir
bid_nme_use = bid_nme{3};

drv_dte.(bid_nme_use) = { '20260318' '20260113' };

red_cap_fle.(bid_nme_use) = 'CAPES_Covariates_ejk_stopgap_FINAL.csv';

%% IDEAS %%%%%% | bidsdir
bid_nme_use = bid_nme{4};

drv_dte.(bid_nme_use) = { '20260204' '20260204' };

red_cap_fle.(bid_nme_use) = 'Metadata_IDEAS_Combined.csv';

%% HEP %%%%%% | bidsdir
bid_nme_use = bid_nme{5};

drv_dte.(bid_nme_use) = { '20260317' '20260317' };

red_cap_fle.(bid_nme_use) = '';

%% DLBS_dataset %%%%%% | bidsdir
bid_nme_use = bid_nme{6};

drv_dte.(bid_nme_use) = { '' '' };

red_cap_fle.(bid_nme_use) = '';

%% MTL_proj %%%%%% | bidsdir
bid_nme_use = bid_nme{7};

drv_dte.(bid_nme_use) = { '' '' };

red_cap_fle.(bid_nme_use) = '';

%% VU_dataset %%%%%% | bidsdir
bid_nme_use = bid_nme{8};

drv_dte.(bid_nme_use) = { '' '' };

red_cap_fle.(bid_nme_use) = '';



