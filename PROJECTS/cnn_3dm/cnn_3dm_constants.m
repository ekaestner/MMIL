%%
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm';

dta_dir = [ prj_dir '/' 'Data' '/' '2023_03_27' '/'];

nii_dir = [ dta_dir '/' 'saliency_maps' '/' ];

out_dir = [ prj_dir '/' 'Output' '/' '2023_03_27' '/'];

%% Demographics
dem_dir = [ prj_dir '/' 'Data/demographics' ];

%% Performance 
dta_fll_grp_fle = 'results_Feb17.csv';
dta_sub_grp_fle = 'sampled_results.csv';



%% cnn_3dm_saliency
cat_nii     = { '^cn-ep_'          '^cn_'             '^ep_' };
cat_nii_nme = { 'cn_vs_ep'         'cn'               'ep' };
cat_col     = { rgb('bright teal') rgb('light green') rgb('maroon') };

atl_dir = [ prj_dir '/' 'Data/atlas/aal_bonilha'];
atl_nme = 'aal';

sal_grp = { '3d_model_original_data' '3d_model_ComBat' '2d_model_original_data' };
sal_nme = { 'fcn_3dm_org'            'fcn_3dm_cbt'     'fcn_2dm_org' };


