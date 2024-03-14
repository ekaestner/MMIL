%%
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Lead/cnn_3dm';

dta_dir = [ prj_dir '/' 'Data' '/' '2023_03_27' '/'];

nii_dir = [ prj_dir '/' 'Data' '/' '2023_05_15' '/' 'saliency_maps_updated' '/' ];

out_dir = [ prj_dir '/' 'Output' '/' '2023_03_27_flip_Reihaneh' '/'];

%% Demographics
dem_dir = [ prj_dir '/' 'Data/demographics' ];

%% Performance 
dta_fll_grp_fle = 'results_Feb17.csv';
dta_sub_grp_fle = 'sampled_results.csv';
dta_fll_grp_rev_fle = 'results_Jan9.csv';


%% cnn_3dm_saliency
cat_nii     = { '^cn-ep_'          '^cn_'             '^ep_' };
cat_nii_nme = { 'cn_vs_ep'         'cn'               'ep' };
cat_col     = { rgb('bright teal') rgb('light green') rgb('maroon') };
cat_typ     = [ 0                  1                  1 ]; % 0: Remove middle %, 1: Remove bottom % 

nor_dta_nme = { 'org_dta' 'zsc_dta' 'min_max_dta' };
kep_dta_nme = { 'keep' 'remove' 'remove_mid' };

low_pct = 75;

atl_dir = [ prj_dir '/' 'Data/atlas/aal_bonilha'];
atl_nme = 'aal';

sal_grp = { '3D'          '3D_combat'   '2D' }; % { '3d_model_original_data' '3d_model_ComBat' '2d_model_original_data' };
sal_nme = { 'fcn_3dm_org' 'fcn_3dm_cbt' 'fcn_2dm_org' }; %{ 'fcn_3dm_org'            'fcn_3dm_cbt'     'fcn_2dm_org' };


