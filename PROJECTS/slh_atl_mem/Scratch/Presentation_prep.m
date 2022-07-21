out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore';

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

%% Subject #'s
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_dev_QC.csv'];
fcfg.dta_col = 2;
[ dti_dev_dta, dti_dev_dta_sbj, dti_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'cort_thick_ctx_aparc_annot_dev_QC.csv'];
fcfg.dta_col = 2;
[ mri_dev_dta, mri_dev_dta_sbj, mri_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_238_QC.csv'];
fcfg.dta_col = 2;
[ dti_238_dta, dti_238_dta_sbj, dti_238_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'cort_thick_ctx_aparc_annot_238_QC.csv'];
fcfg.dta_col = 2;
[ mri_238_dta, mri_238_dta_sbj, mri_238_dta_col] = ejk_dta_frm( fcfg );

% Pre LM2 & DTI
fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_3T_allSurg_all,1)) & ~cellfun(@isnan, dti_dev_dta( grp.controls_pre_3T_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_left,1))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_3T_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_right,1))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_3T_allSurg_right,5)) ) );

fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_allT_allSurg_all,1)) & ~cellfun(@isnan, dti_dev_dta( grp.controls_pre_allT_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_left,1))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_allT_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_right,1))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_allT_allSurg_right,5)) ) );

% Pre LM2 & MRI
fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_3T_allSurg_all,1)) & ~cellfun(@isnan, mri_dev_dta( grp.controls_pre_3T_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_left,1))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_3T_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_right,1))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_3T_allSurg_right,5)) ) );

fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_allT_allSurg_all,1)) & ~cellfun(@isnan, mri_dev_dta( grp.controls_pre_allT_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_left,1))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_allT_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_right,1))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_allT_allSurg_right,5)) ) );

% Pre VPA & DTI
fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_3T_allSurg_all,2)) & ~cellfun(@isnan, dti_dev_dta( grp.controls_pre_3T_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_left,2))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_3T_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_right,2))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_3T_allSurg_right,5)) ) );

fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_allT_allSurg_all,2)) & ~cellfun(@isnan, dti_dev_dta( grp.controls_pre_allT_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_left,2))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_allT_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_right,2))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_pre_allT_allSurg_right,5)) ) );

% Pre VPA & MRI
fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_3T_allSurg_all,2)) & ~cellfun(@isnan, mri_dev_dta( grp.controls_pre_3T_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_left,2))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_3T_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_3T_allSurg_right,2))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_3T_allSurg_right,5)) ) );

fprintf('%d / %d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.controls_pre_allT_allSurg_all,2)) & ~cellfun(@isnan, mri_dev_dta( grp.controls_pre_allT_allSurg_all,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_left,2))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_allT_allSurg_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_pre_allT_allSurg_right,2))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_pre_allT_allSurg_right,5)) ) );

% Post LM2 & DTI
fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_left,7))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_3T_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_right,7))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_3T_ATLonly_right,5)) ) );

fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_left,7))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_allT_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_right,7))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_allT_ATLonly_right,5)) ) );

% Post LM2 & MRI
fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_left,7))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_3T_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_right,7))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_3T_ATLonly_right,5)) ) );

fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_left,7))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_allT_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_right,7))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_allT_ATLonly_right,5)) ) );

% Post VPA & DTI
fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_left,8))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_3T_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_right,8))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_3T_ATLonly_right,5)) ) );

fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_left,8))     & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_allT_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_right,8))    & ~cellfun(@isnan, dti_dev_dta( grp.tle_post_allT_ATLonly_right,5)) ) );

% Post VPA & MRI
fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_left,8))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_3T_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_3T_ATLonly_right,8))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_3T_ATLonly_right,5)) ) );

fprintf('%d / %d \n', ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_left,8))     & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_allT_ATLonly_left,5)) ), ...
    sum( ~cellfun(@isnan, cog_dta( grp.tle_post_allT_ATLonly_right,8))    & ~cellfun(@isnan, mri_dev_dta( grp.tle_post_allT_ATLonly_right,5)) ) );


grp.controls_pre_3T_allSurg_all
grp.controls_pre_allT_allSurg_all
grp.tle_pre_3T_allSurg_left
grp.tle_pre_allT_allSurg_left
grp.tle_pre_3T_allSurg_right
grp.tle_pre_allT_allSurg_right
grp.tle_post_3T_ATLonly_left
grp.tle_post_allT_ATLonly_left
grp.tle_post_3T_ATLonly_right
grp.tle_post_allT_ATLonly_right