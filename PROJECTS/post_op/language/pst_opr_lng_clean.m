cog_sbj_cln = { 'epd_ucsf010' };
cog_roi_cln = { {'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst'} };

dti_fib_sbj_cln = { 'epd_ucsf007' };
dti_fib_roi_cln = { {'all_roi'} };

dti_wmp_sbj_cln = { };
dti_wmp_roi_cln = { };

mri_thk_sbj_cln = { };
mri_thk_roi_cln = { };

mri_vol_sbj_cln = { };
mri_vol_roi_cln = { };

rsf_ctx_sbj_cln = { 'epd_ucsf003' };
rsf_ctx_roi_cln = { {'all_roi'} };

rsf_vol_sbj_cln = { 'epd_ucsf003' };
rsf_vol_roi_cln = { {'all_roi'} };

%% Cognitive Data
cog_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv']);

fcfg = [];

fcfg.dta     = cog_dta;
fcfg.dta_col = 2:9;
fcfg.sbj_col = 1;

fcfg.sbj_nme = cog_sbj_cln;
fcfg.roi_nme = cog_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,2:end));
fcfg.dta_lbl = cln_dta(1,2:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'Cognitive' '/'];
fcfg.out_pre_fix = 'Cognitive';
ejk_qc_roi(fcfg)

%% DTI fiber FA
dti_fib_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '.csv']);

fcfg = [];

fcfg.dta     = dti_fib_dta;
fcfg.dta_col = 5:46;
fcfg.sbj_col = 1;

fcfg.sbj_nme = dti_fib_sbj_cln;
fcfg.roi_nme = dti_fib_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'fiber_FA' '/'];
fcfg.out_pre_fix = 'fiber_FA';
ejk_qc_roi(fcfg)

%% DTI wmparc WM
dti_wmp_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot' '.csv']);

fcfg = [];

fcfg.dta     = dti_wmp_dta;
fcfg.dta_col = 5:size(dti_wmp_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = dti_wmp_sbj_cln;
fcfg.roi_nme = dti_wmp_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'wmparc_FA_wm_aparc_annot' '/'];
fcfg.out_pre_fix = 'wmparc_FA_wm_aparc_annot';
ejk_qc_roi(fcfg)

%% MRI Cortical Thickness
mri_thk_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'cort_thick_ctx_aparc_annot' '.csv']);

fcfg = [];

fcfg.dta     = mri_thk_dta;
fcfg.dta_col = 5:size(mri_thk_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = mri_thk_sbj_cln;
fcfg.roi_nme = mri_thk_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'cort_thick_ctx_aparc_annot' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'cort_thick_ctx_aparc_annot' '/'];
fcfg.out_pre_fix = 'cort_thick_ctx_aparc_annot';
ejk_qc_roi(fcfg)

%% MRI Volumes
mri_vol_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor' '.csv']);

fcfg = [];

fcfg.dta     = mri_vol_dta;
fcfg.dta_col = 5:size(mri_vol_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = mri_vol_sbj_cln;
fcfg.roi_nme = mri_vol_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'subcort_vol_ICV_cor' '/'];
fcfg.out_pre_fix = 'subcort_vol_ICV_cor';
ejk_qc_roi(fcfg)

%% Resting State - Cortex
rsf_ctx_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'var_ctx_aparc_annot' '.csv']);

fcfg = [];

fcfg.dta     = rsf_ctx_dta;
fcfg.dta_col = 5:size(rsf_ctx_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = rsf_ctx_sbj_cln;
fcfg.roi_nme = rsf_ctx_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'var_ctx_aparc_annot' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'var_ctx_aparc_annot' '/'];
fcfg.out_pre_fix = 'var_ctx_aparc_annot';
ejk_qc_roi(fcfg)

%% Resting State - Volume
rsf_vol_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'var_vol' '.csv']);

fcfg = [];

fcfg.dta     = rsf_vol_dta;
fcfg.dta_col = 5:size(rsf_vol_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = rsf_vol_sbj_cln;
fcfg.roi_nme = rsf_vol_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'var_vol' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'var_vol' '/'];
fcfg.out_pre_fix = 'var_vol';
ejk_qc_roi(fcfg)

