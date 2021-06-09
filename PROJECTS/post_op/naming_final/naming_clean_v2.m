cog_sbj_cln = { 'epd_ucsf010'                             'epd_ucsf018'                      'fc078' };
cog_roi_cln = { {'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst'} {'bnt_raw_scr' 'ant_mem_raw_scr' } {'bnt_raw_scr' 'ant_mem_raw_scr' } };

dti_fib_tfa_sbj_cln = { 'epd_ucsf007' 'epd_ucsf018' 'fc041'}; % epd046 & epd041?
dti_fib_tfa_roi_cln = { {'all_roi'}   {'all_roi'}   {'all_roi'} };

dti_wmp_wfa_sbj_cln = { 'epd_ucsf007' 'epd_ucsf018'  'fc041' };
dti_wmp_wfa_roi_cln = { {'all_roi'}   {'all_roi'}    {'all_roi'} };

mri_vol_sbj_cln = { };
mri_vol_roi_cln = { };


%% Cognitive Data
cog_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv']);

fcfg = [];

fcfg.dta     = cog_dta;
fcfg.dta_col = 2:size(cog_dta,2);
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

fcfg.sbj_nme = dti_fib_tfa_sbj_cln;
fcfg.roi_nme = dti_fib_tfa_roi_cln;

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

fcfg.sbj_nme = dti_wmp_wfa_sbj_cln;
fcfg.roi_nme = dti_wmp_wfa_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'wmparc_FA_wm_aparc_annot' '/'];
fcfg.out_pre_fix = 'wmparc_FA_wm_aparc_annot';
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

%% MRI Volumes - No ICV-cor
mri_vol_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol' '.csv']);

fcfg = [];

fcfg.dta     = mri_vol_dta;
fcfg.dta_col = 5:size(mri_vol_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = mri_vol_sbj_cln;
fcfg.roi_nme = mri_vol_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'subcort_vol' '/'];
fcfg.out_pre_fix = 'subcort_vol_ICV_cor';
ejk_qc_roi(fcfg)
