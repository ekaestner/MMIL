cog_sbj_cln = { 'epd_ucsf010'                             'epd_ucsf018'                      'epd024'                          'epd004' };
cog_roi_cln = { {'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst'} {'bnt_raw_scr' 'ant_mem_raw_scr' } {'bnt_raw_scr' 'cat_flu_nor_scr'} {'bnt_raw_scr_pst'} };

dti_fib_tfa_sbj_cln = { 'epd_ucsf007' 'epd_ucsf018' };
dti_fib_tfa_roi_cln = { {'all_roi'}   {'all_roi'} };

dti_wmp_wfa_sbj_cln = { 'epd_ucsf007' 'epd_ucsf018'  };
dti_wmp_wfa_roi_cln = { {'all_roi'}   {'all_roi'} };

mri_thk_sbj_cln = { };
mri_thk_roi_cln = { };

mri_vol_sbj_cln = { };
mri_vol_roi_cln = { };

alc_sbj_cln = { 'epd058'      'epd091' 'epd_ucsf009' 'epd_ucsf013' 'epd_ucsf020' 'epd_ucsf033' 'epd_ucsf034' ...
                'epd_ucsf039' 'fc007'  'fc008'       'fc009'       'fc010'       'fc011'       'fc016' ...
                'fc021'       'fc024'  'fc026'       'fc027'       'fc029'       'fc030'       'fc031' ...
                'fc032'       'fc033'  'fc040'       'fc041'       'fc042'       'fc051'       'fc055' ...
                'fc056'       'fc059'  'fc079'       'fc081'       'fc094' ...
                'epd_ucsf011' 'epd_ucsf016' 'fc071' };
alc_roi_cln = { {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} ...
                {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} ...
                {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} ...
                {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} ...
                {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} {'all_roi'} ...
                {'all_roi'} {'all_roi'} {'all_roi'} };

rsf_ctx_sbj_cln = {  };
rsf_ctx_roi_cln = {  };

rsf_vol_sbj_cln = {  };
rsf_vol_roi_cln = {  };

rsf_net_sbj_cln = {  };
rsf_net_roi_cln = {  };

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'cort_thick_ctx_aparc_annot' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

%% Functional MRI - Alicia
alc_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'alicia' '.csv']);

fcfg = [];

fcfg.dta     = alc_dta;
fcfg.dta_col = 5:size(alc_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = alc_sbj_cln;
fcfg.roi_nme = alc_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'alicia' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'alicia' '/'];
fcfg.out_pre_fix = 'alicia';
ejk_qc_roi(fcfg)

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'alicia' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'var_ctx_aparc_annot' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

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

% LI
[ dta_out , dta_lbl ] = ejk_create_laterality_index( cln_dta(2:end,:) , cln_dta(1,:));
cln_dta_lat = [ cln_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'var_vol' '_' 'LI' '_' 'QC' '.csv'], cln_dta_lat)

%% Resting State - Networks
rsf_vol_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'network' '.csv']);

fcfg = [];

fcfg.dta     = rsf_vol_dta;
fcfg.dta_col = 5:size(rsf_vol_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = rsf_vol_sbj_cln;
fcfg.roi_nme = rsf_vol_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'network' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
fcfg.dta_lbl = cln_dta(1,5:end);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'network' '/'];
fcfg.out_pre_fix = 'network';
ejk_qc_roi(fcfg)

