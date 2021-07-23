
%% L-TLE - ANT
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor' '/' 'LogisticRegression' '/' ]; ejk_chk_dir( out_dir );

cat_cut = -1.28;

%% Load Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
% fcfg.all_num = 1;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ mri_dta, mri_dta_sbj, mri_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ fib_dta, fib_dta_sbj, fib_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ wmp_dta, wmp_dta_sbj, wmp_dta_col] = ejk_dta_frm( fcfg );

prd_dta_sbj = [ cog_dta_sbj ];

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta ];

bnt_scr_col = find(strcmpi( prd_dta_col, 'bnt_raw_scr' ));
ant_scr_col = find(strcmpi( prd_dta_col, 'ant_mem_raw_scr' ));
edu_col     = find(strcmpi( prd_dta_col, 'Educ' ));
lft_hip_col = find(strcmpi( prd_dta_col, 'Left_Hippocampus' ));
rgh_hip_col = find(strcmpi( prd_dta_col, 'Right_Hippocampus' ));
lft_ilf_col = find(strcmpi( prd_dta_col, 'L_ILF' ));
rgh_ilf_col = find(strcmpi( prd_dta_col, 'R_ILF' ));
lft_ifo_col = find(strcmpi( prd_dta_col, 'L_IFO' ));
rgh_ifo_col = find(strcmpi( prd_dta_col, 'R_IFO' ));
lft_fus_col = find(strcmpi( prd_dta_col, 'lh_fusiform' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'rh_fusiform' ));

%% Check/Setup Groups
% BNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_sbj_bnt = find(~isnan(cell2mat(cog_dta(:,3))));

% ReCat
cog_cat_bnt = repmat({''},numel(cog_dta(:,3)),1);
cog_cat_bnt( cell2mat(cog_dta(:,3)) <= cat_cut) = {'Impaired'};
cog_cat_bnt( cell2mat(cog_dta(:,3)) > cat_cut) = {'NoChange'};

% ANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_sbj_ant = find(~isnan(cell2mat(cog_dta(:,4))));

% ReCat
cog_cat_ant = repmat({''},numel(cog_dta(:,4)),1);
cog_cat_ant( cell2mat(cog_dta(:,4)) <= cat_cut) = {'Impaired'};
cog_cat_ant( cell2mat(cog_dta(:,4)) > cat_cut) = {'NoChange'};

%% L-TLE - BNT
inc_sbj = intersect( grp.tle_post_3T_ATLonly_left ,cog_sbj_bnt );

% Setup Explore %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_bnt(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ bnt_scr_col edu_col lft_ilf_col rgh_fus_col rgh_ifo_col lft_hip_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'}             prd_dta_col( 1,   [ bnt_scr_col edu_col lft_ilf_col rgh_fus_col rgh_ifo_col lft_hip_col ]) ];

lcfg.mdl     = { { 'bnt_raw_scr' 'Educ' }  ...
                 { 'L_ILF' 'rh_fusiform' } ...
                 { 'L_ILF' 'R_IFO' }       ...
                 { 'Left_Hippocampus' 'rh_fusiform' } }; 
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' ...
                 'WhiteMatter_tracts' ...
                 'Fusiform_Hippocampus' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 3 4] }; % 3
lcfg.mld_cmp_col = { { rgb('teal') rgb('burnt orange') rgb('light brown')  rgb('gold')} }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'BNT' '_' 'Logistic' '_' 'Explore_inc'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Setup Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_bnt(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ bnt_scr_col edu_col lft_ilf_col rgh_fus_col rgh_ifo_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'}             prd_dta_col( 1,   [ bnt_scr_col edu_col lft_ilf_col rgh_fus_col rgh_ifo_col ]) ];

lcfg.mdl     = { { 'bnt_raw_scr' 'Educ' }  ...
                 { 'L_ILF' 'rh_fusiform' } }; 
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 ] }; % 3
lcfg.mld_cmp_col = { { rgb('teal') rgb('burnt orange') } }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'BNT' '_' 'Logistic' '_' 'Current_inc'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Scatter Check - White %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),5)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),5)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),4)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),4)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'right fusiform FA'  };
fcfg.ylb = { 'left ILF' };

fcfg.ttl = ['BNT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/BNT_Logistic_Current'];
fcfg.out_nme = 'scatter_white';

ejk_scatter(fcfg)

% Scatter Check - Clinical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),3)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),3)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),2)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),2)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'education'  };
fcfg.ylb = { 'preoperative score' };

fcfg.ttl = ['BNT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/BNT_Logistic_Current'];
fcfg.out_nme = 'scatter_clinical';

ejk_scatter(fcfg)

%% L-TLE - ANT
inc_sbj = intersect( grp.tle_post_3T_ATLonly_left ,cog_sbj_ant );

% Setup Explore %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_ant(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ ant_scr_col edu_col lft_ilf_col rgh_fus_col lft_ifo_col lft_hip_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'}             prd_dta_col( 1,   [ ant_scr_col edu_col lft_ilf_col rgh_fus_col lft_ifo_col lft_hip_col ]) ];

lcfg.mdl     = { { 'ant_mem_raw_scr' 'Educ' }  ...
                 { 'L_ILF' 'rh_fusiform' } ...
                 { 'L_ILF' 'L_IFO' }       ...
                 { 'Left_Hippocampus' 'rh_fusiform' } }; 
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' ...
                 'WhiteMatter_tracts' ...
                 'Fusiform_Hippocampus' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 3 4] }; % 3
lcfg.mld_cmp_col = { { rgb('teal') rgb('burnt orange') rgb('light brown')  rgb('gold')} }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'ANT' '_' 'Logistic' '_' 'Explore_inc'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Setup Current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_ant(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ ant_scr_col edu_col lft_ilf_col lft_ifo_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'}             prd_dta_col( 1,   [ ant_scr_col edu_col lft_ilf_col lft_ifo_col ]) ];

lcfg.mdl     = { { 'ant_mem_raw_scr' 'Educ' }  ...
                 { 'L_ILF'           'L_IFO' } }; 
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 ] }; % 3
lcfg.mld_cmp_col = { { rgb('teal') rgb('light brown') } }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'ANT' '_' 'Logistic' '_' 'Current_inc'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Scatter Check - White %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),5)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),5)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),4)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),4)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'left IFOF'  };
fcfg.ylb = { 'left ILF' };

fcfg.ttl = ['ANT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/ANT_Logistic_Current'];
fcfg.out_nme = 'scatter_white';

ejk_scatter(fcfg)

% Scatter Check - Clinical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),3)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),3)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),2)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),2)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'education'  };
fcfg.ylb = { 'preoperative score' };

fcfg.ttl = ['ANT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/ANT_Logistic_Current'];
fcfg.out_nme = 'scatter_clinical';

ejk_scatter(fcfg)

%% R-TLE
