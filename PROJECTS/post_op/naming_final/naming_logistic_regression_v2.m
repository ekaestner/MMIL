load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor' '/' 'LogisticRegression' '/' ]; ejk_chk_dir( out_dir );

cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];
cln_dta_nme = [ dta_dir '/' 'Clinical'                          '.csv'];
mri_dta_nme = [ dta_dir '/' 'subcort_vol_ICV_cor'      '_' 'QC' '.csv'];
fib_dta_nme = [ dta_dir '/' 'fiber_FA'                 '_' 'QC' '.csv'];
wmp_dta_nme = [ dta_dir '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'];

cat_cut = -1.5;

%% Load Data
cog_dta = mmil_readtext(cog_dta_nme);
cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

cln_dta = mmil_readtext(cln_dta_nme);
cln_dta_col = ejk_fix_column_names(cln_dta(1,2:end));
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

mri_dta = mmil_readtext(mri_dta_nme);
mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
mri_dta_sbj = mri_dta(2:end,1);
mri_dta     = mri_dta(2:end,5:end);

fib_dta = mmil_readtext(fib_dta_nme);
fib_dta_col = ejk_fix_column_names(fib_dta(1,5:end));
fib_dta_sbj = fib_dta(2:end,1);
fib_dta     = fib_dta(2:end,5:end);

wmp_dta = mmil_readtext(wmp_dta_nme);
wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
wmp_dta_sbj = wmp_dta(2:end,1);
wmp_dta     = wmp_dta(2:end,5:end);

prd_dta_sbj = [ cog_dta_sbj ];

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta ];

%% Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = cog_tst_inc
    sbj_nme_hld = [ sbj_nme_hld ; sbj_cog.sbj_nme(~isnan(sbj_cog.(cog_tst_nme{iT}))) ];
end
sbj_nme = unique(sbj_nme_hld);

clear sbj_cog

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

%
fcfg = [];
fcfg.rci = 1;
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

%% Check/Setup Groups
% BNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat_hld = pst_cog_cat(2:end,8);
[ num2cell(pst_cog_dta.bnt_raw_scr_pst(grp.tle_post_3T_ATLonly_left)) cat_hld(grp.tle_post_3T_ATLonly_left) ] 
tabulate(cat_hld(grp.tle_post_3T_ATLonly_left))

cog_sbj_bnt = find(~isnan(pst_cog_dta.bnt_raw_scr_pst));

% ReCat
cog_cat_bnt = repmat({''},numel(pst_cog_dta.bnt_raw_scr_pst),1);
cog_cat_bnt(pst_cog_dta.bnt_raw_scr_pst <= cat_cut) = {'Impaired'};
cog_cat_bnt(pst_cog_dta.bnt_raw_scr_pst > cat_cut) = {'NoChange'};

% ANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat_hld = pst_cog_cat(2:end,9);
[ num2cell(pst_cog_dta.ant_mem_raw_scr_pst(grp.tle_post_3T_ATLonly_left)) cat_hld(grp.tle_post_3T_ATLonly_left) ] 
tabulate(cat_hld(grp.tle_post_3T_ATLonly_left))

cog_sbj_ant = find(~isnan(pst_cog_dta.ant_mem_raw_scr_pst));

% ReCat
cog_cat_ant = repmat({''},numel(pst_cog_dta.ant_mem_raw_scr_pst),1);
cog_cat_ant(pst_cog_dta.ant_mem_raw_scr_pst <= cat_cut) = {'Impaired'};
cog_cat_ant(pst_cog_dta.ant_mem_raw_scr_pst > cat_cut) = {'NoChange'};

%% L-TLE - BNT
cat_hld = pst_cog_cat(2:end,[1 8 9]);

pre_scr_col = find(strcmpi( prd_dta_col, 'xbnt_raw_scr' ));
lft_hip_col = find(strcmpi( prd_dta_col, 'educ' ));
lft_ilf_col = find(strcmpi( prd_dta_col, 'xL_ILF' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'xrh_fusiform' ));

inc_sbj = intersect( grp.tle_post_3T_ATLonly_left ,cog_sbj_bnt );

% Model 1 - 4 Var
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_bnt(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ pre_scr_col lft_hip_col lft_ilf_col rgh_fus_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'} prd_dta_col( 1, [ pre_scr_col lft_hip_col lft_ilf_col rgh_fus_col ]) ];

lcfg.mdl     = { { 'xbnt_raw_scr' 'educ' } ...
                 { 'xL_ILF' 'xrh_fusiform' } ...
                 { 'xL_ILF' 'xL_IFO' } ...
                 }; %{ 'xbnt_raw_scr' 'xLeft_Hippocampus' 'xL_ILF' 'xrh_fusiform' } };
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' ...
                 'WhiteMatter_tracts' ...
                 }; %'Clinical+WhiteMatter' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 3] }; % 3
lcfg.mld_cmp_col = { { rgb('teal') rgb('burnt orange') rgb('light brown') } }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'BNT' '_' 'Logistic'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Scatter Check - White
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),5)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),5)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),4)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),4)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'right fusiform FA'  };
fcfg.ylb = { 'left ILF' };

fcfg.ttl = ['BNT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/BNT_Logistic'];
fcfg.out_nme = 'scatter_white';

ejk_scatter(fcfg)

% Scatter Check - Clinical
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),3)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),3)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),2)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),2)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'education'  };
fcfg.ylb = { 'preoperative score' };

fcfg.ttl = ['BNT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/BNT_Logistic'];
fcfg.out_nme = 'scatter_clinical';

ejk_scatter(fcfg)

%% L-TLE - ANT
cat_hld = pst_cog_cat(2:end,[1 8 9]);

pre_scr_col = find(strcmpi( prd_dta_col, 'xant_mem_raw_scr' ));
lft_hip_col = find(strcmpi( prd_dta_col, 'educ' ));
lft_ilf_col = find(strcmpi( prd_dta_col, 'xL_ILF' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'xrh_fusiform' ));

inc_sbj = intersect( grp.tle_post_3T_ATLonly_left ,cog_sbj_ant );

% Model 1 - 4 Var
lcfg = [];

lcfg.sbj_grp = [  prd_dta_sbj( inc_sbj,1) cog_cat_ant(inc_sbj,1) ];
lcfg.lbl_ord = { 'Impaired' 'NoChange' };

lcfg.dta     = [ prd_dta_sbj( inc_sbj,1) prd_dta( inc_sbj, [ pre_scr_col lft_hip_col lft_ilf_col rgh_fus_col ] ) ];
lcfg.dta_lbl = [ {'sbj_nme'} prd_dta_col( 1, [ pre_scr_col lft_hip_col lft_ilf_col rgh_fus_col ]) ];

lcfg.mdl     = { { 'xant_mem_raw_scr' 'educ' } ...
                 { 'xL_ILF' 'xrh_fusiform' } ...
                 { 'xL_ILF' 'xL_IFO' } ...
                 }; %{ 'xant_mem_raw_scr' 'xLeft_Hippocampus' 'xL_ILF' 'xrh_fusiform' } };
lcfg.mdl_nme = { 'Clinical' ...
                 'WhiteMatter' ...
                 'WhiteMatter_tracts' ...
                 }; %'Clinical+WhiteMatter' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2 3] }; %
lcfg.mld_cmp_col = { { rgb('teal') rgb('burnt orange') rgb('light brown')} }; %rgb('purple')
lcfg.mdl_cmp_nme = { [ 'ANT' '_' 'Logistic'] };

lcfg.out_dir = out_dir;

mmil_log_reg(lcfg)

% Scatter Check - White
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),5)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),5)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),4)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),4)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'right fusiform FA'  };
fcfg.ylb = { 'left ILF' };

fcfg.ttl = ['ANT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/ANT_Logistic'];
fcfg.out_nme = 'scatter_white';

ejk_scatter(fcfg)

% Scatter Check - Clinical
fcfg = [];
        
fcfg.xdt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),3)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),3)) };
fcfg.ydt     = { cell2mat(lcfg.dta( strcmpi(lcfg.sbj_grp(:,2),'Impaired'),2)) cell2mat(lcfg.dta(strcmpi(lcfg.sbj_grp(:,2),'NoChange'),2)) };

fcfg.trd_lne = [ 0 0 ];

fcfg.fce_col = { rgb('red')   rgb('blue') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlb = { 'education'  };
fcfg.ylb = { 'preoperative score' };

fcfg.ttl = ['ANT'];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor/LogisticRegression/Models/LogReg/ANT_Logistic'];
fcfg.out_nme = 'scatter_clinical';

ejk_scatter(fcfg)

[ lcfg.sbj_grp prd_dta( inc_sbj,5) lcfg.dta ]

%% R-TLE
