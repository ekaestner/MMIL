clear; clc;

%%
% Folders
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'PostOperative/';

red_fle = 'sbj000_total_2019_03_27.csv';

grp_fle = 'Language_Impairment_Connectome_SubjectList.csv'; % 'allsubjects_agnostic.csv';

% Subject Names
grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};

% grp_fle_con = grp_fle(strcmpi(grp_fle(:,2),'HC'),:);
    
%% Get Redcap Data
cfg = [];
cfg.prj_dir = prj_dir;
cfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);

% Get Post-Op Categorizations
cfg = [];
cfg.out_dir = [prj_dir '/' 'OUTPUT' '/' prj_nme];
[ pst_cog_dta , ~ ] = ejk_post_cognitive(cfg,sbj_cog,sbj_scn);

% Choose Subjects with Post-Operative Data
pst_cog_fld_nme = fieldnames(pst_cog_dta); pst_cog_fld_nme(1) = [];
sbj_nme = {};
for iPC = 1:numel(pst_cog_fld_nme)
    sbj_nme = [ sbj_nme{:} pst_cog_dta.sbj_nme( ~isnan(pst_cog_dta.(pst_cog_fld_nme{iPC})) ,1)' ];
end

%% Get Redcap Data
cfg = [];
cfg.prj_dir = prj_dir;
cfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);

% Get Post-Op Categorizations
cfg = [];
cfg.out_dir = [prj_dir '/' 'OUTPUT' '/' prj_nme];
[ pst_cog_dta , ~ ] = ejk_post_cognitive(cfg,sbj_cog,sbj_scn);

% Choose Subjects with Post-Operative Data
pst_cog_fld_nme = fieldnames(pst_cog_dta); pst_cog_fld_nme(1) = [];
sbj_nme = {};
for iPC = 1:numel(pst_cog_fld_nme)
    sbj_nme = [ sbj_nme{:} pst_cog_dta.sbj_nme( ~isnan(pst_cog_dta.(pst_cog_fld_nme{iPC})) ,1)' ];
end
grp_fle = [grp_fle(1,:) ; grp_fle(ismember(grp_fle(:,1),unique(sbj_nme)),:)];

% Get laterality
for iGF = 1:size(grp_fle,1)
    grp_fle{iGF,5} = sbj_sze.sbj_sde_ons{strcmpi(sbj_sze.sbj_nme(:,1),grp_fle{iGF,1}),1}; % ,3
end
grp_fle{1,5} = 'OnsetLaterality'; % ,3

% Add in Controls
grp_fle = [ grp_fle ; [grp_fle_con repmat(grp_fle_con(:,2),1,1)] ]; %[ grp_fle ; [ grp_fle_con  repmat(grp_fle_con(:,2),1,4)] ];

% Reload Redcap Data
cfg = [];
cfg.prj_dir = prj_dir;
cfg.red_fle = red_fle;
cfg.sbj_nme = grp_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);

% Get Post-Op Categorizations
cfg = [];
cfg.out_dir = [prj_dir '/' 'OUTPUT' '/' prj_nme];
cfg.neg_oly = 1;
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(cfg,sbj_cog,sbj_scn);

grp_fle = [ grp_fle [ pst_cog_cat(1:end-size(grp_fle_con,1),[8 9 10]) ; repmat(grp_fle_con(:,2),1,3)] ]; % Language
grp_fle = [ grp_fle [ pst_cog_cat(1:end-size(grp_fle_con,1),[2 3 4 5 6 7]) ; repmat(grp_fle_con(:,2),1,6)] ];

% Make Clinical Continuous Data Structure
sbj_cln_cnt.sbj_nme = sbj_dem.sbj_nme;
sbj_cln_cnt.sbj_age = sbj_dem.sbj_age;
sbj_cln_cnt.sbj_edu = sbj_dem.sbj_edu;

sbj_cln_cnt.sbj_age_ons = sbj_sze.sbj_age_ons;
sbj_cln_cnt.sbj_sze_frq = sbj_sze.sbj_sze_frq;
sbj_cln_cnt.sbj_aed_num = sbj_sze.sbj_aed_num;

%
grp_fle = [ grp_fle(2:end,1) repmat({'EPD'},27,1) grp_fle(2:end,2:end) ];

%% Get Neuroimaging Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);

fcfg.dta_typ = 'fib_tfa'; % TRACT FA
tfa_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'wmp_wmd_des'; % TRACT FA
wmd_des_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'gry_thk_des'; % TRACT FA
gry_thk_des_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'vol_dta'; % Alicia fMRI Number of Voxels
vol_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'fmr_lng_ali_num_vox'; % Alicia fMRI Number of Voxels
fmr_alc_dta = ejk_load_mcd_data(fcfg);

%% Data Structures
sbj_cog_use.log_mem_nor_scr_one = sbj_cog.log_mem_nor_scr_one;
sbj_cog_use.log_mem_nor_scr_two = sbj_cog.log_mem_nor_scr_two;
sbj_cog_use.cvl_lfr_nor_scr_pst = sbj_cog.cvl_lfr_nor_scr_pst;
sbj_cog_use.vp1_nor_scr         = sbj_cog.vp1_nor_scr;
sbj_cog_use.vp2_nor_scr         = sbj_cog.vp2_nor_scr;
sbj_cog_use.bnt_nor_scr         = sbj_cog.bnt_nor_scr;
sbj_cog_use.ant_mem_raw_scr     = sbj_cog.ant_mem_raw_scr;
sbj_cog_use.cat_flu_nor_scr     = sbj_cog.cat_flu_nor_scr;
sbj_cog_use.cvl_tot_nor_scr     = sbj_cog.cvl_tot_nor_scr;
sbj_cog_use.ltr_tot_nor_scr     = sbj_cog.ltr_tot_nor_scr;
sbj_cog_use.swt_cor_nor_scr     = sbj_cog.swt_cor_nor_scr;
sbj_cog_use.swt_acc_nor_scr     = sbj_cog.swt_acc_nor_scr;

%% Make Neuroimaging Scatter Plots
% Pre-Operative
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'wmd_dta_pre' 'vol_dta_pre' 'sbj_cln_cnt_pre' 'tfa_dta_pre' 'fmr_alc_dta_pre' };
fcfg.xdt     = { sbj_cog_use   sbj_cog_use   sbj_cog_use       sbj_cog_use   sbj_cog_use   }; % Secondary Folders
fcfg.ydt     = { wmd_des_dta   vol_dta       sbj_cln_cnt       tfa_dta       fmr_alc_dta   }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1            4 ];
fcfg.grp_nme = { { 'EPD' }    { 'L'    'R' } };
fcfg.grp_clr = { { 'orange' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)

% Post-Operative
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'wmd_des_dta' 'vol_dta'   'sbj_cln_cnt' 'tfa_dta'   'fmr_alc_dta' 'gry_thk_des_dta' };
fcfg.xdt     = { pst_cog_dta   pst_cog_dta pst_cog_dta   pst_cog_dta pst_cog_dta   pst_cog_dta       }; % Secondary Folders
fcfg.ydt     = { wmd_des_dta   vol_dta     sbj_cln_cnt   tfa_dta     fmr_alc_dta   gry_thk_des_dta   }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1            4 ];
fcfg.grp_nme = { { 'EPD' }    { 'L'    'R' } };
fcfg.grp_clr = { { 'orange' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)

%% Make Neuroimaging Bar Plots
% Language 
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'wmd_des_dta' 'vol_dta'   'tfa_dta'   'fmr_alc_dta' 'gry_thk_des_dta' }; % 'sbj_cln_cnt'
fcfg.ydt     = { wmd_des_dta   vol_dta     tfa_dta     fmr_alc_dta   gry_thk_des_dta   }; % sbj_cln_cnt

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    3                              4                             5 ];
fcfg.nme_col = grp_fle(1,fcfg.grp_col+1);
fcfg.grp_nme = { { 'HC'    'EPD'}   { 'HC' 'Impaired' 'NoChange' } { 'HC' 'Impaired' 'NoChange'} { 'HC' 'Impaired' 'NoChange'} };
fcfg.grp_clr = { { 'black' 'orange' } { 'black' 'red' 'blue' }     { 'black' 'red' 'blue' }      { 'black' 'red' 'blue' } };
fcfg.xdt     = { [ 1       2 ]        [ 1       3     4 ]          [ 1       3     4 ]           [ 1       3     4 ] };

fcfg.plt_cmp = { { [1 2] }            {[1 2] [2 3] [1 3]}          {[1 2] [2 3] [1 3]}           {[1 2] [2 3] [1 3]} };
fcfg.plt_anv = { []                   [ 1 2 3 ]                    [ 1 2 3 ]                     [ 1 2 3 ] };

ejk_roi_bar(fcfg)

% Memory
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'wmd_des_dta' 'vol_dta'   'tfa_dta'   'fmr_alc_dta' 'gry_thk_des_dta' }; % 'sbj_cln_cnt'
fcfg.ydt     = { wmd_des_dta   vol_dta     tfa_dta     fmr_alc_dta   gry_thk_des_dta   }; % sbj_cln_cnt

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    6                            7                             8                            9 10 11];
fcfg.nme_col = grp_fle(1,fcfg.grp_col+1);
fcfg.grp_nme = { { 'HC'    'EPD'}   { 'HC' 'Impaired' 'NoChange' } { 'HC' 'Impaired' 'NoChange'} { 'HC' 'Impaired' 'NoChange'} };
fcfg.grp_clr = { { 'black' 'orange' } { 'black' 'red' 'blue' }     { 'black' 'red' 'blue' }      { 'black' 'red' 'blue' } };
fcfg.xdt     = { [ 1       2 ]        [ 1       3     4 ]          [ 1       3     4 ]           [ 1       3     4 ] };

fcfg.plt_cmp = { { [1 2] }            {[1 2] [2 3] [1 3]}          {[1 2] [2 3] [1 3]}           {[1 2] [2 3] [1 3]} };
fcfg.plt_anv = { []                   [ 1 2 3 ]                    [ 1 2 3 ]                     [ 1 2 3 ] };

ejk_roi_bar(fcfg)

% Cognitive Control



%% Make Neuroimaging Group Surface Plots
fcfg = [];

fcfg.prj_dir = '/home/ekaestne/PROJECTS/';
fcfg.prj_nme = 'Tests';

fcfg.mes_typ     = 'wmparc_md'; 

fcfg.sbj_grp     = grp_fle;
fcfg.sbj_grp_col = grp_fle(1,2:end);

fcfg.sbj_grp_ind = { [1 2]   [1 2 3]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4]             [2 3 4] };
fcfg.sbj_grp_cmp = { {[1 2]} {[1 2] [1 3] [2 3]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} {[2 3] [2 4] [3 4]} };

fcfg.hms     = {'lhs' 'rhs'};    

ejk_surf_group_avg(fcfg)

%% Make Neuroimaging Correlation Surface Plots
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.sbj_grp     = grp_fle;
fcfg.sbj_grp_col = grp_fle(1,2:end);
fcfg.sbj_grp_ind = { [1] [2 3] [3 4] [3 4] [3 4] [3 4] [3 4] [3 4] [3 4] [3 4] [3 4] };

fcfg.mes_typ     = 'wmparc_md'; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

fcfg.cov         = { pst_cog_dta sbj_sze sbj_dem }; % { sbj_cog      sbj_cog           sbj_cog        sbj_cog               sbj_cog };

fcfg.hms     = {'lhs' 'rhs'};

ejk_surf_group_cor(fcfg)

%%









