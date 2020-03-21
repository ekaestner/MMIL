clear; clc;

%%
% Folders
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'PostOperative/';

red_fle = 'sbj000_total_2019_03_27.csv';

grp_fle = 'allsubjects_agnostic.csv';

% Subject Names
grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};

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
grp_fle = grp_fle(ismember(grp_fle(:,1),unique(sbj_nme)),:);

% Get laterality
for iGF = 1:size(grp_fle,1)
    grp_fle{iGF,3} = sbj_sze.sbj_sde_ons{strcmpi(sbj_sze.sbj_nme(:,1),grp_fle{iGF,1}),1};
end

% Reload Redcap Data
cfg = [];
cfg.prj_dir = prj_dir;
cfg.red_fle = red_fle;
cfg.sbj_nme = grp_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);

% Get Post-Op Categorizations
cfg = [];
cfg.out_dir = [prj_dir '/' 'OUTPUT' '/' prj_nme];
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(cfg,sbj_cog,sbj_scn);

% Make Clinical Continuous Data Structure
sbj_cln_cnt.sbj_nme = sbj_dem.sbj_nme;
sbj_cln_cnt.sbj_age = sbj_dem.sbj_age;
sbj_cln_cnt.sbj_edu = sbj_dem.sbj_edu;

sbj_cln_cnt.sbj_age_ons = sbj_sze.sbj_age_ons;
sbj_cln_cnt.sbj_sze_frq = sbj_sze.sbj_sze_frq;
sbj_cln_cnt.sbj_aed_num = sbj_sze.sbj_aed_num;

%% Get Neuroimaging Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = grp_fle(:,1);

fcfg.dta_typ = 'fib_tfa'; % TRACT FA
tfa_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'wmp_wmd_des'; % TRACT FA
wmd_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'vol_dta'; % Alicia fMRI Number of Voxels
vol_dta = ejk_load_mcd_data(fcfg);

fcfg.dta_typ = 'fmr_lng_ali_num_vox'; % Alicia fMRI Number of Voxels
fmr_alc_dta = ejk_load_mcd_data(fcfg);

%% Data STruc

%% Make Neuroimaging Scatter Plots
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'vol_dta'   'sbj_cln_cnt' 'tfa_dta'   'fmr_alc_dta' };
fcfg.xdt     = { pst_cog_dta pst_cog_dta   pst_cog_dta pst_cog_dta }; % Secondary Folders
fcfg.ydt     = { vol_dta     sbj_cln_cnt   tfa_dta     fmr_alc_dta }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1            2 ];
fcfg.grp_nme = { { 'EPD' }    { 'L'    'R' } };
fcfg.grp_clr = { { 'orange' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)




