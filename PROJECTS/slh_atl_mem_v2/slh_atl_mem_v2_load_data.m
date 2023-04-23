
% Demographics
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_cph ];
fcfg.dta_col = 2;
[ dem_nme_dta, dem_nme_row, dem_nme_col] = ejk_dta_frm( fcfg );

%% Load Demographics & Clinical
% Initial Redcap Load
fcfg = [];
fcfg.red_fle = [ dta_dir '/' red_cap_fle ];
fcfg.sep     = '|';
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

% Reassign
red_cap_dta.sbj_dem = sbj_dem;
red_cap_dta.sbj_sze = sbj_sze;
red_cap_dta.sbj_scn = sbj_scn;
red_cap_dta.sbj_cog = sbj_cog;
red_cap_dta.sbj_emo = sbj_emo;
red_cap_dta.sbj_srg = sbj_srg;

% Initial Cognitive Emory Load
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' emy_cog_fle ];
fcfg.dta_col = 2;
[ emy_cog_dta, emy_cog_sbj, emy_cog_col] = ejk_dta_frm( fcfg );

% Initial Demographic Emory Load
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' emy_dem_fle ];
fcfg.dta_col = 2;
[ emy_dem_dta, emy_dem_sbj, emy_dem_col] = ejk_dta_frm( fcfg );

%% Combine Data & Recode
slh_atl_mem_v2_combine_emory

%% Load Imaging Data
slh_atl_mem_v2_load_neurobio

%% Sub-select subjects for study


%% Save



