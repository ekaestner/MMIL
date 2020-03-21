%%
clear; clc;

sbj_grp = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/tmp/PhenotypeConnectome/Language_Impairment_Connectome_SubjectList.csv');

% Folders
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'PhenotypeConnectome';

%
red_fle = 'sbj000_total_2019_03_27.csv';

cfg = [];
cfg.prj_dir = prj_dir;
cfg.red_fle = red_fle;
cfg.sbj_nme = sbj_grp(:,1);
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);

%% Table 1: Demographics Table
% Two Groups
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;
fcfg.tbl_nme = 'two_group';

fcfg.sbj_nme = sbj_grp;
fcfg.grp_col = 1; % 1; % 2; % 3; % 4;
fcfg.grp_nme = {'languageimpaired' 'nonimpaired'};  % { 'Language Impairment' 'No Language Impairment'  }; % { 'left' 'right' } % {} % { 'HC' 'EPD' }

fcfg.mes_sbj_nme = { sbj_dem.sbj_nme sbj_dem.sbj_nme sbj_dem.sbj_nme     sbj_sze.sbj_nme  sbj_sze.sbj_nme     sbj_sze.sbj_nme     sbj_sze.sbj_nme      sbj_sze.sbj_nme     sbj_cog.sbj_nme      sbj_cog.sbj_nme         sbj_cog.sbj_nme   }; % { sbj_dem.sbj_nme sbj_dem.sbj_nme };
fcfg.mes_nme     = { 'Education'     'Sex: M/F'      'Handedness: L/R/A' 'MTS: Yes/No'    'Side: L/R/B'       'Age of Onset'      'Number of AEDs'     'Seizure Frequency' 'Boston Naming Test' 'Auditory Naming Test'  'D-KEFS Category'           }; % { 'Education'     'Sex: M/F' };
fcfg.mes_lbl     = { sbj_dem.sbj_edu sbj_dem.sbj_sex sbj_dem.sbj_hnd     sbj_sze.sbj_mts  sbj_sze.sbj_sde_ons sbj_sze.sbj_age_ons sbj_sze.sbj_aed_num  sbj_sze.sbj_sze_frq sbj_cog.bnt_nor_scr  sbj_cog.ant_mem_raw_scr sbj_cog.cat_flu_nor_scr  }; % { sbj_dem.sbj_edu sbj_dem.sbj_sex };
fcfg.mes_typ     = { 'con'           'ord'           'ord'               'ord'            'ord'               'con'               'con'                'con'               'con'                'con'                   'con' }; % { 'con'           'ord' };
fcfg.ord_ord     = { {}              { 'M' 'F' }     {'L' 'R' }          {'L' 'R' 'N/A' } {'L' 'R'}           {}                  {}                   {}                  {}                   {}                      {} }; % { {}              { 'M' 'F' } };

ejk_clinical_table(fcfg)

% Four Groups
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;
fcfg.tbl_nme = 'four_group';

fcfg.sbj_nme = sbj_grp;
fcfg.grp_col = 3; % 1; % 2; % 3; % 4;
fcfg.grp_nme = { 'languageimpaired_UCSD' 'nonimpaired_UCSD' 'languageimpaired_UCSF' 'nonimpaired_UCSF' };  % { 'Language Impairment' 'No Language Impairment'  }; % { 'left' 'right' } % {} % { 'HC' 'EPD' }

fcfg.mes_sbj_nme = { sbj_dem.sbj_nme sbj_dem.sbj_nme sbj_dem.sbj_nme     sbj_sze.sbj_nme  sbj_sze.sbj_nme     sbj_sze.sbj_nme     sbj_sze.sbj_nme      sbj_sze.sbj_nme     sbj_cog.sbj_nme      sbj_cog.sbj_nme         sbj_cog.sbj_nme   }; % { sbj_dem.sbj_nme sbj_dem.sbj_nme };
fcfg.mes_nme     = { 'Education'     'Sex: M/F'      'Handedness: L/R/A' 'MTS: Yes/No'    'Side: L/R/B'       'Age of Onset'      'Number of AEDs'     'Seizure Frequency' 'Boston Naming Test' 'Auditory Naming Test'  'D-KEFS Category'           }; % { 'Education'     'Sex: M/F' };
fcfg.mes_lbl     = { sbj_dem.sbj_edu sbj_dem.sbj_sex sbj_dem.sbj_hnd     sbj_sze.sbj_mts  sbj_sze.sbj_sde_ons sbj_sze.sbj_age_ons sbj_sze.sbj_aed_num  sbj_sze.sbj_sze_frq sbj_cog.bnt_nor_scr  sbj_cog.ant_mem_raw_scr sbj_cog.cat_flu_nor_scr  }; % { sbj_dem.sbj_edu sbj_dem.sbj_sex };
fcfg.mes_typ     = { 'con'           'ord'           'ord'               'ord'            'ord'               'con'               'con'                'con'               'con'                'con'                   'con' }; % { 'con'           'ord' };
fcfg.ord_ord     = { {}              { 'M' 'F' }     {'L' 'R' }          {'L' 'R' 'N/A' } {'L' 'R'}           {}                  {}                   {}                  {}                   {}                      {} }; % { {}              { 'M' 'F' } };

ejk_clinical_table(fcfg)

%% Table 2: Overall Performance Table

%% Table 3: Leave 2 out Performance Table
phn_cnn_tbl2

%% Figure 1: 

%% Figure 2: Model Framework

%% Figure 3: ROC Curves
phn_cnn_figure3

%% Figure 4: Feature Importances

%% Figure 5: Distributions

%% Figure 6: Connectome Plot





