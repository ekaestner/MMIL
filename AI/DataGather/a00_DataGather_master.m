addpath(genpath('/home/ekaestne/gitrep/MMIL/'));

clear; clc;

a01_DataGather_constants
a01_DataGather_setup_covariates

a01_DataGather_investigate_Redcap_v2

%% Data Gathering: Imaging
a02_DataGather_characterize_BIDS

a03_DataGather_characterize_CAT12

%% Data Gathering: Covariates
a04_DataGather_characterize_Redcap

a04_DataGather_tbl_ini

%% Create Groups
a05_DataGather_grp

a05_DataGather_pickle_prep

a05_DataGather_text_snippet_prep

%% Data Gathering: Atlases
a06_DataGather_atlas

%% Scan Info: QC & Parameters


%% Data Checks

