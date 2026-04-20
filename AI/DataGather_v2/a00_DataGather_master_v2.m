addpath(genpath('/home/ekaestne/gitrep/MMIL/'));

clear; clc;

a01_1_DataGather_constants_v2
a01_2_DataGather_setup_covariates_v2

a01_3_DataGather_investigate_Redcap_v3

%% Data Gathering: Imaging
a02_DataGather_characterize_Imaging_v2

%% Data Gathering: Covariates
a03_DataGather_characterize_Covariates_v2

a04_DataGather_tbl_ini_v2

%% Create Groups
a05_1_DataGather_grp_v2

a05_2_PicklePrep_CMIG

a05_DataGather_text_snippet_prep_v2

%% Data Gathering: Atlases
a06_DataGather_atlas

%% Scan Info: QC & Biovariates
a06_biovariate_gather

%% Data Checks

