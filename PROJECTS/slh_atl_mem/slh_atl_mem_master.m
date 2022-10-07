clear; clc;

slh_atl_mem_constants

%% Load Data, QC, Define Groups
slh_atl_mem_load_data

slh_atl_mem_define_groups_v2

slh_atl_mem_qal_con

slh_atl_mem_load_nuisance

%% Initial Analysis
% Clinical Stats & Tables
slh_atl_mem_cln

% Cognitive Tables
slh_atl_mem_cog

% Initial Neuroimaging exploration
slh_atl_mem_neu_cog_cor % Initial Cognitive/Neuroimaging correlations & tables

%% Second pass
slh_atl_mem_cor_v2

%% Figures & Tables
slh_atl_mem_fig1
slh_atl_mem_fig2
slh_atl_LI_surf

%% Misc / Scratch
slh_atl_mem_cog_inv_ini % Initial Cognitive exploration

slh_atl_mem_robust_regression

slh_atl_mem_load_nuisance