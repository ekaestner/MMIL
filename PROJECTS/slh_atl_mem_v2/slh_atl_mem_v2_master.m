% addpath(genpath('/home/ekaestner/gitrep/MMIL/')); rmpath(genpath('/home/ekaestner/gitrep/MMIL/EXTERNAL/spm12/external/fieldtrip/compat/matlablt2010b')); rmpath(genpath('/home/ekaestner/gitrep/MMIL/EXTERNAL/fieldtrip-20201023/compat/matlablt2010b/'))

clear; clc;

slh_atl_mem_v2_constants

%% Load Data
slh_atl_mem_v2_load_data

%% Sample & Groups
slh_atl_mem_v2_create_groups

slh_atl_mem_v2_QC

%% Demographics/Clinical
slh_atl_mem_v2_table1

%% Cognitive
slh_atl_mem_v2_table2

%% Correlations
slh_atl_mem_v2_table3

slh_atl_mem_v2_fig2

slh_atl_mem_v2_table4

%% Investigate
slh_atl_mem_v2_nuisance