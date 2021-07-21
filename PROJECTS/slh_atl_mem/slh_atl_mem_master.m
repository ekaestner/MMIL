clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'slh_atl_mem';

%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_06_22.csv';

cog_tst_nme = { '' '' ...
                '' '' };
            
roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

%% Load Data
slh_atl_mem_master