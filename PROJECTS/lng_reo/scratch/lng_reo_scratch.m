%%
clear; clc;

%  Overall Data Locations
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'LanguageReorganization';

% File Names
grp_fle = 'Lng_Reo_subject_list.csv';
red_fle = 'sbj000_total_2019_03_27.csv';

%% Load Data1
% Lists %%%%%%%%%%%%%%%%%%%%%
fsr_sbj_nme = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv']);
bld_sbj_nme = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv']);
grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' grp_fle]);

% Redcap Data %%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
fcfg.sbj_nme = grp_fle(:,1);
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);
