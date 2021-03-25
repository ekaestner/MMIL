%%
clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_final_sample_2.csv';
    sbj_nme = mmil_readtext( [prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' sbj_nme ]);

cov_nme = 'Epilepsy_Aging_Final_Sample_ASC_edit.csv';


%% Load ROIs
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.sbj_nme = sbj_nme;

fcfg.dta_nme = [ 'MRI' '_' 'thickness' ];
fcfg.prc_nme = '';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'ROI'];

ejk_project_pull_roi(fcfg)

