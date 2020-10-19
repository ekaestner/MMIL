clear; clc;

%% Put Together Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_dti_dfa = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Data/cognitive.csv');

% Data
dta_sbj_nme = dta_dti_dfa(2:end,1);
dta_roi_nme = dta_dti_dfa(1,2:end);
dta_dta     = cell2mat(dta_dti_dfa(2:end,2:end));

%% Data Cleaning 
% Check data for missing values
sum(isnan(dta_dta(:)))

%% LPA - Original
% Fit 
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.dta     = dta_dta;
fcfg.dta_nme = dta_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/';

ejk_fit_lpa( fcfg );

% LPA Model
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.mdl_num = 3;

fcfg.dta     = dta_dta;
fcfg.dta_nme = dta_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/';

ejk_lpa( fcfg );

% LPA Model
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.mdl_num = 4;

fcfg.dta     = dta_dta;
fcfg.dta_nme = dta_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/';

ejk_lpa( fcfg );

% LPA Model
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.mdl_num = 4;

fcfg.dta     = dta_dta;
fcfg.dta_nme = dta_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive/';

ejk_lpa( fcfg );

%% Kmeans - Original
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.dta     = dta_dta;
fcfg.dta_nme = dta_roi_nme;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_cognitive_kmeans/';

ejk_kmeans( fcfg );

