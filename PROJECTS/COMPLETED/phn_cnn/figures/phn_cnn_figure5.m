clear; clc;

%% Leave2Out
dta_loc = '/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/Leave2Out/';

% Clinical Variables

% Tracts
trc_nme = '005_TRACTS_Outcomes_all3T_temporal_toall_performance.csv';
trc_dta = cell2mat(mmil_readtext([ dta_loc '/' trc_nme ]));

% Connectome
cnn_nme = '004_CONNECTOME_Outcomes_all3T_temporal_toall_performance.csv';
cnn_dta = cell2mat(mmil_readtext([ dta_loc '/' cnn_nme ]));

%% ACC Random Distribution
figure();
hold on;

hist(trc_dta(:,1),5);

hist(cnn_dta(:,1),5);

%% AUC
figure();
hold on;

hist(trc_dta(:,6),50,rgb('maroon'));

hist(cnn_dta(:,6),50);





