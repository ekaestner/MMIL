clear; clc;

dta_loc     = '/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/data/';
dta_prf_loc = '/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/';
dta_l2o_loc = '/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/Leave2Out/';

%% Clinical Variables
trc_dta_nme = '003_CLINICAL_Outcomes_all3T_importances_labels.csv';
trc_dta     = mmil_readtext([ dta_prf_loc '/' trc_dta_nme ]);
trc_dta_lbl = trc_dta';

trc_prf_nme = '003_CLINICAL_Outcomes_all3T_importances.csv';
trc_prf_dta = cell2mat(mmil_readtext([ dta_prf_loc '/' trc_prf_nme ]));

trc_l2o_nme = '006_CLINICAL_Outcomes_all3T_temporal_toall_importances.csv';
trc_l2o_dta = cell2mat(mmil_readtext([ dta_l2o_loc '/' trc_l2o_nme ]));

[ trc_dta(1,2:end)' num2cell(trc_prf_dta) num2cell(mean(trc_l2o_dta,2)) ];

% Put together
trc_prf_dta = roundsd(trc_prf_dta,2);

[ ~ , trc_ind ] = sort( trc_prf_dta );

hld_lat = trc_prf_dta(trc_ind(7));

trc_prf_dta = trc_prf_dta(trc_ind([6 8:15]));
trc_dta_lbl = trc_dta_lbl(trc_ind([6 8:15]));
trc_ind     = 1:10;

% Plot
fcfg = [];

for iTR = 1:9
    fcfg.xdt{iTR}     = trc_prf_dta(trc_ind(iTR));
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('greyish blue');
    fcfg.ylb{iTR}     = trc_dta_lbl{trc_ind(iTR)};
end

sex_hld = fcfg.xdt{5};

fcfg.xdt{5} = fcfg.xdt{4} + hld_lat;
fcfg.xdt{4} = sex_hld;

fcfg.ylb{1} = 'Handedness';
fcfg.ylb{2} = 'MTS Status';
fcfg.ylb{4} = 'Sex';
fcfg.ylb{5} = 'Seizure Laterality';

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure4/';
fcfg.out_nme = 'ClinicalImportance_v2';

ejk_bar_vert(fcfg)

%% Tracts
trc_dta_nme = 'language_tracts.csv';
trc_dta     = mmil_readtext([ dta_loc '/' trc_dta_nme ]);
trc_dta_lbl = trc_dta(1,2:end);

trc_prf_nme = '002_TRACTS_Outcomes_all3T_importances.csv';
trc_prf_dta = cell2mat(mmil_readtext([ dta_prf_loc '/' trc_prf_nme ]));

trc_l2o_nme = '005_TRACTS_Outcomes_all3T_temporal_toall_importances.csv';
trc_l2o_dta = cell2mat(mmil_readtext([ dta_l2o_loc '/' trc_l2o_nme ]));

[ trc_dta(1,2:end)' num2cell(trc_prf_dta) num2cell(mean(trc_l2o_dta,2)) ];

% Put together
trc_prf_dta = roundsd(trc_prf_dta,2);

[ ~ , trc_ind ] = sort( trc_prf_dta );

% Plot
fcfg = [];

for iTR = 1:numel(trc_ind)
    fcfg.xdt{iTR}     = trc_prf_dta(trc_ind(iTR));
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('dark maroon');
    fcfg.ylb{iTR}     = trc_dta_lbl{trc_ind(iTR)}(1:end-3);
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure4/';
fcfg.out_nme = 'TractImportance';

ejk_bar_vert(fcfg)

%% Connectome
cnn_dta_lbl = strcat( repmat({'PC-'},1,40) , cellfun(@num2str,num2cell(1:40),'uni',0) );

cnn_prf_nme = '001_CONNECTOME_Outcomes_all3T_temporal_toall_importances.csv';
cnn_prf_dta = cell2mat(mmil_readtext([ dta_prf_loc '/' cnn_prf_nme ]));

cnn_l2o_nme = '004_CONNECTOME_Outcomes_all3T_temporal_toall_importances.csv';
cnn_l2o_dta = cell2mat(mmil_readtext([ dta_l2o_loc '/' cnn_l2o_nme ]));

dta_hld = roundsd(mean(cnn_l2o_dta,2),3);
dta_hld(dta_hld<.01) = 0;

% Put together
cnn_prf_dta = roundsd(cnn_prf_dta,2);

[ ~ , cnn_ind ] = sort( cnn_prf_dta );

cnn_prf_dta = cnn_prf_dta(cnn_ind(32:40));
cnn_dta_lbl = cnn_dta_lbl(cnn_ind(32:40));
cnn_ind     = 1:9;

% Plot
fcfg = [];

fcfg.xdt{1}     = 0;
fcfg.ydt{1}     = 0;
fcfg.fce_col{1} = rgb('purple');
fcfg.ylb{1}     = 'Other PCs';

for iTR = 1:9
    fcfg.xdt{iTR+1}     = cnn_prf_dta(cnn_ind(iTR));
    fcfg.ydt{iTR+1}     = iTR;
    fcfg.fce_col{iTR+1} = rgb('purple');
    fcfg.ylb{iTR+1}     = cnn_dta_lbl{cnn_ind(iTR)};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure4/';
fcfg.out_nme = 'ConnectomeImportance';

ejk_bar_vert(fcfg)

