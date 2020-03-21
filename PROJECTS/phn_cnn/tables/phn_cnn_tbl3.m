dta_loc = '/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/Leave2Out/';

%% Clinical Variables
trc_nme = '006_CLINICAL_Outcomes_all3T_temporal_toall_performance.csv';

cln_dta = cell2mat(mmil_readtext([ dta_loc '/' trc_nme ]));

acc_men_hld     = [ num2str(roundsd(mean(cln_dta(:,1)),2)) ];
    acc_std_hld = [ num2str(roundsd(std(cln_dta(:,1)),2)) ];
ppv_men_hld     = [ num2str(roundsd(mean(cln_dta(:,2)),2)) ];
    ppv_std_hld = [ num2str(roundsd(std(cln_dta(:,2)),2)) ];
npv_men_hld     = [ num2str(roundsd(mean(cln_dta(:,3)),2)) ];
    npv_std_hld = [ num2str(roundsd(std(cln_dta(:,3)),2)) ];
sns_men_hld     = [ num2str(roundsd(mean(cln_dta(:,4)),2)) ];
    sns_std_hld = [ num2str(roundsd(std(cln_dta(:,4)),2)) ];
spe_men_hld     = [ num2str(roundsd(mean(cln_dta(:,5)),2)) ];
    spe_std_hld = [ num2str(roundsd(std(cln_dta(:,5)),2)) ];
auc_men_hld     = [ num2str(roundsd(mean(cln_dta(:,6)),2)) ];
    auc_std_hld = [ num2str(roundsd(std(cln_dta(:,6)),2)) ];

[ acc_men_hld(2:end) ' +/- ' acc_std_hld(2:end)]
[ ppv_men_hld(2:end) ' +/- ' ppv_std_hld(2:end)]
[ npv_men_hld(2:end) ' +/- ' npv_std_hld(2:end)]
[ sns_men_hld(2:end) ' +/- ' sns_std_hld(2:end)]
[ spe_men_hld(2:end) ' +/- ' spe_std_hld(2:end)]
[ auc_men_hld(2:end) ' +/- ' auc_std_hld(2:end)]


%% Tracts
trc_nme = '005_TRACTS_Outcomes_all3T_temporal_toall_performance.csv';

trc_dta = cell2mat(mmil_readtext([ dta_loc '/' trc_nme ]));

acc_men_hld     = [ num2str(roundsd(mean(trc_dta(:,1)),2)) ];
    acc_std_hld = [ num2str(roundsd(std(trc_dta(:,1)),2)) ];
ppv_men_hld     = [ num2str(roundsd(mean(trc_dta(:,2)),2)) ];
    ppv_std_hld = [ num2str(roundsd(std(trc_dta(:,2)),2)) ];
npv_men_hld     = [ num2str(roundsd(mean(trc_dta(:,3)),2)) ];
    npv_std_hld = [ num2str(roundsd(std(trc_dta(:,3)),2)) ];
sns_men_hld     = [ num2str(roundsd(mean(trc_dta(:,4)),2)) ];
    sns_std_hld = [ num2str(roundsd(std(trc_dta(:,4)),2)) ];
spe_men_hld     = [ num2str(roundsd(mean(trc_dta(:,5)),2)) ];
    spe_std_hld = [ num2str(roundsd(std(trc_dta(:,5)),2)) ];
auc_men_hld     = [ num2str(roundsd(mean(trc_dta(:,6)),2)) ];
    auc_std_hld = [ num2str(roundsd(std(trc_dta(:,6)),2)) ];

[ acc_men_hld(2:end) ' +/- ' acc_std_hld(2:end)]
[ ppv_men_hld(2:end) ' +/- ' ppv_std_hld(2:end)]
[ npv_men_hld(2:end) ' +/- ' npv_std_hld(2:end)]
[ sns_men_hld(2:end) ' +/- ' sns_std_hld(2:end)]
[ spe_men_hld(2:end) ' +/- ' spe_std_hld(2:end)]
[ auc_men_hld(2:end) ' +/- ' auc_std_hld(2:end)]

%% Connectome
cnn_nme = '004_CONNECTOME_Outcomes_all3T_temporal_toall_performance.csv';

cnn_dta = cell2mat(mmil_readtext([ dta_loc '/' cnn_nme ]));

acc_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,1)),2)) ];
    acc_std_hld = [ num2str(roundsd(std(cnn_dta(:,1)),2)) ];
ppv_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,2)),2)) ];
    ppv_std_hld = [ num2str(roundsd(std(cnn_dta(:,2)),2)) ];
npv_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,3)),2)) ];
    npv_std_hld = [ num2str(roundsd(std(cnn_dta(:,3)),2)) ];
sns_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,4)),2)) ];
    sns_std_hld = [ num2str(roundsd(std(cnn_dta(:,4)),2)) ];
spe_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,5)),2)) ];
    spe_std_hld = [ num2str(roundsd(std(cnn_dta(:,5)),2)) ];
auc_men_hld     = [ num2str(roundsd(mean(cnn_dta(:,6)),2)) ];
    auc_std_hld = [ num2str(roundsd(std(cnn_dta(:,6)),2)) ];

[ acc_men_hld(2:end) ' +/- ' acc_std_hld(2:end)]
[ ppv_men_hld(2:end) ' +/- ' ppv_std_hld(2:end)]
[ npv_men_hld(2:end) ' +/- ' npv_std_hld(2:end)]
[ sns_men_hld(2:end) ' +/- ' sns_std_hld(2:end)]
[ spe_men_hld(2:end) ' +/- ' spe_std_hld(2:end)]
[ auc_men_hld(2:end) ' +/- ' auc_std_hld(2:end)]

%% Stat Comparisons
% ACC CNN VS CLN
[ ~ , pvl ] = ttest2( cnn_dta(:,1) , cln_dta(:,1) );

% ACC CNN VS TRC
[ ~ , pvl ] = ttest2( cnn_dta(:,1) , trc_dta(:,1) );

% AUC CNN VS CLN
[ ~ , pvl ] = ttest2( cnn_dta(:,6) , cln_dta(:,6) );

% AUC CNN VS TRC
[ ~ , pvl ] = ttest2( cnn_dta(:,6) , trc_dta(:,6) );









