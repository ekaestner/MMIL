clear; clc;

out_plt_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/figures/png/Figure4';

dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/performance_final';

%% Load Data & Setup
cut_off = 50;

% CNN
cnn_acc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_test.csv' ]))*100;
    cnn_acc = cnn_acc(:,cut_off);
cnn_con     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_con_mat.csv' ]));
cnn_auc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_auc.csv' ]))*100;

cnn_roc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_tpr.csv' ]));

% shuffle-CNN
log_acc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_log.csv' ]))*100;
log_con = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_con_mat_log.csv' ]));
log_auc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_auc_log.csv' ]))*100;

log_roc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_tpr_log.csv' ]));

% Sensitivity
lft_sen = [ 1 2 ];
rgh_sen = [ 4 3 ];
 
cnn_lft_sen = cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))*100;
cnn_rgh_sen = cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))*100;

log_lft_sen = log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))*100;
log_rgh_sen = log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))*100;

% Subtractions
sub_acc = cnn_acc - log_acc;
sub_auc = cnn_auc - log_auc;

sub_cnn_sen = cnn_lft_sen - cnn_rgh_sen;
sub_log_sen = log_lft_sen - log_rgh_sen;

% ROC
cnn_roc_lne = mean(cnn_roc)*100;
log_roc_lne = mean(log_roc)*100;

cnn_roc_top_std = cnn_roc_lne + std(cnn_roc*100);
cnn_roc_dwn_std = cnn_roc_lne - std(cnn_roc*100);

log_roc_top_std = log_roc_lne + std(log_roc*100);
log_roc_dwn_std = log_roc_lne - std(log_roc*100);

%% Subplot 1: Accuracies
fcfg = [];

fcfg.xdt = { 1       2 };
fcfg.ydt = { cnn_acc log_acc };

fcfg.fce_col     = { rgb('light purple')       rgb('light blue')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark purple') rgb('dark blue') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'cnn' 'log' };
fcfg.xlm = [ 0.5 6.5 ];
fcfg.ylb = {'Accuracy (%)'};
fcfg.ylm = [ 50 100 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Subplot1_Accuracy.png';

ejk_scatter(fcfg)

%% Subplot 2: AUC
fcfg = [];

fcfg.xdt = { 1       2 };
fcfg.ydt = { cnn_auc log_auc };

fcfg.fce_col     = { rgb('light purple')       rgb('light blue')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark purple') rgb('dark blue') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'cnn' 'log' };
fcfg.xlm = [ 0.5 6.5 ];
fcfg.ylb = {'AUC'};
fcfg.ylm = [ 50 100 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Subplot2_AUC.png';

ejk_scatter(fcfg)

%% Subplot 3: Accuracy, AUC Subtraction
fcfg = [];

fcfg.xdt = { 1       3 };
fcfg.ydt = { sub_acc sub_auc };

fcfg.fce_col     = { rgb('grey')        rgb('grey')  };
fcfg.edg_col     = { [0 0 0]            [0 0 0] };
fcfg.box_plt_col = { rgb('black')       rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'acc' 'auc' };
fcfg.xlm = [ 0.5 6.5 ];
fcfg.ylb = {'Subtraction (CNN - Shuffle-CNN)'};
fcfg.ylm = [ -20 20 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Subplot3_Subtraction.png';

ejk_scatter(fcfg)

%% Subplot 4: Sensitivty
fcfg = [];

fcfg.xdt = { 1           2           4           5 };
fcfg.ydt = { cnn_lft_sen cnn_rgh_sen log_lft_sen log_rgh_sen};

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('reddish grey') rgb('bluish grey') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  rgb('reddish grey')  rgb('bluish grey') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'clf' 'crg' 'slf' 'srg' };
fcfg.xlm = [ 0.5 6.5 ];
fcfg.ylb = {'Sensitivity (%)'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Subplot4_sensitivity.png';

ejk_scatter(fcfg)

%% Subplot 5: Sensitivity Subtraction
fcfg = [];

fcfg.xdt = { 1           3 };
fcfg.ydt = { sub_cnn_sen sub_log_sen };

fcfg.fce_col     = { rgb('light purple')       rgb('light blue')  };
fcfg.edg_col     = { [0 0 0]             [0 0 0] };
fcfg.box_plt_col = { rgb('dark purple') rgb('dark blue') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'cnn' 'log' };
fcfg.xlm = [ 0.5 6.5 ];
fcfg.ylb = {'Sensitivity Subtraction (L-TLE - R-TLE)'};
fcfg.ylm = [ -60 60 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'Subplot5_Subtraction.png';

ejk_scatter(fcfg)

%% Subplot 6: ROC curve average
figure('Visible','off'); xlim([0 100]); ylim([0 100]); hold on;

plot([0 100],[0 100],'k--','LineWidth',2)

patch( [ 1:100           fliplr(1:100) ], ...
       [ cnn_roc_top_std fliplr(cnn_roc_dwn_std) ], ...
       rgb('light purple'), ...
       'FaceAlpha',0.3,'EdgeColor','none')
patch( [ 1:100           fliplr(1:100) ], ...
       [ log_roc_top_std fliplr(log_roc_dwn_std) ], ...
       rgb('light blue'), ...
       'FaceAlpha',0.3,'EdgeColor','none')
   
plot(1:100,cnn_roc_lne,'Color',rgb('dark purple'),'LineWidth',3)
plot(1:100,log_roc_lne,'Color',rgb('dark blue'),'LineWidth',3)

print( gcf, [ out_plt_dir '/' 'Subplot6_ROC.png' ], '-dpng')
close all;

