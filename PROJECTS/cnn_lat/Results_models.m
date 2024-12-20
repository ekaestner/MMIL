clear; clc;

dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/performance_final';

%% Load
lft_sbj_nme = mmil_readtext([ dta_loc '/' 'fn_left.csv' ]);
rgh_sbj_nme = mmil_readtext([ dta_loc '/' 'fn_right.csv' ]);

% CNN %%%%%%%%%%%%%%%%%%%%%%%
cnn_acc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_test.csv' ]))*100;
cnn_acc_val = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_val.csv' ]))*100;
cnn_con     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_con_mat.csv' ]));
cnn_auc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_auc.csv' ]))*100;

cnn_lft_sbj = mmil_readtext([ dta_loc '/' 'exp_6_scores_left_sbj.csv' ]);
cnn_rgh_sbj = mmil_readtext([ dta_loc '/' 'exp_6_scores_right_sbj.csv' ]);

cut_off = 50;

% Shuffled-CNN %%%%%%%%%%%%%%%%%%%%%%%
shf_acc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_score_test_rand.csv' ]))*100;
shf_con = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_con_mat_rand.csv' ]));
shf_auc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_auc_rand.csv' ]))*100;

% Logistic Hippocampus %%%%%%%%%%%%%%%%%%%%%%%
log_acc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_log.csv' ]))*100;
log_con = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_con_mat_log.csv' ]));
log_auc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_auc_log.csv' ]))*100;

%% Sub-categories
lft_sen = [ 1 2 ];
lft_pre = [ 1 3 ];
rgh_sen = [ 4 3 ];
rgh_pre = [ 4 2 ];

%% Summary Statistics %%%%%%%%%%%%%%%%%%%%%%%
[~, bst_cnn_mdl] = max(cnn_acc(:,cut_off));
[~, wst_cnn_mdl] = min(cnn_acc(:,cut_off));

fprintf ('\n\nCNN\n')
fprintf ('CNN - Accuracy, mean: %g\n',roundsd(mean(cnn_acc(:,cut_off)),3))
fprintf ('CNN - Accuracy, median: %g\n',median(cnn_acc(:,cut_off)))
fprintf ('CNN - Accuracy, std: %g\n',roundsd(std(cnn_acc(:,cut_off)),3))
fprintf ('CNN - Accuracy, max: %g\n',roundsd(cnn_acc(bst_cnn_mdl,cut_off),3))
fprintf ('CNN - Accuracy, min: %g\n',roundsd(cnn_acc(wst_cnn_mdl,cut_off),3))
fprintf ('\n')
fprintf ('CNN - AUC, mean: %g\n',roundsd(mean(cnn_auc),3))
fprintf ('CNN - AUC, median: %g\n',median(cnn_auc))
fprintf ('CNN - AUC, std: %g\n',roundsd(std(cnn_auc),3))
fprintf ('CNN - AUC, max: %g\n',roundsd(cnn_auc(bst_cnn_mdl),3))
fprintf ('CNN - AUC, min: %g\n',roundsd(cnn_auc(wst_cnn_mdl),3))
fprintf ('\n')
fprintf ('CNN - FDCI (shuffled-cnn): %g\n',sum(cnn_acc(:,cut_off)>shf_acc(:,cut_off)))
fprintf ('CNN - Diff (shuffled-cnn): %g\n',roundsd(mean(cnn_acc(:,cut_off)-shf_acc(:,cut_off)),3))
fprintf ('CNN - Diff (shuffled-cnn): %g\n',roundsd(median(cnn_acc(:,cut_off)-shf_acc(:,cut_off)),3))
fprintf ('\n')
fprintf ('CNN - FDCI (logistic): %g\n',sum(cnn_acc(:,cut_off)>log_acc))
fprintf ('CNN - Diff (logistic): %g\n',roundsd(mean(cnn_acc(:,cut_off)-log_acc),3))
fprintf ('CNN - Diff (logistic): %g\n',roundsd(median(cnn_acc(:,cut_off)-log_acc),3))
fprintf ('\n')
fprintf ('CNN - L+: %g\n',roundsd(mean(cnn_con(:,1)),3))
fprintf ('CNN - L-: %g\n',roundsd(mean(cnn_con(:,2)),3))
fprintf ('CNN - R-: %g\n',roundsd(mean(cnn_con(:,3)),3))
fprintf ('CNN - R+: %g\n',roundsd(mean(cnn_con(:,4)),3))
fprintf ('\n')
fprintf ('CNN - Left Sensitivity, mean: %g\n',roundsd(mean( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))),3))
fprintf ('CNN - Left Sensitivity, median: %g\n',roundsd(median( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))),3))
fprintf ('CNN - Left Sensitivity, std: %g\n',roundsd(std( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))),3))
fprintf ('CNN - Left Sensitivity, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,lft_sen(1)) ./ (cnn_con(bst_cnn_mdl,lft_sen(1))+cnn_con(bst_cnn_mdl,lft_sen(2))),3))
fprintf ('CNN - Left Sensitivity, min: %g\n',roundsd( cnn_con(wst_cnn_mdl,lft_sen(1)) ./ (cnn_con(wst_cnn_mdl,lft_sen(1))+cnn_con(wst_cnn_mdl,lft_sen(2))),3))
fprintf ('\n')
fprintf ('CNN - Left Predictive, mean: %g\n',roundsd(mean( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, median: %g\n',roundsd(median( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, std: %g\n',roundsd(std( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,lft_pre(1)) ./ (cnn_con(bst_cnn_mdl,lft_pre(1))+cnn_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('CNN - Left Predictive, min: %g\n',roundsd( cnn_con(wst_cnn_mdl,lft_pre(1)) ./ (cnn_con(wst_cnn_mdl,lft_pre(1))+cnn_con(wst_cnn_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('CNN - Right Sensitivity, mean: %g\n',roundsd(mean( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, median: %g\n',roundsd(median( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, std: %g\n',roundsd(std( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,rgh_sen(1)) ./ (cnn_con(bst_cnn_mdl,rgh_sen(1))+cnn_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('CNN - Right Sensitivity, min: %g\n',roundsd( cnn_con(wst_cnn_mdl,rgh_sen(1)) ./ (cnn_con(wst_cnn_mdl,rgh_sen(1))+cnn_con(wst_cnn_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('CNN - Right Predictive, mean: %g\n',roundsd(mean( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, median: %g\n',roundsd(median( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, std: %g\n',roundsd(std( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,rgh_pre(1)) ./ (cnn_con(bst_cnn_mdl,rgh_pre(1))+cnn_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('CNN - Right Predictive, min: %g\n',roundsd( cnn_con(wst_cnn_mdl,rgh_pre(1)) ./ (cnn_con(wst_cnn_mdl,rgh_pre(1))+cnn_con(wst_cnn_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('CNN - FDCI (left-sensitivity), mean: %g\n',sum(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))>( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))))
fprintf ('CNN - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))-( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))),3)*100)
fprintf ('CNN - Diff (left-sensitivity), median: %g\n',roundsd(median(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))-( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('CNN - FDCI (left-predictive), mean: %g\n',sum(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))>( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))))
fprintf ('CNN - Diff (left-predictive), mean: %g\n',roundsd(mean(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))-( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))),3)*100)
fprintf ('CNN - Diff (left-predictive), median: %g\n',roundsd(median(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))-( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))),3)*100)



[~, bst_shf_mdl] = max(shf_acc(:,cut_off));
[~, wst_shf_mdl] = min(shf_acc(:,cut_off));

fprintf ('\n\nShuffled-CNN\n')
fprintf ('Shuffled-CNN - Accuracy, mean: %g\n',roundsd(mean(shf_acc(:,cut_off)),3))
fprintf ('Shuffled-CNN - Accuracy, median: %g\n',median(shf_acc(:,cut_off)))
fprintf ('Shuffled-CNN - Accuracy, std: %g\n',roundsd(std(shf_acc(:,cut_off)),3))
fprintf ('Shuffled-CNN - Accuracy, max (cnn indice): %g\n',roundsd(shf_acc(bst_cnn_mdl,cut_off),3))
fprintf ('Shuffled-CNN - Accuracy, min (cnn indice): %g\n',roundsd(shf_acc(wst_cnn_mdl,cut_off),3))
fprintf ('Shuffled-CNN - Accuracy, max: %g\n',roundsd(shf_acc(bst_shf_mdl,cut_off),3))
fprintf ('Shuffled-CNN - Accuracy, min: %g\n',roundsd(shf_acc(wst_shf_mdl,cut_off),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - AUC, mean: %g\n',roundsd(mean(shf_auc),3))
fprintf ('Shuffled-CNN - AUC, median: %g\n',median(shf_auc))
fprintf ('Shuffled-CNN - AUC, std: %g\n',roundsd(std(shf_auc),3))
fprintf ('Shuffled-CNN - AUC, max (cnn indice): %g\n',roundsd(shf_auc(bst_cnn_mdl),3))
fprintf ('Shuffled-CNN - AUC, min (cnn indice): %g\n',roundsd(shf_auc(wst_cnn_mdl),3))
fprintf ('Shuffled-CNN - AUC, max: %g\n',roundsd(shf_auc(bst_shf_mdl),3))
fprintf ('Shuffled-CNN - AUC, min: %g\n',roundsd(shf_auc(wst_shf_mdl),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - L+: %g\n',roundsd(mean(shf_con(:,1)),3))
fprintf ('Shuffled-CNN - L-: %g\n',roundsd(mean(shf_con(:,2)),3))
fprintf ('Shuffled-CNN - R-: %g\n',roundsd(mean(shf_con(:,3)),3))
fprintf ('Shuffled-CNN - R+: %g\n',roundsd(mean(shf_con(:,4)),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Left Sensitivity, mean: %g\n',roundsd(mean( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, median: %g\n',roundsd(median( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, std: %g\n',roundsd(std( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, max (cnn indice): %g\n',roundsd( shf_con(bst_cnn_mdl,lft_sen(1)) ./ (shf_con(bst_cnn_mdl,lft_sen(1))+shf_con(bst_cnn_mdl,lft_sen(2))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, min (cnn indice): %g\n',roundsd( shf_con(wst_cnn_mdl,lft_sen(1)) ./ (shf_con(wst_cnn_mdl,lft_sen(1))+shf_con(wst_cnn_mdl,lft_sen(2))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, max: %g\n',roundsd( shf_con(bst_shf_mdl,lft_sen(1)) ./ (shf_con(bst_shf_mdl,lft_sen(1))+shf_con(bst_shf_mdl,lft_sen(2))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, min: %g\n',roundsd( shf_con(wst_shf_mdl,lft_sen(1)) ./ (shf_con(wst_shf_mdl,lft_sen(1))+shf_con(wst_shf_mdl,lft_sen(2))),3))

fprintf ('\n')
fprintf ('Shuffled-CNN - Left Predictive, mean: %g\n',roundsd(mean( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, median: %g\n',roundsd(median( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, std: %g\n',roundsd(std( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, max (cnn indice): %g\n',roundsd( shf_con(bst_cnn_mdl,lft_pre(1)) ./ (shf_con(bst_cnn_mdl,lft_pre(1))+shf_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('Shuffled-CNN - Left Predictive, min (cnn indice): %g\n',roundsd( shf_con(wst_cnn_mdl,lft_pre(1)) ./ (shf_con(wst_cnn_mdl,lft_pre(1))+shf_con(wst_cnn_mdl,lft_pre(2))),3))
fprintf ('Shuffled-CNN - Left Predictive, max: %g\n',roundsd( shf_con(bst_shf_mdl,lft_pre(1)) ./ (shf_con(bst_shf_mdl,lft_pre(1))+shf_con(bst_shf_mdl,lft_pre(2))),3))
fprintf ('Shuffled-CNN - Left Predictive, min: %g\n',roundsd( shf_con(wst_shf_mdl,lft_pre(1)) ./ (shf_con(wst_shf_mdl,lft_pre(1))+shf_con(wst_shf_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Right Sensitivity, mean: %g\n',roundsd(mean( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, median: %g\n',roundsd(median( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, std: %g\n',roundsd(std( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, max (cnn indice): %g\n',roundsd( shf_con(bst_cnn_mdl,rgh_sen(1)) ./ (shf_con(bst_cnn_mdl,rgh_sen(1))+shf_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, min (cnn indice): %g\n',roundsd( shf_con(wst_cnn_mdl,rgh_sen(1)) ./ (shf_con(wst_cnn_mdl,rgh_sen(1))+shf_con(wst_cnn_mdl,rgh_sen(2))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, max: %g\n',roundsd( shf_con(bst_shf_mdl,rgh_sen(1)) ./ (shf_con(bst_shf_mdl,rgh_sen(1))+shf_con(bst_shf_mdl,rgh_sen(2))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, min: %g\n',roundsd( shf_con(wst_shf_mdl,rgh_sen(1)) ./ (shf_con(wst_shf_mdl,rgh_sen(1))+shf_con(wst_shf_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Right Predictive, mean: %g\n',roundsd(mean( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, median: %g\n',roundsd(median( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, std: %g\n',roundsd(std( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, max (cnn indice): %g\n',roundsd( shf_con(bst_cnn_mdl,rgh_pre(1)) ./ (shf_con(bst_cnn_mdl,rgh_pre(1))+shf_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('Shuffled-CNN - Right Predictive, min (cnn indice): %g\n',roundsd( shf_con(wst_cnn_mdl,rgh_pre(1)) ./ (shf_con(wst_cnn_mdl,rgh_pre(1))+shf_con(wst_cnn_mdl,rgh_pre(2))),3))
fprintf ('Shuffled-CNN - Right Predictive, max: %g\n',roundsd( shf_con(bst_shf_mdl,rgh_pre(1)) ./ (shf_con(bst_shf_mdl,rgh_pre(1))+shf_con(bst_shf_mdl,rgh_pre(2))),3))
fprintf ('Shuffled-CNN - Right Predictive, min: %g\n',roundsd( shf_con(wst_shf_mdl,rgh_pre(1)) ./ (shf_con(wst_shf_mdl,rgh_pre(1))+shf_con(wst_shf_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - FDCI (left-sensitivity), mean: %g\n',sum(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))>( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))))
fprintf ('Shuffled-CNN - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))-( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))),3)*100)
fprintf ('Shuffled-CNN - Diff (left-sensitivity), median: %g\n',roundsd(median(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))-( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('Shuffled-CNN - FDCI (left-predictive), mean: %g\n',sum(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))>( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))))
fprintf ('Shuffled-CNN - Diff (left-predictive), mean: %g\n',roundsd(mean(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))-( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))),3)*100)
fprintf ('Shuffled-CNN - Diff (left-predictive), median: %g\n',roundsd(median(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))-( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))),3)*100)



[~, bst_log_mdl] = max(log_acc);
[~, wst_log_mdl] = min(log_acc);

fprintf ('\n\nLogistic\n')
fprintf ('Logistic - Accuracy, mean: %g\n',roundsd(mean(log_acc),3))
fprintf ('Logistic - Accuracy, median: %g\n',median(log_acc))
fprintf ('Logistic - Accuracy, std: %g\n',roundsd(std(log_acc),3))
fprintf ('Logistic - Accuracy, max (cnn indice): %g\n',roundsd(log_acc(bst_cnn_mdl),3))
fprintf ('Logistic - Accuracy, min (cnn indice): %g\n',roundsd(log_acc(wst_cnn_mdl),3))
fprintf ('Logistic - Accuracy, max: %g\n',roundsd(log_acc(bst_log_mdl),3))
fprintf ('Logistic - Accuracy, min: %g\n',roundsd(log_acc(wst_log_mdl),3))
fprintf ('\n')
fprintf ('Logistic - AUC, mean: %g\n',roundsd(mean(log_auc),3))
fprintf ('Logistic - AUC, median: %g\n',median(log_auc))
fprintf ('Logistic - AUC, std: %g\n',roundsd(std(log_auc),3))
fprintf ('Logistic - AUC, max (cnn indice): %g\n',roundsd(log_auc(bst_cnn_mdl),3))
fprintf ('Logistic - AUC, min (cnn indice): %g\n',roundsd(log_auc(wst_cnn_mdl),3))
fprintf ('Logistic - AUC, max: %g\n',roundsd(log_auc(bst_log_mdl),3))
fprintf ('Logistic - AUC, min: %g\n',roundsd(log_auc(wst_log_mdl),3))
fprintf ('\n')
fprintf ('Logistic - L+: %g\n',roundsd(mean(log_con(:,1)),3))
fprintf ('Logistic - L-: %g\n',roundsd(mean(log_con(:,2)),3))
fprintf ('Logistic - R-: %g\n',roundsd(mean(log_con(:,3)),3))
fprintf ('Logistic - R+: %g\n',roundsd(mean(log_con(:,4)),3))
fprintf ('\n')
fprintf ('Logistic - Left Sensitivity, mean: %g\n',roundsd(mean( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, median: %g\n',roundsd(median( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, std: %g\n',roundsd(std( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, max (cnn indice): %g\n',roundsd( log_con(bst_cnn_mdl,lft_sen(1)) ./ (log_con(bst_cnn_mdl,lft_sen(1))+log_con(bst_cnn_mdl,lft_sen(2))),3))
fprintf ('Logistic - Left Sensitivity, min (cnn indice): %g\n',roundsd( log_con(wst_cnn_mdl,lft_sen(1)) ./ (log_con(wst_cnn_mdl,lft_sen(1))+log_con(wst_cnn_mdl,lft_sen(2))),3))
fprintf ('Logistic - Left Sensitivity, max: %g\n',roundsd( log_con(bst_log_mdl,lft_sen(1)) ./ (log_con(bst_log_mdl,lft_sen(1))+log_con(bst_log_mdl,lft_sen(2))),3))
fprintf ('Logistic - Left Sensitivity, min: %g\n',roundsd( log_con(wst_log_mdl,lft_sen(1)) ./ (log_con(wst_log_mdl,lft_sen(1))+log_con(wst_log_mdl,lft_sen(2))),3))
fprintf ('\n')
fprintf ('Logistic - Left Predictive, mean: %g\n',roundsd(mean( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, median: %g\n',roundsd(median( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, std: %g\n',roundsd(std( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, max (cnn indice): %g\n',roundsd( log_con(bst_cnn_mdl,lft_pre(1)) ./ (log_con(bst_cnn_mdl,lft_pre(1))+log_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('Logistic - Left Predictive, min (cnn indice): %g\n',roundsd( log_con(wst_cnn_mdl,lft_pre(1)) ./ (log_con(wst_cnn_mdl,lft_pre(1))+log_con(wst_cnn_mdl,lft_pre(2))),3))
fprintf ('Logistic - Left Predictive, max: %g\n',roundsd( log_con(bst_log_mdl,lft_pre(1)) ./ (log_con(bst_log_mdl,lft_pre(1))+log_con(bst_log_mdl,lft_pre(2))),3))
fprintf ('Logistic - Left Predictive, min: %g\n',roundsd( log_con(wst_log_mdl,lft_pre(1)) ./ (log_con(wst_log_mdl,lft_pre(1))+log_con(wst_log_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('Logistic - Right Sensitivity, mean: %g\n',roundsd(mean( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, median: %g\n',roundsd(median( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, std: %g\n',roundsd(std( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, max (cnn indice): %g\n',roundsd( log_con(bst_cnn_mdl,rgh_sen(1)) ./ (log_con(bst_cnn_mdl,rgh_sen(1))+log_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('Logistic - Right Sensitivity, min (cnn indice): %g\n',roundsd( log_con(wst_cnn_mdl,rgh_sen(1)) ./ (log_con(wst_cnn_mdl,rgh_sen(1))+log_con(wst_cnn_mdl,rgh_sen(2))),3))
fprintf ('Logistic - Right Sensitivity, max: %g\n',roundsd( log_con(bst_log_mdl,rgh_sen(1)) ./ (log_con(bst_log_mdl,rgh_sen(1))+log_con(bst_log_mdl,rgh_sen(2))),3))
fprintf ('Logistic - Right Sensitivity, min: %g\n',roundsd( log_con(wst_log_mdl,rgh_sen(1)) ./ (log_con(wst_log_mdl,rgh_sen(1))+log_con(wst_log_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('Logistic - Right Predictive, mean: %g\n',roundsd(mean( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, median: %g\n',roundsd(median( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, std: %g\n',roundsd(std( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, max (cnn indice): %g\n',roundsd( log_con(bst_cnn_mdl,rgh_pre(1)) ./ (log_con(bst_cnn_mdl,rgh_pre(1))+log_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('Logistic - Right Predictive, min (cnn indice): %g\n',roundsd( log_con(wst_cnn_mdl,rgh_pre(1)) ./ (log_con(wst_cnn_mdl,rgh_pre(1))+log_con(wst_cnn_mdl,rgh_pre(2))),3))
fprintf ('Logistic - Right Predictive, max: %g\n',roundsd( log_con(bst_log_mdl,rgh_pre(1)) ./ (log_con(bst_log_mdl,rgh_pre(1))+log_con(bst_log_mdl,rgh_pre(2))),3))
fprintf ('Logistic - Right Predictive, min: %g\n',roundsd( log_con(wst_log_mdl,rgh_pre(1)) ./ (log_con(wst_log_mdl,rgh_pre(1))+log_con(wst_log_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('Logistic - FDCI (left-sensitivity), mean: %g\n',sum(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))>( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))))
fprintf ('Logistic - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))-( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))),3)*100)
fprintf ('Logistic - Diff (left-sensitivity), median: %g\n',roundsd(median(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))-( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('Logistic - FDCI (left-predictive), mean: %g\n',sum(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))>( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))))
fprintf ('Logistic - Diff (left-predictive), mean: %g\n',roundsd(mean(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))-( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))),3)*100)
fprintf ('Logistic - Diff (left-predictive), median: %g\n',roundsd(median(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))-( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))),3)*100)

%% Table
clear tbl_out

tbl_out{1,1} = [ num2str(roundsd(mean(cnn_acc(:,cut_off)),3)) '% (' num2str(roundsd(std(cnn_acc(:,cut_off)),3)) '%)' ];
tbl_out{1,2} = [ num2str(roundsd(mean(cnn_auc/100),3)) ' (' num2str(roundsd(std(cnn_auc/100),3)) ')' ];
tbl_out{1,3} = [ num2str(roundsd(mean( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))),3)*100) '% (' num2str(roundsd(std( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2)))),3)*100) '%)' ];
tbl_out{1,4} = [ num2str(roundsd(mean( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3)*100) '% (' num2str(roundsd(std( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3)*100) '%)' ];

tbl_out{2,1} = [ num2str(roundsd(cnn_acc(bst_cnn_mdl,cut_off),3)) '%' ];
tbl_out{2,2} = [ num2str(roundsd(cnn_auc(bst_cnn_mdl)/100,3))  ];
tbl_out{2,3} = [ num2str(roundsd( cnn_con(bst_cnn_mdl,lft_sen(1)) ./ (cnn_con(bst_cnn_mdl,lft_sen(1))+cnn_con(bst_cnn_mdl,lft_sen(2))),3)*100) '%' ];
tbl_out{2,4} = [ num2str(roundsd( cnn_con(bst_cnn_mdl,rgh_sen(1)) ./ (cnn_con(bst_cnn_mdl,rgh_sen(1))+cnn_con(bst_cnn_mdl,rgh_sen(2))),3)*100) '%' ];

tbl_out{3,1} = [ num2str(roundsd(cnn_acc(wst_cnn_mdl,cut_off),3)) '%' ];
tbl_out{3,2} = [ num2str(roundsd(cnn_auc(wst_cnn_mdl)/100,3))  ];
tbl_out{3,3} = [ num2str(roundsd( cnn_con(wst_cnn_mdl,lft_sen(1)) ./ (cnn_con(wst_cnn_mdl,lft_sen(1))+cnn_con(wst_cnn_mdl,lft_sen(2))),3)*100) '%' ];
tbl_out{3,4} = [ num2str(roundsd( cnn_con(wst_cnn_mdl,rgh_sen(1)) ./ (cnn_con(wst_cnn_mdl,rgh_sen(1))+cnn_con(wst_cnn_mdl,rgh_sen(2))),3)*100) '%' ];

tbl_out{4,1} = [ num2str(roundsd(mean(shf_acc(:,cut_off)),3)) '% (' num2str(roundsd(std(shf_acc(:,cut_off)),3)) '%)' ];
tbl_out{4,2} = [ num2str(roundsd(mean(shf_auc/100),3)) ' (' num2str(roundsd(std(shf_auc/100),3)) ')' ];
tbl_out{4,3} = [ num2str(roundsd(mean( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3)) '% (' num2str(roundsd(std( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3)) '%)' ];
tbl_out{4,4} = [ num2str(roundsd(mean( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3)) '% (' num2str(roundsd(std( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3)) '%)' ];

tbl_out{5,1} = [ num2str(roundsd(shf_acc(bst_cnn_mdl,cut_off),3)) '%' ];
tbl_out{5,2} = [ num2str(roundsd(shf_auc(bst_cnn_mdl)/100,3))  ];
tbl_out{5,3} = [ num2str(roundsd( shf_con(bst_cnn_mdl,lft_sen(1)) ./ (shf_con(bst_cnn_mdl,lft_sen(1))+shf_con(bst_cnn_mdl,lft_sen(2))),3)) '%' ];
tbl_out{5,4} = [ num2str(roundsd( shf_con(bst_cnn_mdl,rgh_sen(1)) ./ (shf_con(bst_cnn_mdl,rgh_sen(1))+shf_con(bst_cnn_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{6,1} = [ num2str(roundsd(shf_acc(wst_cnn_mdl,cut_off),3)) '%' ];
tbl_out{6,2} = [ num2str(roundsd(shf_auc(wst_cnn_mdl)/100,3))  ];
tbl_out{6,3} = [ num2str(roundsd( shf_con(wst_cnn_mdl,lft_sen(1)) ./ (shf_con(wst_cnn_mdl,lft_sen(1))+shf_con(wst_cnn_mdl,lft_sen(2))),3)) '%' ];
tbl_out{6,4} = [ num2str(roundsd( shf_con(wst_cnn_mdl,rgh_sen(1)) ./ (shf_con(wst_cnn_mdl,rgh_sen(1))+shf_con(wst_cnn_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{7,1} = [ num2str(roundsd(shf_acc(bst_shf_mdl,cut_off),3)) '%' ];
tbl_out{7,2} = [ num2str(roundsd(shf_auc(bst_shf_mdl)/100,3))  ];
tbl_out{7,3} = [ num2str(roundsd( shf_con(bst_shf_mdl,lft_sen(1)) ./ (shf_con(bst_shf_mdl,lft_sen(1))+shf_con(bst_shf_mdl,lft_sen(2))),3)) '%' ];
tbl_out{7,4} = [ num2str(roundsd( shf_con(bst_shf_mdl,rgh_sen(1)) ./ (shf_con(bst_shf_mdl,rgh_sen(1))+shf_con(bst_shf_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{8,1} = [ num2str(roundsd(shf_acc(wst_shf_mdl,cut_off),3)) '%' ];
tbl_out{8,2} = [ num2str(roundsd(shf_auc(wst_shf_mdl)/100,3))  ];
tbl_out{8,3} = [ num2str(roundsd( shf_con(wst_shf_mdl,lft_sen(1)) ./ (shf_con(wst_shf_mdl,lft_sen(1))+shf_con(wst_shf_mdl,lft_sen(2))),3)) '%' ];
tbl_out{8,4} = [ num2str(roundsd( shf_con(wst_shf_mdl,rgh_sen(1)) ./ (shf_con(wst_shf_mdl,rgh_sen(1))+shf_con(wst_shf_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{9,1} = [ num2str(roundsd(mean(log_acc),3)) '% (' num2str(roundsd(std(log_acc),3)) '%)' ];
tbl_out{9,2} = [ num2str(roundsd(mean(log_auc),3)) ' (' num2str(roundsd(std(log_auc),3)) ')' ];
tbl_out{9,3} = [ num2str(roundsd(mean( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3)) '% (' num2str(roundsd(std( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3)) '%)' ];
tbl_out{9,4} = [ num2str(roundsd(mean( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3)) '% (' num2str(roundsd(std( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3)) '%)' ];

tbl_out{10,1} = [ num2str(roundsd(log_acc(bst_cnn_mdl),3)) '%' ];
tbl_out{10,2} = [ num2str(roundsd(log_auc(bst_cnn_mdl),3))  ];
tbl_out{10,3} = [ num2str(roundsd( log_con(bst_cnn_mdl,lft_sen(1)) ./ (log_con(bst_cnn_mdl,lft_sen(1))+log_con(bst_cnn_mdl,lft_sen(2))),3)) '%' ];
tbl_out{10,4} = [ num2str(roundsd( log_con(bst_cnn_mdl,rgh_sen(1)) ./ (log_con(bst_cnn_mdl,rgh_sen(1))+log_con(bst_cnn_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{11,1} = [ num2str(roundsd(log_acc(wst_cnn_mdl),3)) '%' ];
tbl_out{11,2} = [ num2str(roundsd(log_auc(wst_cnn_mdl),3))  ];
tbl_out{11,3} = [ num2str(roundsd( log_con(wst_cnn_mdl,lft_sen(1)) ./ (log_con(wst_cnn_mdl,lft_sen(1))+log_con(wst_cnn_mdl,lft_sen(2))),3)) '%' ];
tbl_out{11,4} = [ num2str(roundsd( log_con(wst_cnn_mdl,rgh_sen(1)) ./ (log_con(wst_cnn_mdl,rgh_sen(1))+log_con(wst_cnn_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{12,1} = [ num2str(roundsd(log_acc(bst_log_mdl),3)) '%' ];
tbl_out{12,2} = [ num2str(roundsd(log_auc(bst_log_mdl),3))  ];
tbl_out{12,3} = [ num2str(roundsd( log_con(bst_log_mdl,lft_sen(1)) ./ (log_con(bst_log_mdl,lft_sen(1))+log_con(bst_log_mdl,lft_sen(2))),3)) '%' ];
tbl_out{12,4} = [ num2str(roundsd( log_con(bst_log_mdl,rgh_sen(1)) ./ (log_con(bst_log_mdl,rgh_sen(1))+log_con(bst_log_mdl,rgh_sen(2))),3)) '%' ];

tbl_out{13,1} = [ num2str(roundsd(log_acc(wst_log_mdl),3)) '%' ];
tbl_out{13,2} = [ num2str(roundsd(log_auc(wst_log_mdl),3))  ];
tbl_out{13,3} = [ num2str(roundsd( log_con(wst_log_mdl,lft_sen(1)) ./ (log_con(wst_log_mdl,lft_sen(1))+log_con(wst_log_mdl,lft_sen(2))),3)) '%' ];
tbl_out{13,4} = [ num2str(roundsd( log_con(wst_log_mdl,rgh_sen(1)) ./ (log_con(wst_log_mdl,rgh_sen(1))+log_con(wst_log_mdl,rgh_sen(2))),3)) '%' ];

cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/submission/2022_09_Neurology/Reviews/Figures/NewPieces' '/' 'table2_revision.csv'],tbl_out)






