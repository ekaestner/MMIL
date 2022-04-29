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

fprintf ('\n\nCNN\n')
fprintf ('CNN - Accuracy, mean: %g\n',roundsd(mean(cnn_acc(:,cut_off)),3))
fprintf ('CNN - Accuracy, median: %g\n',median(cnn_acc(:,cut_off)))
fprintf ('CNN - Accuracy, std: %g\n',roundsd(std(cnn_acc(:,cut_off)),3))
fprintf ('CNN - Accuracy, max: %g\n',roundsd(max(cnn_acc(:,cut_off)),3))
fprintf ('\n')
fprintf ('CNN - AUC, mean: %g\n',roundsd(mean(cnn_auc),3))
fprintf ('CNN - AUC, median: %g\n',median(cnn_auc))
fprintf ('CNN - AUC, std: %g\n',roundsd(std(cnn_auc),3))
fprintf ('CNN - AUC, max: %g\n',roundsd(max(cnn_auc),3))
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
fprintf ('\n')
fprintf ('CNN - Left Predictive, mean: %g\n',roundsd(mean( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, median: %g\n',roundsd(median( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, std: %g\n',roundsd(std( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2)))),3))
fprintf ('CNN - Left Predictive, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,lft_pre(1)) ./ (cnn_con(bst_cnn_mdl,lft_pre(1))+cnn_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('CNN - Right Sensitivity, mean: %g\n',roundsd(mean( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, median: %g\n',roundsd(median( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, std: %g\n',roundsd(std( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2)))),3))
fprintf ('CNN - Right Sensitivity, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,rgh_sen(1)) ./ (cnn_con(bst_cnn_mdl,rgh_sen(1))+cnn_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('CNN - Right Predictive, mean: %g\n',roundsd(mean( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, median: %g\n',roundsd(median( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, std: %g\n',roundsd(std( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2)))),3))
fprintf ('CNN - Right Predictive, max: %g\n',roundsd( cnn_con(bst_cnn_mdl,rgh_pre(1)) ./ (cnn_con(bst_cnn_mdl,rgh_pre(1))+cnn_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('CNN - FDCI (left-sensitivity), mean: %g\n',sum(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))>( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))))
fprintf ('CNN - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))-( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))),3)*100)
fprintf ('CNN - Diff (left-sensitivity), median: %g\n',roundsd(median(( cnn_con(:,lft_sen(1)) ./ (cnn_con(:,lft_sen(1))+cnn_con(:,lft_sen(2))))-( cnn_con(:,rgh_sen(1)) ./ (cnn_con(:,rgh_sen(1))+cnn_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('CNN - FDCI (left-predictive), mean: %g\n',sum(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))>( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))))
fprintf ('CNN - Diff (left-predictive), mean: %g\n',roundsd(mean(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))-( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))),3)*100)
fprintf ('CNN - Diff (left-predictive), median: %g\n',roundsd(median(( cnn_con(:,lft_pre(1)) ./ (cnn_con(:,lft_pre(1))+cnn_con(:,lft_pre(2))))-( cnn_con(:,rgh_pre(1)) ./ (cnn_con(:,rgh_pre(1))+cnn_con(:,rgh_pre(2))))),3)*100)



fprintf ('\n\nShuffled-CNN\n')
fprintf ('Shuffled-CNN - Accuracy, mean: %g\n',roundsd(mean(shf_acc(:,cut_off)),3))
fprintf ('Shuffled-CNN - Accuracy, median: %g\n',median(shf_acc(:,cut_off)))
fprintf ('Shuffled-CNN - Accuracy, std: %g\n',roundsd(std(shf_acc(:,cut_off)),3))
fprintf ('Shuffled-CNN - Accuracy, max: %g\n',roundsd(max(shf_acc(:,cut_off)),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - AUC, mean: %g\n',roundsd(mean(shf_auc),3))
fprintf ('Shuffled-CNN - AUC, median: %g\n',median(shf_auc))
fprintf ('Shuffled-CNN - AUC, std: %g\n',roundsd(std(shf_auc),3))
fprintf ('Shuffled-CNN - AUC, max: %g\n',roundsd(max(shf_auc),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - L+: %g\n',roundsd(mean(shf_con(:,1)),3))
fprintf ('Shuffled-CNN - L-: %g\n',roundsd(mean(shf_con(:,2)),3))
fprintf ('Shuffled-CNN - R-: %g\n',roundsd(mean(shf_con(:,3)),3))
fprintf ('Shuffled-CNN - R+: %g\n',roundsd(mean(shf_con(:,4)),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Left Sensitivity, mean: %g\n',roundsd(mean( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, median: %g\n',roundsd(median( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, std: %g\n',roundsd(std( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2)))),3))
fprintf ('Shuffled-CNN - Left Sensitivity, max: %g\n',roundsd( shf_con(bst_cnn_mdl,lft_sen(1)) ./ (shf_con(bst_cnn_mdl,lft_sen(1))+shf_con(bst_cnn_mdl,lft_sen(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Left Predictive, mean: %g\n',roundsd(mean( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, median: %g\n',roundsd(median( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, std: %g\n',roundsd(std( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2)))),3))
fprintf ('Shuffled-CNN - Left Predictive, max: %g\n',roundsd( shf_con(bst_cnn_mdl,lft_pre(1)) ./ (shf_con(bst_cnn_mdl,lft_pre(1))+shf_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Right Sensitivity, mean: %g\n',roundsd(mean( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, median: %g\n',roundsd(median( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, std: %g\n',roundsd(std( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2)))),3))
fprintf ('Shuffled-CNN - Right Sensitivity, max: %g\n',roundsd( shf_con(bst_cnn_mdl,rgh_sen(1)) ./ (shf_con(bst_cnn_mdl,rgh_sen(1))+shf_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - Right Predictive, mean: %g\n',roundsd(mean( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, median: %g\n',roundsd(median( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, std: %g\n',roundsd(std( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2)))),3))
fprintf ('Shuffled-CNN - Right Predictive, max: %g\n',roundsd( shf_con(bst_cnn_mdl,rgh_pre(1)) ./ (shf_con(bst_cnn_mdl,rgh_pre(1))+shf_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('Shuffled-CNN - FDCI (left-sensitivity), mean: %g\n',sum(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))>( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))))
fprintf ('Shuffled-CNN - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))-( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))),3)*100)
fprintf ('Shuffled-CNN - Diff (left-sensitivity), median: %g\n',roundsd(median(( shf_con(:,lft_sen(1)) ./ (shf_con(:,lft_sen(1))+shf_con(:,lft_sen(2))))-( shf_con(:,rgh_sen(1)) ./ (shf_con(:,rgh_sen(1))+shf_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('Shuffled-CNN - FDCI (left-predictive), mean: %g\n',sum(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))>( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))))
fprintf ('Shuffled-CNN - Diff (left-predictive), mean: %g\n',roundsd(mean(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))-( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))),3)*100)
fprintf ('Shuffled-CNN - Diff (left-predictive), median: %g\n',roundsd(median(( shf_con(:,lft_pre(1)) ./ (shf_con(:,lft_pre(1))+shf_con(:,lft_pre(2))))-( shf_con(:,rgh_pre(1)) ./ (shf_con(:,rgh_pre(1))+shf_con(:,rgh_pre(2))))),3)*100)




fprintf ('\n\nLogistic\n')
fprintf ('Logistic - Accuracy, mean: %g\n',roundsd(mean(log_acc),3))
fprintf ('Logistic - Accuracy, median: %g\n',median(log_acc))
fprintf ('Logistic - Accuracy, std: %g\n',roundsd(std(log_acc),3))
fprintf ('Logistic - Accuracy, max: %g\n',roundsd(max(log_acc),3))
fprintf ('\n')
fprintf ('Logistic - AUC, mean: %g\n',roundsd(mean(log_auc),3))
fprintf ('Logistic - AUC, median: %g\n',median(log_auc))
fprintf ('Logistic - AUC, std: %g\n',roundsd(std(log_auc),3))
fprintf ('Logistic - AUC, max: %g\n',roundsd(max(log_auc),3))
fprintf ('\n')
fprintf ('Logistic - L+: %g\n',roundsd(mean(log_con(:,1)),3))
fprintf ('Logistic - L-: %g\n',roundsd(mean(log_con(:,2)),3))
fprintf ('Logistic - R-: %g\n',roundsd(mean(log_con(:,3)),3))
fprintf ('Logistic - R+: %g\n',roundsd(mean(log_con(:,4)),3))
fprintf ('\n')
fprintf ('Logistic - Left Sensitivity, mean: %g\n',roundsd(mean( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, median: %g\n',roundsd(median( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, std: %g\n',roundsd(std( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2)))),3))
fprintf ('Logistic - Left Sensitivity, max: %g\n',roundsd( log_con(bst_cnn_mdl,lft_sen(1)) ./ (log_con(bst_cnn_mdl,lft_sen(1))+log_con(bst_cnn_mdl,lft_sen(2))),3))
fprintf ('\n')
fprintf ('Logistic - Left Predictive, mean: %g\n',roundsd(mean( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, median: %g\n',roundsd(median( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, std: %g\n',roundsd(std( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2)))),3))
fprintf ('Logistic - Left Predictive, max: %g\n',roundsd( log_con(bst_cnn_mdl,lft_pre(1)) ./ (log_con(bst_cnn_mdl,lft_pre(1))+log_con(bst_cnn_mdl,lft_pre(2))),3))
fprintf ('\n')
fprintf ('Logistic - Right Sensitivity, mean: %g\n',roundsd(mean( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, median: %g\n',roundsd(median( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, std: %g\n',roundsd(std( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2)))),3))
fprintf ('Logistic - Right Sensitivity, max: %g\n',roundsd( log_con(bst_cnn_mdl,rgh_sen(1)) ./ (log_con(bst_cnn_mdl,rgh_sen(1))+log_con(bst_cnn_mdl,rgh_sen(2))),3))
fprintf ('\n')
fprintf ('Logistic - Right Predictive, mean: %g\n',roundsd(mean( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, median: %g\n',roundsd(median( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, std: %g\n',roundsd(std( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2)))),3))
fprintf ('Logistic - Right Predictive, max: %g\n',roundsd( log_con(bst_cnn_mdl,rgh_pre(1)) ./ (log_con(bst_cnn_mdl,rgh_pre(1))+log_con(bst_cnn_mdl,rgh_pre(2))),3))
fprintf ('\n')
fprintf ('Logistic - FDCI (left-sensitivity), mean: %g\n',sum(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))>( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))))
fprintf ('Logistic - Diff (left-sensitivity), mean: %g\n',roundsd(mean(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))-( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))),3)*100)
fprintf ('Logistic - Diff (left-sensitivity), median: %g\n',roundsd(median(( log_con(:,lft_sen(1)) ./ (log_con(:,lft_sen(1))+log_con(:,lft_sen(2))))-( log_con(:,rgh_sen(1)) ./ (log_con(:,rgh_sen(1))+log_con(:,rgh_sen(2))))),3)*100)
fprintf ('\n')
fprintf ('Logistic - FDCI (left-predictive), mean: %g\n',sum(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))>( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))))
fprintf ('Logistic - Diff (left-predictive), mean: %g\n',roundsd(mean(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))-( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))),3)*100)
fprintf ('Logistic - Diff (left-predictive), median: %g\n',roundsd(median(( log_con(:,lft_pre(1)) ./ (log_con(:,lft_pre(1))+log_con(:,lft_pre(2))))-( log_con(:,rgh_pre(1)) ./ (log_con(:,rgh_pre(1))+log_con(:,rgh_pre(2))))),3)*100)










