ttt = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/leavitrunning.csv');

bin_cen = 0.5:1:99.5;

rel_avg = round(cell2mat(ttt(1,:)')*100);
rel_std = round(cell2mat(ttt(2,:)')*100);
shf_avg = round(cell2mat(ttt(3,:)')*100);
shf_std = round(cell2mat(ttt(4,:)')*100);

figure()
subplot(2,1,1)
hist(rel_avg,bin_cen)
line( [ mean(rel_avg) mean(rel_avg)] , [0 40] , 'Color',  'r', 'LineWidth', 2)
title('DAGNet : Original Sample Accuracy')
ylabel('iterations')

subplot(2,1,2)
hist(shf_avg,bin_cen)
line( [ mean(shf_avg) mean(shf_avg)] , [0 40] , 'Color',  'r', 'LineWidth', 2)
title('DAGNet : Shuffled Sample Accuracy')
ylabel('iterations')
xlabel('accuracy (%)')

tightfig()


