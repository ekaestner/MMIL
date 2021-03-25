% size 227,60
% sbj num = 73 79 70

ent_roi = [5 39];
con_row = find(strcmpi(sbj_nme(:,2),'HC'))-1;
epd_row = find(strcmpi(sbj_nme(:,2),'EPD_Old'))-1;
mci_row = find(strcmpi(sbj_nme(:,2),'MCI'))-1;


thk_roi_nme(ent_roi)

subplot(3,1,1)
scatter( thk_dta( con_row, ent_roi(1) ), thk_dta( con_row, ent_roi(2) ) )
subplot(3,1,2)
scatter( thk_dta( epd_row, ent_roi(1) ), thk_dta( epd_row, ent_roi(2) ) )
subplot(3,1,3)
scatter( thk_dta( mci_row, ent_roi(1) ), thk_dta( mci_row, ent_roi(2) ) )

figure()
subplot(3,1,1); hold on;
plot(cor_mtx{1}{1}(5,:))
plot(cor_mtx{1}{1}(39,:),'r')
subplot(3,1,2); hold on;
plot(cor_mtx{1}{2}(5,:))
plot(cor_mtx{1}{2}(39,:),'r')
subplot(3,1,3); hold on;
plot(cor_mtx{1}{3}(5,:))
plot(cor_mtx{1}{3}(39,:),'r')