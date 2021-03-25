%%
ttt = load('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/Graph_Scratch/Thickness_Desikan_Modified/GraphThickMatrix.mat');

[ min(ttt.cor_mtx{1}{1}(:)) ...
  max(ttt.cor_mtx{1}{1}(:)) ]

[ min(ttt.cor_mtx{1}{2}(:)) ...
  max(ttt.cor_mtx{1}{2}(:)) ]

[ min(ttt.cor_mtx{1}{3}(:)) ...
  max(ttt.cor_mtx{1}{3}(:)) ]

figure();
hist( ttt.cor_mtx{1}{3}(:),1000);

%% Add Colorbar
lin_spc = linspace(-.2,.7,10);
lin_col = hot(70);
lin_col = [ repmat(lin_col(1,:),30,1) ; lin_col ];

figure()
axes('OuterPosition',[.9 .05 .05 .9],'visible','off')
colormap(lin_col)
clb = colorbar('v','Position',[.9 .05 .05 .9]);
clb.TickLength = 0;
clb.TickLabels = cellfun(@num2str,num2cell(roundsd(lin_spc,2)),'uni',0);

print(gcf,[ '/home/ekaestne/PROJECTS/OUTPUT/epd_age_cnn' '/' 'colorbar' '.png'],'-dpng','-r200')
close all
