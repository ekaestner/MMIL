clear; clc;

%% ROC Curves
% Clinical
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/003_CLINICAL_Outcomes_all3T_roc.csv');

xax{1} = cell2mat(dta_hld(1,:));
yax{1} = cell2mat(dta_hld(2,:));
mld_cmp_col(1,:) = rgb('greyish blue');

% Tracts
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/002_TRACTS_Outcomes_all3T_roc.csv');

xax{2} = cell2mat(dta_hld(1,:));
yax{2} = cell2mat(dta_hld(2,:));
mld_cmp_col(2,:) = rgb('dark maroon');

% Connectome
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_roc.csv');

xax{3} = cell2mat(dta_hld(1,:));
yax{3} = cell2mat(dta_hld(2,:));
mld_cmp_col(3,:) = rgb('purple');

% Plots
figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
hold on;

plot([0 1],[0 1],'Color',rgb('grey'),'LineWidth',3,'LineStyle','--')

for iSP = 1:numel(xax)
    plot(xax{iSP},yax{iSP},'Color',rgb('dark grey'),'LineWidth',10)
    plot(xax{iSP},yax{iSP},'Color',rgb('light grey'),'LineWidth',9)
    plot(xax{iSP},yax{iSP},'Color',mld_cmp_col(iSP,:),'LineWidth',7)
end

set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'LineWidth',4)
set(gca,'FontSize',26)

tightfig();
set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/Figure3' '/' 'Plot_ROC_FIGURE.png'],'-dpng')
close all

%% AUC Bar Plot
% Clinical
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/003_CLINICAL_Outcomes_all3T_performance.csv');

yax{1} = dta_hld{6};
mld_cmp_col(1,:) = rgb('greyish blue');

% Tracts
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/002_TRACTS_Outcomes_all3T_performance.csv');

yax{2} = dta_hld{6};
mld_cmp_col(2,:) = rgb('dark maroon');

% Connectome
dta_hld = mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_performance.csv');

yax{3} = dta_hld{6};
mld_cmp_col(3,:) = rgb('purple');

% Plot
figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
hold on;

for iSP = 1:numel(yax)
    bar_hld(iSP) = barh(numel(yax) - (iSP-1),yax{iSP});
    bar_hld(iSP).FaceColor = mld_cmp_col(iSP,:);
    bar_hld(iSP).EdgeColor = rgb('black');
    bar_hld(iSP).LineWidth = 3.5;
end
set(gca,'LineWidth',4)
set(gca,'FontSize',26)
set(gca,'Xlim',[0.5 0.85])
set(gca,'XTick',[0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85])
set(gca,'YTick',1:numel(yax))
tightfig();
set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

tightfig();
set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);

print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/Figure3' '/' 'Plot_AUC_FIGURE.png'],'-dpng')
close all