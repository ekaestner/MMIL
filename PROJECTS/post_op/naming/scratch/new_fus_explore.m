clear; clc;

%% Setup subjects
out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/NewWmparcHold/Checks/';

load( [ '/home/ekaestne/PROJECTS/OUTPUT' '/' 'PostOperative/Naming' '/' 'groups.mat' ] );

cog_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'Cognitive'                '_' 'QC' '.csv'];;
cog_dta = mmil_readtext(cog_dta_nme);
    cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
    cog_dta_sbj = cog_dta(2:end,1);
    cog_dta     = cell2mat(cog_dta(2:end,2:end));

wmp_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'];
wmp_dta = mmil_readtext(wmp_dta_nme);
    wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
    wmp_dta_sbj = wmp_dta(2:end,1);
    wmp_dta_rcn = wmp_dta(2:end,2);
    wmp_dta     = wmp_dta(2:end,5:end);

new_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/NewWmparcHold/new_roi.csv'];
new_dta = mmil_readtext(new_dta_nme);
    new_dta_col = ejk_fix_column_names(new_dta(1,3:end));
    new_dta_sbj = new_dta(2:end,1);
    new_dta_rcn = new_dta(2:end,2);
    new_dta     = new_dta(2:end,3:end);
    
%% Check ROIs
% Original ROIs %%%%%%%%%%%%%
org_roi_ind = [ 4 11 14 23 30 37 44 47 56 63];
new_roi_ind = [ 1 2  3  4  5  10 11 12 13 14 ];
[ wmp_dta_col(org_roi_ind)' new_dta_col(new_roi_ind)'];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)-1
    subplot(3,3,iP)
    scatter( cell2mat(wmp_dta( grp.tle_post_3T_ATLonly_left, org_roi_ind(iP) )), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title(mmil_spec_char(wmp_dta_col{org_roi_ind(iP)},{'_'},{' '}))
end
tightfig();
print(gcf,[ out_dir '/' 'original_ROIs.png'],'-dpng');
close all

% New ROIs %%%%%%%%%%%%%
org_roi_ind = [ 6 6 6 39 39 39 ];
new_roi_ind = [ 6 7 8 15 16 17 ];
[ wmp_dta_col(org_roi_ind)' new_dta_col(new_roi_ind)'];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)
    subplot(3,3,iP)
    scatter( cell2mat(wmp_dta( grp.tle_post_3T_ATLonly_left, org_roi_ind(iP) )), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title(mmil_spec_char(new_dta_col{new_roi_ind(iP)},{'_'},{' '}))
end
tightfig();
print(gcf,[ out_dir '/' 'new_ROIs.png'],'-dpng');
close all

%% Check stats
% New ROIs/BNT-POST %%%%%%%%%%%%%
cog_roi_ind = [ 4 4 4 4  4  4 ];
new_roi_ind = [ 6 7 8 15 16 17 ];
[ wmp_dta_col(cog_roi_ind)' new_dta_col(new_roi_ind)'];

org_roi_ind = [ 6 39 ];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)
    subplot(3,3,iP)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(new_dta_col{new_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

for iP = 1:numel(org_roi_ind)
    subplot(3,2,iP+4)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(wmp_dta_col{org_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

tightfig();
print(gcf,[ out_dir '/' 'new_ROIs_BNT_post.png'],'-dpng');
close all

% New ROIs/ANT-POST %%%%%%%%%%%%%
cog_roi_ind = [ 5 5 5 5  5  5 ];
new_roi_ind = [ 6 7 8 15 16 17 ];
[ wmp_dta_col(cog_roi_ind)' new_dta_col(new_roi_ind)'];

org_roi_ind = [ 6 39 ];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)
    subplot(3,3,iP)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(new_dta_col{new_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

for iP = 1:numel(org_roi_ind)
    subplot(3,2,iP+4)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(wmp_dta_col{org_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

tightfig();
print(gcf,[ out_dir '/' 'new_ROIs_ANT_post.png'],'-dpng');
close all

% New ROIs/BNT-PRE %%%%%%%%%%%%%
cog_roi_ind = [ 1 1 1 1  1  1 ];
new_roi_ind = [ 6 7 8 15 16 17 ];
[ wmp_dta_col(cog_roi_ind)' new_dta_col(new_roi_ind)'];

org_roi_ind = [ 6 39 ];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)
    subplot(3,3,iP)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(new_dta_col{new_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

for iP = 1:numel(org_roi_ind)
    subplot(3,2,iP+4)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(wmp_dta_col{org_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

tightfig();
print(gcf,[ out_dir '/' 'new_ROIs_BNT_pre.png'],'-dpng');
close all

% New ROIs/ANT-PRE %%%%%%%%%%%%%
cog_roi_ind = [ 2 2 2 2  2  2 ];
new_roi_ind = [ 6 7 8 15 16 17 ];
[ wmp_dta_col(cog_roi_ind)' new_dta_col(new_roi_ind)'];

org_roi_ind = [ 6 39 ];

figure('Position',[0 0 1080 1080],'Units','pixels','Visible','off')
for iP = 1:numel(new_roi_ind)
    subplot(3,3,iP)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(new_dta(:,new_roi_ind(iP))) );
    title( sprintf( [ mmil_spec_char(new_dta_col{new_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

for iP = 1:numel(org_roi_ind)
    subplot(3,2,iP+4)
    [rvl_hld, pvl_hld] = corr(cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))), 'type', 'spearman', 'row', 'complete');
    scatter( cog_dta( grp.tle_post_3T_ATLonly_left, cog_roi_ind(iP) ), cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,org_roi_ind(iP))) );
    title( sprintf([ mmil_spec_char(wmp_dta_col{org_roi_ind(iP)},{'_'},{' '}) '\n r=' num2str(roundsd(rvl_hld,2)) ' p=' num2str(roundsd(pvl_hld,2))] ))
end

tightfig();
print(gcf,[ out_dir '/' 'new_ROIs_ANT_pre.png'],'-dpng');
close all

