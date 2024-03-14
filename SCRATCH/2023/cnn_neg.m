clear; clc;

dta_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/05_Developement/cnn_neg';

%% Load
fcfg = [];
fcfg.dta_loc = [ dta_fld '/' 'Sz_Outcomes.csv'];
[ sze_dta, sze_sbj, sze_col] = ejk_dta_frm(fcfg);

load([ dta_fld '/' 'demographics.mat']);
load([ dta_fld '/' 'prediction.mat']);

%% Combine in Seizure outcome
cmb_col{end+1} = 'sze_out';

sze_out = nan(numel(cmb_sbj),1);
for iS = 1:numel(sze_sbj)
    ind_sbj = strcmpi(cmb_sbj,sze_sbj{iS});
    if any(ind_sbj)
        sze_out(ind_sbj,1) = sze_dta{iS,strcmpi(sze_col,'Outcome')};
    end
end

cmb_dta(:,end+1) = num2cell(sze_out);

%% Create Groups
sze_val = cell2mat(cmb_dta(:,strcmpi(cmb_col,'sze_out')));
rad_val = cmb_dta(:,strcmpi(cmb_col,'MTS Status'));

grp.rad_out.mri_neg_sze_fre = find( strcmpi(rad_val,'MRI-') & sze_val==1 & ~isnan(sze_val) );
grp.rad_out.mri_neg_sze_bad = find( strcmpi(rad_val,'MRI-') & sze_val~=1 & ~isnan(sze_val) );
grp.rad_out.mri_pos_sze_fre = find( strcmpi(rad_val,'MRI+') & sze_val==1 & ~isnan(sze_val) );
grp.rad_out.mri_pos_sze_bad = find( strcmpi(rad_val,'MRI+') & sze_val~=1 & ~isnan(sze_val) );

%% Check Performance
fprintf('%s\n',[ 'MRI-, Seizure Free     | n: ' num2str(numel(grp.rad_out.mri_neg_sze_fre)) ' | mean: ' num2str(roundsd(mean(prd_scr(grp.rad_out.mri_neg_sze_fre)),3)) ' , median: ' num2str(roundsd(median(prd_scr(grp.rad_out.mri_neg_sze_fre)),3))])
fprintf('%s\n',[ 'MRI-, Seizures Ongoing | n: ' num2str(numel(grp.rad_out.mri_neg_sze_bad)) ' | mean: ' num2str(roundsd(mean(prd_scr(grp.rad_out.mri_neg_sze_bad)),3)) ' , median: ' num2str(roundsd(median(prd_scr(grp.rad_out.mri_neg_sze_bad)),3))])
fprintf('%s\n',[ 'MRI+, Seizure Free     | n: ' num2str(numel(grp.rad_out.mri_pos_sze_fre)) ' | mean: ' num2str(roundsd(mean(prd_scr(grp.rad_out.mri_pos_sze_fre)),3)) ' , median: ' num2str(roundsd(median(prd_scr(grp.rad_out.mri_pos_sze_fre)),3))])
fprintf('%s\n',[ 'MRI+, Seizures Ongoing | n: ' num2str(numel(grp.rad_out.mri_pos_sze_bad)) ' | mean: ' num2str(roundsd(mean(prd_scr(grp.rad_out.mri_pos_sze_bad)),3)) ' , median: ' num2str(roundsd(median(prd_scr(grp.rad_out.mri_pos_sze_bad)),3))])

numel(sze_sbj)

neg_pvl = ranksum( prd_scr(grp.rad_out.mri_neg_sze_fre), prd_scr(grp.rad_out.mri_neg_sze_bad) )
[~, neg_pvl] = ttest2( prd_scr(grp.rad_out.mri_neg_sze_fre), prd_scr(grp.rad_out.mri_neg_sze_bad) )

%% Save out
cell2csv( [ dta_fld '/' 'mri_negative_and_seizure_free_subjects.csv' ], [ cmb_sbj(grp.rad_out.mri_neg_sze_fre) cmb_dta( grp.rad_out.mri_neg_sze_fre, ismember(cmb_col,{'sze_out' 'MTS Status'})) ]);
