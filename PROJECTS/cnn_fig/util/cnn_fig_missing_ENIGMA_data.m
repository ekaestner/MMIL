
out_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/new_data/enigma_conglom/cat12/investigation';

iD = 1;

%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'];
[ bth_dta, bth_sbj, bth_col ] = ejk_dta_frm(fcfg);

kep_bth_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Imaging_Name')));
bth_dta = bth_dta(kep_bth_ind,:);
bth_sbj = bth_sbj(kep_bth_ind,:);

%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'missing_Imaging' '_' dte_str '.csv'];
[ mss_img_dta, mss_img_sbj, mss_img_col ] = ejk_dta_frm(fcfg);

cov_nme = fieldnames(cde_bok.(bid_nme{iD}));
cov_val = cell(1,numel(cov_nme));
for iCV = 1:numel(cov_nme)
    cov_val{iCV} = cde_bok.(bid_nme{iD}).(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
end

fcfg = [];

fcfg.dta     = mss_img_dta;
fcfg.dta_col = mss_img_col;

fcfg.rcd_nme = cov_nme;
fcfg.rcd_val = cov_val;

fcfg.swt_val = 1;

mss_img_dta_rcd = ejk_recode(fcfg);

%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'missing_Redcap' '_' dte_str '.csv'];
[ mss_cov_dta, mss_cov_sbj, mss_cov_col ] = ejk_dta_frm(fcfg);

%% Initial play
bth_ste = bth_dta(:,strcmpi(bth_col,'Site'));

% Total Numbers = 
% 2455        2208        1034
% 2926        2192        790
% 2838        2168        878
[ numel(bth_sbj) numel(mss_img_sbj) numel(mss_cov_sbj) ]

% Redcap missing imaging = 1574 | 634
sum(cellfun(@isempty,mss_img_dta_rcd(:,strcmpi(mss_img_col,'subjid'))))
sum(~cellfun(@isempty,mss_img_dta_rcd(:,strcmpi(mss_img_col,'subjid'))))

mss_img_non_ste = mss_img_dta_rcd(cellfun(@isempty,mss_img_dta_rcd(:,strcmpi(mss_img_col,'subjid'))),strcmpi(mss_img_col,'site'));
mss_img_yes_ste = mss_img_dta_rcd(~cellfun(@isempty,mss_img_dta_rcd(:,strcmpi(mss_img_col,'subjid'))),strcmpi(mss_img_col,'site'));

% Missing Redcap
mss_cov_ste = mss_cov_dta(:,strcmpi(mss_cov_col,'site'));

%% Make Table of Site
tot_ste = unique( [ bth_ste ; mss_img_yes_ste ; mss_img_non_ste ; mss_cov_ste ] );

num_col{1} = "Both CAT12 & Redcap";
num_col{2} = "Only CAT12";
num_col{3} = "Only Redcap / Has imaging";
num_col{4} = "Only Redcap / No imaging";
num_tbl = cell(numel(tot_ste),4);
for iT = 1:numel(tot_ste)
    
    num_tbl{iT,1} = sum(strcmpi(bth_ste,tot_ste{iT}));
    num_tbl{iT,2} = sum(strcmpi(mss_cov_ste,tot_ste{iT}));
    num_tbl{iT,3} = sum(strcmpi(mss_img_yes_ste,tot_ste{iT}));
    num_tbl{iT,4} = sum(strcmpi(mss_img_non_ste,tot_ste{iT}));

end
cell2csv( [ out_dir '/' 'check_match_table' '_' dte_str '.csv'],[ {'site'} num_col ; tot_ste num_tbl]);

cell2csv( [ out_dir '/' 'check_missing_imaging' '_' dte_str '.csv'], [ {'sbj_nme'} mss_img_col ; mss_img_sbj mss_img_dta_rcd ]);

%% Determine presence of covariates
% tot_ste = unique( bth_ste );
% 
% % [ bth_dta, bth_sbj, bth_col ]
% 
% epd_var = { 'Focal' 'Localization' 'Lateralization' 'MTS' };
% 
% cov_tbl = cell(numel(tot_ste) , numel(bth_col)+2 );
% for iT = 1:numel(tot_ste)
%     cov_tbl{iT,1} = sum(strcmpi(bth_dta(:,strcmpi(bth_col,'Site')),tot_ste{iT}));
%     cov_tbl{iT,2} = sum(strcmpi(bth_dta(strcmpi(bth_dta(:,strcmpi(bth_col,'Site')),tot_ste{iT}),strcmpi(bth_col,'Diagnosis')),'EPD'));
% 
%     for iC = 1:numel(bth_col)
%         if any(ismember(epd_var,bth_col{iC}))
%             cov_tbl{iT,iC+2} = roundsd((sum(~cellfun(@isempty,bth_dta(strcmpi(bth_dta(:,strcmpi(bth_col,'Site')),tot_ste{iT}),iC))) / cov_tbl{iT,2}) * 100,3);
%         else
%             cov_tbl{iT,iC+2} = roundsd((sum(~cellfun(@isempty,bth_dta(strcmpi(bth_dta(:,strcmpi(bth_col,'Site')),tot_ste{iT}),iC))) / cov_tbl{iT,1}) * 100,3);
%         end
%     end
% end
% 
% cell2csv( [ out_dir '/' 'check_site_covariates' '_' dte_str '.csv'], [ {'Site'} {'N-all'} {'N-EPD'} bth_col; tot_ste cov_tbl ] );

%% Investigate Missing Naming
iD = 1;

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'covariates' '/' '/' red_cap_fle.(bid_nme{iD})];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

if strcmpi(bid_nme{iD},'enigma_conglom'); cnn_fig_format_ENIGMA_redcap; end

cov_nme = fieldnames(cde_bok.(bid_nme{iD}));
cov_val = cell(1,numel(cov_nme));
for iCV = 1:numel(cov_nme)
    cov_val{iCV} = cde_bok.(bid_nme{iD}).(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
end

fcfg = [];

fcfg.dta     = red_cap_dta;
fcfg.dta_col = red_cap_col;

fcfg.rcd_nme = cov_nme;
fcfg.rcd_val = cov_val;

fcfg.swt_val = 1;

red_cap_dta_rcd = ejk_recode(fcfg);

tot_ste = unique(red_cap_dta_rcd(:,strcmpi(red_cap_col,'site')));

mss_nme = cell(numel(tot_ste),3);
for iT = 1:numel(tot_ste)
    mss_nme{iT,1} = tot_ste{iT};
    mss_nme{iT,2} = sum(strcmpi(red_cap_dta_rcd(:,strcmpi(red_cap_col,'site')),tot_ste{iT}));
    mss_nme{iT,3} = sum(cellfun(@isempty,red_cap_dta_rcd(strcmpi(red_cap_dta_rcd(:,strcmpi(red_cap_col,'site')),tot_ste{iT}),strcmpi(red_cap_col,'subjid'))));
end
cell2csv( [ out_dir '/' 'check_imaging_name' '_' dte_str '.csv'], mss_nme );

%% Investigate River of Covariates Flow
iD = 1;

%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'];
[ bth_dta, bth_sbj, bth_col ] = ejk_dta_frm(fcfg);

tot_ste = unique( bth_ste );
cov_tbl = cell(numel(tot_ste)+1 , 13 );

%
hss_img_nme_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Imaging_Name')));
mss_img_nme_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Imaging_Name')));

hss_dxe_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Diagnosis')));
mss_dxe_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Diagnosis')));

%
mss_sex_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Sex'))); 
hss_sex_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Sex')));

mss_age_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Age'))); 
hss_age_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Age')));

mss_hnd_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Handedness')));
hss_hnd_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Handedness')));

%
con_ind = strcmpi(bth_dta(:,strcmpi(bth_col,'Diagnosis')),'HC');
epd_ind = strcmpi(bth_dta(:,strcmpi(bth_col,'Diagnosis')),'EPD');

hss_foc_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Focal')));
mss_foc_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Focal')));

mss_loc_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Localization')));
hss_loc_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Localization')));

tle_ind = strcmpi(bth_dta(:,strcmpi(bth_col,'Localization')),'temporal') | strcmpi(bth_dta(:,strcmpi(bth_col,'Localization')),'temporal_plus');
ext_ind = strcmpi(bth_dta(:,strcmpi(bth_col,'Localization')),'ex-tle-frontal');

mss_lat_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Lateralization')));
hss_lat_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Lateralization')));

mss_mts_ind = cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'MTS')));
hss_mts_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'MTS')));

% CAT12 -> "subjid" -> Redcap
cov_tbl{1,1} = sprintf('%i / %i',sum(hss_img_nme_ind),sum(mss_img_nme_ind)); 
    cov_out_col{1} = 'Has Redcap Match (y/n)'; 
cov_tbl{1,2} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind),sum(hss_img_nme_ind & mss_dxe_ind));
    cov_out_col{2} = 'Has Diagnosis (y/n)'; 

cov_tbl{1,3} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & hss_sex_ind), sum(hss_img_nme_ind & hss_dxe_ind & mss_sex_ind));
    cov_out_col{3} = 'Has Sex (y/n)'; 
cov_tbl{1,4} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & hss_age_ind), sum(hss_img_nme_ind & hss_dxe_ind & mss_age_ind));
    cov_out_col{4} = 'Has Age (y/n)'; 
cov_tbl{1,5} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & hss_hnd_ind), sum(hss_img_nme_ind & hss_dxe_ind & mss_hnd_ind));
    cov_out_col{5} = 'Has Age (y/n)'; 

cov_tbl{1,6} = sprintf('%i',sum(hss_img_nme_ind & hss_dxe_ind & con_ind));
    cov_out_col{6} = 'HC (N)';
cov_tbl{1,7} = sprintf('%i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind));
    cov_out_col{7} = 'EPD (N)';
cov_tbl{1,8} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind),sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & mss_foc_ind));
    cov_out_col{8} = 'Has Focal (y/n)';
cov_tbl{1,9}  = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind),sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & mss_loc_ind));
    cov_out_col{9} = 'Has Localization (y/n)';
cov_tbl{1,10} = sprintf('%i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & ext_ind));
    cov_out_col{10} = 'ex-TLE (N)';
cov_tbl{1,11} = sprintf('%i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind));
    cov_out_col{11} = 'TLE (N)';
cov_tbl{1,12} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind),sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & mss_lat_ind));
    cov_out_col{12} = 'Has Lateralization (y/n)';
cov_tbl{1,13} = sprintf('%i / %i',sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind & hss_mts_ind),sum(hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind & mss_mts_ind));
    cov_out_col{13} = 'Has MTS (y/n)';  
  
for iT = 1:numel(tot_ste)

    ste_ind = strcmpi(bth_dta(:,strcmpi(bth_col,'Site')),tot_ste{iT});

    cov_tbl{iT+1,1} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind),sum(ste_ind & mss_img_nme_ind));
    cov_tbl{iT+1,2} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind),sum(ste_ind & hss_img_nme_ind & mss_dxe_ind));

    cov_tbl{iT+1,3} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & hss_sex_ind), sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & mss_sex_ind));
    cov_tbl{iT+1,4} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & hss_age_ind), sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & mss_age_ind));
    cov_tbl{iT+1,5} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & hss_hnd_ind), sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & mss_hnd_ind));


    cov_tbl{iT+1,6} = sprintf('%i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & con_ind));
    cov_tbl{iT+1,7} = sprintf('%i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind));
    cov_tbl{iT+1,8} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind),sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & mss_foc_ind));
    cov_tbl{iT+1,9}  = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind),sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & mss_loc_ind));
    cov_tbl{iT+1,10} = sprintf('%i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & ext_ind));
    cov_tbl{iT+1,11} = sprintf('%i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind));
    cov_tbl{iT+1,12} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind),sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & mss_lat_ind));
    cov_tbl{iT+1,13} = sprintf('%i / %i',sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind & hss_mts_ind),sum(ste_ind & hss_img_nme_ind & hss_dxe_ind & epd_ind & hss_foc_ind & hss_loc_ind & tle_ind & hss_lat_ind & mss_mts_ind));

end

cell2csv([ out_dir '/' 'check_site_covariates' '_' dte_str '.csv'] ,[ {'Site'} cov_out_col ; [ {'Total'} ; tot_ste] cov_tbl]);

%% Investigate Exact Matches
% red_cap_dta_rcd, red_cap_col, red_cap_sbj

iD = 1;

% %%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'missing_Redcap' '_' dte_str '.csv'];
[ mss_cov_dta, mss_cov_sbj, mss_cov_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' bid_nme{iD} '_' 'covariates.csv'];
[ bth_dta, bth_sbj, bth_col ] = ejk_dta_frm(fcfg);

kep_bth_ind = ~cellfun(@isempty,bth_dta(:,strcmpi(bth_col,'Imaging_Name')));
bth_dta = bth_dta(kep_bth_ind,:);
bth_sbj = bth_sbj(kep_bth_ind,:);

ste_sbj = [ mss_cov_sbj mss_cov_dta(:, strcmpi(mss_cov_col,'site')) ; bth_sbj bth_dta(:,strcmpi(bth_col,'Site'))];

% %%%%%%%%%%
img_fle = dir([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/']);
    img_fle = { img_fle(:).name};
    img_fle = img_fle(string_find(img_fle,'sub-'));
    img_fle = strrep(img_fle,'__cat12_mwp1.nii',"");
    img_fle = strrep(img_fle,'_cat12_mwp1.nii',"");
    img_fle = cellfun(@(x) strsplit(x,'_ses-'),img_fle,'uni',0);
    img_fle = cellfun(@(x) x{1}, img_fle,'uni',0);

% %%%%%%%%%%
sbj_nme_out = cell(numel(img_fle),3);
for iS = 1:numel(img_fle)

    sbj_nme_out{iS,1} = img_fle{iS};

    ste_ind = strcmpi(ste_sbj(:,1), img_fle{iS});
    if any(ste_ind)
        sbj_nme_out{iS,2} = ste_sbj{ste_ind,2};
    else
        sbj_nme_out{iS,2} = '';
    end

    mtc_ind = strcmpi(red_cap_dta_rcd(:,strcmpi(red_cap_col,'subjid')), img_fle{iS});
    if sum(mtc_ind)==1
        sbj_nme_out{iS,3} = red_cap_sbj{mtc_ind};
    end

end

fprintf('\nTotal subjid match: %i \n',sum(~cellfun(@isempty,sbj_nme_out(:,3))))
fprintf('\nMissing subjid match: %i \n',sum(cellfun(@isempty,sbj_nme_out(:,3))))
fprintf('\nMissing site: %i \n',sum(cellfun(@isempty,sbj_nme_out(:,2))))

cell2csv([ out_dir '/' 'strict_match' '_' dte_str '.csv'],sbj_nme_out)

% %%%%%%%%%%
tot_ste = unique(sbj_nme_out(:,2));

sbj_nme_tbl = cell(numel(tot_ste),4);
for iT = 1:numel(tot_ste)
    ste_ind = strcmpi(sbj_nme_out(:,2),tot_ste{iT});
    sbj_nme_tbl{iT,1} = tot_ste{iT};
    sbj_nme_tbl{iT,2} = sum(ste_ind);
    sbj_nme_tbl{iT,3} = sum(~cellfun(@isempty,sbj_nme_out(ste_ind,3)));
    sbj_nme_tbl{iT,4} = sum(cellfun(@isempty,sbj_nme_out(ste_ind,3)));
end
sbj_nme_tbl{1,1} = 'UCSD_ENIGMA';
cell2csv([ out_dir '/' 'strict_match_check' '_' dte_str '.csv'],[ {'Site'} {'Total Scans'} {'Scans with matching subjid'} {'Scans without matching subjid'} ; sbj_nme_tbl])

%% Investigate Imaging in Redcap
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'covariates' '/' 'ENIGMAEpilepsy-Subjidcheck_DATA_LABELS_2024-08-28_0925.csv'];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

% %%%%%%%%%%
tot_ste = unique(red_cap_dta(:,strcmpi(red_cap_col,'Institution Name')));

sbj_nme_tbl = cell(numel(tot_ste),4);
for iT = 1:numel(tot_ste)
    ste_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'Institution Name')),tot_ste{iT});
    sbj_nme_tbl{iT,1} = tot_ste{iT};
    sbj_nme_tbl{iT,2} = sum(ste_ind);
    sbj_nme_tbl{iT,3} = sum(~cellfun(@isempty,red_cap_dta(ste_ind,strcmpi(red_cap_col,'BIDS Subject ID'))));
    sbj_nme_tbl{iT,4} = sum(cellfun(@isempty,red_cap_dta(ste_ind,strcmpi(red_cap_col,'BIDS Subject ID'))));
end

cell2csv([ out_dir '/' 'Subjidcheck' '_' dte_str '.csv'],[ {'Site'} {'Total Scans'} {'Scans with matching subjid'} {'Scans without matching subjid'} ; sbj_nme_tbl])

%% Mass cross-check
iD = 1;

% %%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'covariates' '/' 'ENIGMAEpilepsy-Subjidcheck_DATA_LABELS_2024-08-28_0925.csv'];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);
cnn_fig_format_ENIGMA_redcap_subjid_check

% %%%%%%%%%%
img_fle_org = dir([ bid_dir '/' bid_nme{iD} '' '/' ]);
    img_fle_org = { img_fle_org(:).name};
    img_fle_org = img_fle_org(string_find(img_fle_org,'sub-'));

img_fle_rnm = dir([ bid_dir '/' bid_nme{iD} '_' 'rename' '/' ]);
    img_fle_rnm = { img_fle_rnm(:).name};
    img_fle_rnm = img_fle_rnm(string_find(img_fle_rnm,'sub-'));

tot_img_fle = unique([ img_fle_org img_fle_rnm]);

% %%%%%%%%%%
sbj_nme_bid = cell(numel(tot_img_fle),7); % Redcap ID, BIDS ID, Original BIDS ID, Institution, Original Match, Rename Match, Classification
red_cap_bid_fnd = [];
for iS = 1:numel(tot_img_fle)
    
    % %%%%%%
    if any(strcmpi(img_fle_org,tot_img_fle{iS}))
        mtc_ind_org = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'BIDS Subject ID')),tot_img_fle{iS});
    else
        mtc_ind_org = zeros(1,numel(red_cap_dta(:,strcmpi(red_cap_col,'BIDS Subject ID'))));
    end
    if any(strcmpi(img_fle_rnm,tot_img_fle{iS}))
        mtc_ind_rnm = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'BIDS Subject ID')),tot_img_fle{iS});
    else
        mtc_ind_rnm = zeros(1,numel(red_cap_dta(:,strcmpi(red_cap_col,'BIDS Subject ID'))));
    end

    % %%%%%%
    sbj_nme_bid{iS,2} = tot_img_fle{iS};
    if (any(mtc_ind_org) || any(mtc_ind_rnm)) && ~(sum(mtc_ind_org)>1 || sum(mtc_ind_rnm)>1)

        if any(mtc_ind_org); mtc_ind_org_fnd = find(mtc_ind_org); else mtc_ind_org_fnd = []; end
        if any(mtc_ind_rnm); mtc_ind_rnm_fnd = find(mtc_ind_rnm); else mtc_ind_rnm_fnd = []; end

        if any(mtc_ind_org) && any(mtc_ind_rnm) && (mtc_ind_rnm_fnd ~= mtc_ind_rnm_fnd); error('Mismatch names'); end

        mtc_fnd = unique([ mtc_ind_org_fnd mtc_ind_rnm_fnd]);
        red_cap_bid_fnd = [ red_cap_bid_fnd mtc_fnd];

        sbj_nme_bid{iS,1} = red_cap_sbj{mtc_fnd};
        sbj_nme_bid{iS,3} = red_cap_dta{mtc_fnd,strcmpi(red_cap_col,'Subject id given by site (or subjid from previous redcap)')};
        sbj_nme_bid{iS,4} = red_cap_dta{mtc_fnd,strcmpi(red_cap_col,'Institution Name')};

        if any(mtc_ind_org); sbj_nme_bid{iS,5} = 'yes'; else sbj_nme_bid{iS,5} = 'no'; end
        if any(mtc_ind_rnm); sbj_nme_bid{iS,6} = 'yes'; else sbj_nme_bid{iS,6} = 'no'; end
       
    elseif sum(mtc_ind_org)>1 || sum(mtc_ind_rnm)>1

        sbj_nme_bid{iS,5} = 'Naming Issue';
        sbj_nme_bid{iS,6} = 'Naming Issue';

    else

        sbj_nme_bid{iS,5} = '';
        sbj_nme_bid{iS,6} = '';

    end

    % %%%%%%
    if strcmpi( sbj_nme_bid{iS,5}, 'yes' ) && strcmpi(  sbj_nme_bid{iS,6}, 'yes'  )
        sbj_nme_bid{iS,7} = 'Redcap_BIDS_match';
    elseif strcmpi( sbj_nme_bid{iS,5}, 'yes' ) && strcmpi(  sbj_nme_bid{iS,6}, 'no'  )
        sbj_nme_bid{iS,7} = 'Redcap_BIDS_match_ORIGINAL_ONLY';
    elseif strcmpi( sbj_nme_bid{iS,5}, 'no' ) && strcmpi(  sbj_nme_bid{iS,6}, 'yes'  )
        sbj_nme_bid{iS,7} = 'Redcap_BIDS_match_RENAME_ONLY';
    elseif strcmpi( sbj_nme_bid{iS,5}, '' ) && strcmpi(  sbj_nme_bid{iS,6}, ''  )
        sbj_nme_bid{iS,7} = 'BIDS_Missing_Redcap';
    end

end

% %%%%%%%%%%
% red_cap_sbj_oly(red_cap_bid_fnd,:) = [];
% red_cap_dta_oly(red_cap_bid_fnd,:) = [];
% 
% sbj_nme_bid = [ sbj_nme_bid ; red_cap_sbj_oly ...
%                               red_cap_dta_oly(:,strcmpi(red_cap_col,'BIDS Subject ID')) ...
%                               red_cap_dta_oly(:,strcmpi(red_cap_col,'Subject id given by site (or subjid from previous redcap)')) ...
%                               red_cap_dta_oly(:,strcmpi(red_cap_col,'Institution Name')) ...
%                               repmat({''},numel(red_cap_sbj_oly),1) ...
%                               repmat({''},numel(red_cap_sbj_oly),1) ...
%                               repmat({'Redcap_Missing_BIDS'},numel(red_cap_sbj_oly),1) ];

% %%%%%%%%%%
prb_sbj = sbj_nme_bid;
prb_sbj(strcmpi(prb_sbj(:,7),'Redcap_BIDS_match'),:) = [];
prb_sbj(strcmpi(prb_sbj(:,7),'Redcap_BIDS_match_RENAME_ONLY'),:) = [];
cell2csv([ out_dir '/' 'BIDS_match_to_Redcap_Problem_Subjects' '_' dte_str '.csv'], prb_sbj)

% %%%%%%%%%%
cell2csv([ out_dir '/' 'BIDS_match_to_Redcap' '_' dte_str '.csv'], sbj_nme_bid)
fprintf('\n\nRedcap & BIDS match: %i\n',sum(strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match')))
fprintf('Redcap & BIDS match (Rename Only): %i\n',sum(strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match_RENAME_ONLY')))
fprintf('Redcap & BIDS match (Original Only): %i\n',sum(strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match_ORIGINAL_ONLY')))
fprintf('BIDS, missing Redcap: %i\n',sum(strcmpi(sbj_nme_bid(:,7),'BIDS_Missing_Redcap')))
fprintf('Redcap, missing BIDS: %i\n',sum(strcmpi(sbj_nme_bid(:,7),'Redcap_Missing_BIDS')))

% %%%%%%%%%%
prb_sbj_exp = mmil_readtext([ out_dir '/' 'BIDS_match_to_Redcap_Problem_Subjects_2024_08_21_ejk.csv']);
prb_typ = unique(prb_sbj_exp(:,8));

prb_tbl = cell(numel(prb_typ),2);
for iP = 1:numel(prb_typ)
    prb_tbl{iP,1} = prb_typ{iP};
    prb_tbl{iP,2} = sum(strcmpi(prb_sbj_exp(:,8),prb_typ{iP}));
end
prb_tbl

ste_nme = unique(red_cap_dta(:,strcmpi(red_cap_col,'Institution Name')));
ste_tbl = cell(numel(ste_nme),3);
for iST = 1:numel(ste_nme)
    ste_tbl{iST,1} = ste_nme{iST};
    ste_tbl{iST,2} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & (strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match') | strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match_RENAME_ONLY') ));
    ste_tbl{iST,3} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,7),'Redcap_Missing_BIDS'));
    ste_tbl{iST,4} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,7),'Redcap_BIDS_match_ORIGINAL_ONLY'));
end
ste_tbl

% %%%%%%%%%%
% img_fle_rnm = mmil_readtext([ out_dir '/' 'BIDS_match_to_Redcap_2024_08_19_renamed.csv' ]);
% img_fle_org = mmil_readtext([ out_dir '/' 'BIDS_match_to_Redcap_2024_08_19_original.csv' ]);
% 
% img_fle_rnm_mtc = img_fle_rnm(~cellfun(@isempty,img_fle_rnm(:,2)),1);
% img_fle_org_mtc = img_fle_org(~cellfun(@isempty,img_fle_org(:,2)),1);
% 
% numel(img_fle_rnm(:,1))
% numel(img_fle_org(:,1))
% 
% cell2csv([ out_dir '/' 'renamed_and_original_difference.csv'],setxor( img_fle_rnm(:,1), img_fle_org(:,1)))
% 
% setxor(img_fle_rnm_mtc,img_fle_org_mtc)
