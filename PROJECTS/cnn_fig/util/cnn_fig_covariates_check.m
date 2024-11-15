%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'];
[ bth_dta, bth_sbj, bth_col ] = ejk_dta_frm(fcfg);

bth_ste = bth_dta(:,strcmpi(bth_col,'Site'));
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

cell2csv([ prj_dir '/' 'new_data' '/' 'covariates' '/' 'site_covariate_check' '_' dte_str '.csv'] ,[ {'Site'} cov_out_col ; [ {'Total'} ; tot_ste] cov_tbl]);
