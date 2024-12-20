

out_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/new_data/enigma_conglom/cat12/investigation';

eng_nme = { 'enigma_conglom' 'enigma_conglom_rename' };

%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{1} '/' 'covariates' '/' 'ENIGMAEpilepsy-Subjidcheck_DATA_LABELS_2024-08-28_0925.csv'];
[ red_sbj_dta, red_sbj_sbj, red_sbj_col ] = ejk_dta_frm(fcfg);
cnn_fig_format_ENIGMA_redcap_subjid_check

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{1} '/' 'covariates' '/' 'ENIGMAEpilepsy-CNN_DATA_2024-08-28_0926.csv'];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);
cnn_fig_format_ENIGMA_redcap

cov_nme = fieldnames(cde_bok.(bid_nme{1}));
cov_val = cell(1,numel(cov_nme));
for iCV = 1:numel(cov_nme)
    cov_val{iCV} = cde_bok.(bid_nme{1}).(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
end

fcfg = [];

fcfg.dta     = red_cap_dta;
fcfg.dta_col = red_cap_col;

fcfg.rcd_nme = cov_nme;
fcfg.rcd_val = cov_val;

fcfg.swt_val = 1;

red_cap_dta = ejk_recode(fcfg);

%% BIDS/Redcap Match
for iD = 1:numel(eng_nme)

    % %%%%
    img_fle = dir([ bid_dir '/' eng_nme{iD} '' '/' ]);
    img_fle = { img_fle(:).name};
    img_fle = img_fle(string_find(img_fle,'sub-'));

    % %%%%
    sbj_nme_bid = cell(numel(img_fle),6); % Redcap ID, BIDS ID, Original BIDS ID, Institution, Match, Classification
    red_sbj_bid_fnd = [];

    for iS = 1:numel(img_fle)

        % %%%%%%
        mtc_ind = strcmpi(red_sbj_dta(:,strcmpi(red_sbj_col,'BIDS Subject ID')),img_fle{iS});

        % %%%%%%
        sbj_nme_bid{iS,2} = img_fle{iS};
        if any(mtc_ind) && ~(sum(mtc_ind)>1)

            mtc_ind_fnd = find(mtc_ind);

            red_sbj_bid_fnd = [ red_sbj_bid_fnd ; mtc_fnd];

            sbj_nme_bid{iS,1} = red_sbj_sbj{mtc_ind_fnd};
            sbj_nme_bid{iS,3} = red_sbj_dta{mtc_ind_fnd,strcmpi(red_sbj_col,'Subject id given by site (or subjid from previous redcap)')};
            sbj_nme_bid{iS,4} = red_sbj_dta{mtc_ind_fnd,strcmpi(red_sbj_col,'Institution Name')};

            sbj_nme_bid{iS,5} = 'yes';

        elseif sum(mtc_ind)>1

            red_sbj_bid_fnd = [ red_sbj_bid_fnd ; mtc_fnd];
            sbj_nme_bid{iS,5} = 'no';
            sbj_nme_bid{iS,6} = 'Naming Issue';

        else

            sbj_nme_bid{iS,5} = '';
            sbj_nme_bid{iS,6} = 'No Redcap Match';

        end

        % %%%%%%
        if strcmpi( sbj_nme_bid{iS,5}, 'yes' )
            sbj_nme_bid{iS,6} = 'Redcap_BIDS_match';
        elseif strcmpi( sbj_nme_bid{iS,5}, '' )
            sbj_nme_bid{iS,6} = 'BIDS_Missing_Redcap';
        end

    end

    % %%%%%%%%%%
    red_cap_sbj_oly = red_sbj_sbj;
    red_cap_dta_oly = red_sbj_dta;

    red_cap_sbj_oly(red_cap_bid_fnd,:) = [];
    red_cap_dta_oly(red_cap_bid_fnd,:) = [];
    
    sbj_nme_bid = [ sbj_nme_bid ; red_cap_sbj_oly ...
                                  red_cap_dta_oly(:,strcmpi(red_sbj_col,'BIDS Subject ID')) ...
                                  red_cap_dta_oly(:,strcmpi(red_sbj_col,'Subject id given by site (or subjid from previous redcap)')) ...
                                  red_cap_dta_oly(:,strcmpi(red_sbj_col,'Institution Name')) ...
                                  repmat({''},numel(red_cap_sbj_oly),1) ...
                                  repmat({'Redcap Missing BIDS'},numel(red_cap_sbj_oly),1) ];

    % %%%%%%%%%%
    ste_nme = unique(red_sbj_dta(:,strcmpi(red_sbj_col,'Institution Name')));
    ste_tbl = cell(numel(ste_nme),5);
    for iST = 1:numel(ste_nme)
        ste_tbl{iST,1} = ste_nme{iST};
        ste_tbl{iST,2} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,6),'Redcap_BIDS_match'));
        ste_tbl{iST,3} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,6),'Redcap Missing BIDS'));
        ste_tbl{iST,4} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,6),'Naming Issue'));
        ste_tbl{iST,5} = sum( strcmpi(sbj_nme_bid(:,4),ste_nme{iST}) & strcmpi(sbj_nme_bid(:,6),'BIDS_Missing_Redcap'));
    end

    tot_tbl{1,1} = 'Total';
    tot_tbl{1,2} = sum( strcmpi(sbj_nme_bid(:,6),'Redcap_BIDS_match'));
    tot_tbl{1,3} = sum( strcmpi(sbj_nme_bid(:,6),'Redcap Missing BIDS'));
    tot_tbl{1,4} = sum( strcmpi(sbj_nme_bid(:,6),'Naming Issue'));
    tot_tbl{1,5} = sum( strcmpi(sbj_nme_bid(:,6),'BIDS_Missing_Redcap'));

    % %%%%%%%%%%
    cell2csv([ out_dir '/' 'BIDS_match_to_Redcap' '_' eng_nme{iD} '_' dte_str '_' 'subjects.csv'], [ [ {'Redcap ID'} {'BIDS ID'} {'Original BIDS ID'} {'Institution'} {'Match'} {'Classification'} ] ; sbj_nme_bid ]);
    cell2csv([ out_dir '/' 'BIDS_match_to_Redcap' '_' eng_nme{iD} '_' dte_str '.csv'], [ {'Site'} {'Matching'} {'Redcap Only'} {'Naming Issue'} {'BIDS only'} ; tot_tbl ; ste_tbl ])

end

%% Covariates Check: Table 1
for iD = 1:numel(eng_nme)

    %
    bth_ste = red_cap_dta(:,strcmpi(red_cap_col,'site'));
    tot_ste = unique( bth_ste );
    cov_tbl = cell(numel(tot_ste)+1 , 13 );

    %
    hss_img_nme_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'subjid')));
    mss_img_nme_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'subjid')));

    hss_dxe_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'dx')));
    mss_dxe_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'dx')));

    %
    mss_sex_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'sex')));
    hss_sex_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'sex')));

    mss_age_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'age')));
    hss_age_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'age')));

    mss_hnd_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'handedness')));
    hss_hnd_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'handedness')));

    %
    con_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'dx')),'HC');
    epd_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'dx')),'EPD');

    hss_foc_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'dx_epi_type')));
    mss_foc_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'dx_epi_type')));

    mss_loc_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'localization')));
    hss_loc_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'localization')));

    tle_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'localization')),'temporal') | strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'localization')),'temporal_plus');
    ext_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'localization')),'ex-tle-frontal');

    mss_lat_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'lateralization')));
    hss_lat_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'lateralization')));

    mss_mts_ind = cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'mts')));
    hss_mts_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'mts')));

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
    cov_out_col{5} = 'Has Handedness (y/n)';

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

        ste_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'Site')),tot_ste{iT});

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

    cell2csv([ out_dir '/' 'Data_flow' '_' eng_nme{iD} '_' dte_str '_' '.csv'] ,[ {'Site'} cov_out_col ; [ {'Total'} ; tot_ste] cov_tbl]);

end
