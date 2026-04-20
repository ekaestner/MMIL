for iDT = 1:numel(dta_typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/'])

    %%
    fcfg = [];
    fcfg.dta_loc = [ new_dta_dir '/' '03_covariates' '/' dta_typ_nme{iDT} '/' '00_' dta_typ_nme{iDT} '_' 'covariates' '_' dte_str '.csv'];
    fcfg.dta_col = 2;
    [ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

    %%
    ste_var = 'Site';           ste_var_col = strcmpi(cov_col,ste_var);
    con_var = 'Diagnosis';      con_var_col = strcmpi(cov_col,con_var);
    loc_var = 'Localization';   loc_var_col = strcmpi(cov_col,loc_var);
    lat_var = 'Lateralization'; lat_var_col = strcmpi(cov_col,lat_var);
    dta_set_var = 'dataset';    dta_set_col = strcmpi(cov_col,dta_set_var);
    ste_var     = 'Site';       ste_var_col = strcmpi(cov_col,ste_var);
    mts_var     = 'MTS';        mts_var_col = strcmpi(cov_col,mts_var);
    foc_var = 'Focal';   foc_var_col = strcmpi(cov_col,foc_var);


    %%
    epd_sbj = strcmpi(cov_dta(:,con_var_col),'EPD');
    con_sbj = strcmpi(cov_dta(:,con_var_col),'HC');

    tle_sbj = strcmpi(cov_dta(:,loc_var_col),'temporal');

    lft_sbj = strcmpi(cov_dta(:,lat_var_col),'left');
    rgh_sbj = strcmpi(cov_dta(:,lat_var_col),'right');

    mri_pls_sbj = strcmpi(cov_dta(:,mts_var_col),'HS_mri') | strcmpi(cov_dta(:,mts_var_col),'HS_pathology');
    mri_neg_sbj = strcmpi(cov_dta(:,mts_var_col),'noHS');

    gge_sbj = strcmpi(cov_dta(:,foc_var_col),'generalized');
    fle_sbj = strcmpi(cov_dta(:,foc_var_col),'focal') & ( strcmpi(cov_dta(:,loc_var_col),'ex-tle-frontal') | ...
        strcmpi(cov_dta(:,loc_var_col),'frontal'));

    dta_set_nme = unique(cov_dta(:,dta_set_col));

    ste_nme = unique(cov_dta(:,ste_var_col)); ste_nme(1) = [];

    %% Total
    % Diagnosis
    grp.total.all.diagnosis.EPD = find(epd_sbj);
    grp.total.all.diagnosis.HC  = find(con_sbj);

    % Temporal
    grp.total.all.temporal.TLE = find(epd_sbj & tle_sbj & ( lft_sbj | rgh_sbj));
    grp.total.all.temporal.HC  = find(con_sbj);

    % Diagnosis MRI+
    grp.total.all.temporal_mri_plus.TLE_mri_plus = find(epd_sbj & tle_sbj & ( lft_sbj | rgh_sbj) & mri_pls_sbj);
    grp.total.all.temporal_mri_plus.HC           = find(con_sbj);

    % Lateralization
    grp.total.all.lateralization.left  = find(lft_sbj & tle_sbj );
    grp.total.all.lateralization.right = find(rgh_sbj & tle_sbj );
    grp.total.all.lateralization.HC    = find(con_sbj);

    % MRI+ subjects
    grp.total.all.lateralization_mri_plus.left_mri_pls  = find(lft_sbj & tle_sbj & mri_pls_sbj );
    grp.total.all.lateralization_mri_plus.right_mri_pls = find(rgh_sbj & tle_sbj & mri_pls_sbj );
    grp.total.all.lateralization_mri_plus.HC            = find(con_sbj);

    % Extratemporal_True | non-TLE subjects (GGE & focal ex-TLE)
    grp.total.all.extratemporal_true.GGE = find(epd_sbj & gge_sbj);
    grp.total.all.extratemporal_true.FLE = find(epd_sbj & fle_sbj);
    grp.total.all.extratemporal_true.HC = find(con_sbj);

    % Extratemporal_TLE | non-TLE subjects (GGE & focal ex-TLE)
    grp.total.all.extratemporal_tle.HC = find(con_sbj);
    grp.total.all.extratemporal_tle.TLE = find(epd_sbj & (fle_sbj | gge_sbj));

    %% Dataset
    for iD = 1:numel(dta_set_nme)

        grp.dataset.(dta_set_nme{iD}).diagnosis.EPD = find(epd_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).diagnosis.HC  = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

        % Temporal
        grp.dataset.(dta_set_nme{iD}).temporal.TLE = find(epd_sbj & tle_sbj & ( lft_sbj | rgh_sbj) & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).temporal.HC  = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

        % Lateralization
        grp.dataset.(dta_set_nme{iD}).lateralization.left  = find(lft_sbj & tle_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).lateralization.right = find(rgh_sbj & tle_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).lateralization.HC    = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

        % MRI+ subjects
        grp.dataset.(dta_set_nme{iD}).lateralization_mri_plus.left_mri_pls  = find(lft_sbj & tle_sbj & mri_pls_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).lateralization_mri_plus.right_mri_pls = find(rgh_sbj & tle_sbj & mri_pls_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
        grp.dataset.(dta_set_nme{iD}).lateralization_mri_plus.HC            = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

    end

    %% Site
    for iS = 1:numel(ste_nme)
        try
            grp.site.(ste_nme{iS}).diagnosis.EPD = find(epd_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
            grp.site.(ste_nme{iS}).diagnosis.HC  = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));

            % Temporal
            grp.site.(ste_nme{iS}).temporal.TLE = find(epd_sbj & tle_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
            grp.site.(ste_nme{iS}).temporal.HC  = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));

            % Lateralization
            grp.site.(ste_nme{iS}).lateralization.left  = find(lft_sbj & tle_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
            grp.site.(ste_nme{iS}).lateralization.right = find(rgh_sbj & tle_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
            grp.site.(ste_nme{iS}).lateralization.HC    = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
        catch; end

        % MRI+ subjects
        grp.site.(ste_nme{iS}).lateralization_mri_plus.left_mri_pls  = find(lft_sbj & tle_sbj & mri_pls_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
        grp.site.(ste_nme{iS}).lateralization_mri_plus.right_mri_pls = find(rgh_sbj & tle_sbj & mri_pls_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
        grp.site.(ste_nme{iS}).lateralization_mri_plus.HC            = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));

    end

    %% Save
    save([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/' 'group.mat'], 'grp')

end
