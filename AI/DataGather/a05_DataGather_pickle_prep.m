for iDT = 1:numel(dta_typ_nme)

    load([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/' 'group.mat'], 'grp')

    %% Setup Tasks
    tsk_nme     = { 'HC_vs_EPD' 'LTLE_vs_RTLE' 'HC_vs_LTLE_vs_RTLE' ...
                    'HC_vs_EPD_MRI_Plus' 'LTLE_vs_RTLE_MRI_Plus' 'HC_vs_LTLE_vs_RTLE_MRI_Plus' ...
                    'Extratemporal_True' 'Extratemporal_TLE' };
    tsk_fst_lvl = { 'total'     'total'        'total'              ...
                    'total'              'total'                 'total' ...
                    'total'                 'total' };
    tsk_scd_lvl = { 'all'       'all'           'all'               ...
                    'all'                'all'                   'all' ...
                    'all'                 'all' };
    tsk_thr_lvl = { 'temporal'  'lateralization' 'lateralization'   ...
                    'temporal_mri_plus' 'lateralization_mri_plus' 'lateralization_mri_plus'...
                    'extratemporal_true' 'extratemporal_tle'};
    tsk_grp     = { { 'HC' 'TLE' } {'left' 'right'} { 'HC' 'left' 'right' } ...
                    { 'HC' 'TLE_mri_plus' } { 'left_mri_pls' 'right_mri_pls' } { 'HC' 'left_mri_pls' 'right_mri_pls' } ...
                    { 'HC' 'GGE' 'FLE' } { 'HC' 'TLE' } };

    cnt_col = { 'Age' };

    %%
    fcfg = [];
    fcfg.dta_loc = [ new_dta_dir '/' '03_covariates' '/' dta_typ_nme{iDT} '/' '00_' dta_typ_nme{iDT} '_' 'covariates' '_' dte_str '.csv'];
    fcfg.dta_col = 2;
    [ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

    loc_sbj = mmil_readtext([ new_dta_dir '/' '02_subjects' '/' dta_typ_nme{iDT} '/' '00_' dta_typ_nme{iDT}  '_' 'subjects' '_' dte_str '.csv']);

    %%
    pkl_col = cell( [ {'sbj_nme'} {'file'} strcat(tsk_nme,'_numeric') strcat(tsk_nme,'_string') ] );
    pkl_dta = cell( size(loc_sbj,1), numel(tsk_nme)*2+1);

    pkl_dta(:,1) = loc_sbj(:,2);
    for iPK = 1:numel(tsk_nme)

        % Numeric setup
        pkl_dta(:,iPK+1) = {NaN};
        for iG = 1:numel(tsk_grp{iPK})
            pkl_dta( grp.(tsk_fst_lvl{iPK}).(tsk_scd_lvl{iPK}).(tsk_thr_lvl{iPK}).(tsk_grp{iPK}{iG}) , iPK+1 ) = {iG};
        end

        % String setup
        pkl_dta(:,iPK+numel(tsk_nme)+1) = {''};
        for iG = 1:numel(tsk_grp{iPK})
            pkl_dta( grp.(tsk_fst_lvl{iPK}).(tsk_scd_lvl{iPK}).(tsk_thr_lvl{iPK}).(tsk_grp{iPK}{iG}) , iPK+numel(tsk_nme)+1 ) = tsk_grp{iPK}(iG);
        end

    end

    %% Alphabetize
    [ ~, new_row ] = sort(lower(cov_sbj));
    pkl_dta = pkl_dta(new_row,:);
    cov_sbj = cov_sbj(new_row,:);
    cov_dta = cov_dta(new_row,:);

    %% Add Continuous Variables
    for iC = 1:numel(cnt_col)
        pkl_dta = [ pkl_dta cov_dta(:,strcmpi(cov_col,cnt_col{iC}))];
    end

    %%
    cell2csv([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/' 'pickle_prep_' dte_str '.csv'], [ pkl_col cnt_col ; cov_sbj pkl_dta ])

    %% Remove duplicates
    [~, uniqueIdx] = unique(cov_sbj);
    cov_sbj = cov_sbj(uniqueIdx,:);
    pkl_dta = pkl_dta(uniqueIdx,:);

    [ ~, new_row ] = sort(lower(cov_sbj));
    pkl_dta = pkl_dta(new_row,:);
    cov_sbj = cov_sbj(new_row,:);

    cell2csv([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/' 'pickle_prep_' dte_str '_no_duplicates.csv'], [ pkl_col cnt_col ; cov_sbj pkl_dta ])

end