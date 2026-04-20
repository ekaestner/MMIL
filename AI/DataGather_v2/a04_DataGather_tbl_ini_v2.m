for iDT = 1:numel(dta_typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '03_Datasets' '/' typ_nme{iDT} '/'])

    %% Load Data & Combine
    tot_cov_dta = cell(0); tot_img_dta = cell(0);
    tot_cov_sbj = cell(0);
    grp_idn = nan(0);
    grp_idn_nme = cell(0);
    for iD = 1:numel(bid_nme)
        if exist([ new_dta_dir '/' '02_covariates' '/' typ_nme{iDT} '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'])

            % Load Imaging Data
            img_fle_pth = [ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' bid_nme{iD} '_' typ_nme{iDT} '_' 'file_characterization_' dte_str '.csv'];
            
            opts = detectImportOptions(img_fle_pth, 'Delimiter', ',', 'FileType', 'text');
            opts.VariableNamesLine = 0;      % Treat headers as data to ensure everything is read
            opts.DataLines = [1 Inf];
            opts = setvartype(opts, 'char'); % Strictly force all columns to be 'char'
            
            T = readtable(img_fle_pth, opts);
            ses_dta = table2cell(T);         % Convert back to cell array structure

            % Load Covariate Data
            fcfg = [];
            fcfg.dta_loc = [ new_dta_dir '/' '02_covariates' '/' typ_nme{iDT} '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'];
            fcfg.dta_col = 2;
            [ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

            % Put together
            tot_cov_dta = [ tot_cov_dta ; cov_dta];
            tot_cov_sbj = [ tot_cov_sbj ; cov_sbj];

            tot_img_dta = [ tot_img_dta ; ses_dta];

            grp_idn = [ grp_idn ; ones(numel(cov_sbj),1)*iD];
            grp_idn_nme = [ grp_idn_nme ; repmat(bid_nme(iD),numel(cov_sbj),1)];
            
        end
    end

    % Combine
    img_col = {'sbj_nme', 'ses_nme'};
    for iR = 1:numel(raw_fle_nme{iDT})
        img_col{end+1} = ['raw_file_' raw_fle_nme{iDT}{iR}];
        for iP = 1:numel(dta_typ_nme{iDT}{iR})
            img_col{end+1} = ['processed_file_' dta_typ_nme{iDT}{iR}{iP}];
        end
    end
    cmb_col_set = [ img_col {'dta_set_str'} {'dta_set_num'} cov_col ];
    cmb_dta_set = [ cmb_col_set ; tot_img_dta grp_idn_nme num2cell(grp_idn)  tot_cov_dta ];

    %% Save
    cell2csv([ new_dta_dir '/' '02_covariates' '/' typ_nme{iDT} '/' '00_' typ_nme{iDT} '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} cov_col {'dataset'};  tot_cov_sbj tot_cov_dta grp_idn_nme])

    cell2csv([ new_dta_dir '/' '03_Datasets' '/' '00_' typ_nme{iDT} '_' 'dataset' '_' dte_str '.csv'], cmb_dta_set);


    %% Make groups
    for iD = 1:numel(bid_nme)
        if any(grp_idn==iD)
            grp.tot_tbl.(bid_nme{iD}) = find(grp_idn==iD);
        end
    end
    grp.tot_tbl.total = [1:numel(grp_idn)]';
    fld_nme = fieldnames(grp.tot_tbl);

    %% Initialize
    tbl_typ = { 'count'    'Diagnosis'       'HC/EPD/other'              1 ; ...        
        'mean/std' 'Age'                     ''                          1 ; ...
        'mean/std' 'Education'               ''                          1 ; ...
        'count'    'Sex'                     'Male/Female/other'         1 ; ...
        'count'    'Handedness'              'Right/Left/Ambi'           1 ; ...
        'mean/std' 'AgeOfOnset'              ''                          1 ; ...
        'mean/std' 'Duration'                ''                          1 ; ...
        'mean/std' 'ASM_number'              ''                          1 ; ...        
        'count'    'Focal'                   'focal/generalized/mixed'   1; ...
        'count'    'Localization'            'temporal/temporal_plus/ex-tle-frontal/unknown'   1 ; ...
        'count'    'Lateralization'          'left/right/bilateral/multifocal/unknown'         1 ; ...
        'count'    'MTS'                     'noHS/HS_mri/HS_pathology/unknown'                1 ; ...
        'count'    'DrugResistance'          'refractory/responsive' 1 ; ...
        'count'    'EngelOutcome'            'I/II/III/IV'           1 ; ...
        'count'    'TonicClonic'             'Focal/GTC'             1 ; ...
        };


    %% Make Table
    % Create Table
    fcfg = [];

    for iR = 1:size(tbl_typ)
        if strcmpi(tbl_typ{iR,1},'mean/std')
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
            end
        elseif  strcmpi(tbl_typ{iR,1},'count')
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
            end
        end
    end

    fcfg.dta = {[ cov_col ; out_tot_dta ]};
    fcfg.grp = grp.tot_tbl;
    tbl_out  = ejk_create_table( fcfg );

    % Add-ons
    num_sbj{1} = 'N';
    tbl_lbl{1} = '';
    for iN = 1:numel(fld_nme)
        num_sbj{iN+1} = numel(grp.tot_tbl.(fld_nme{iN}));
        tbl_lbl{iN+1} = fld_nme{iN};
    end

    %% Save it
    tbl_row = cell(size(tbl_typ,1),1);
    for iR = 1:size(tbl_typ,1)
        if strcmpi(tbl_typ{iR,1},'count')
            tbl_row{iR,1} = [ tbl_typ{iR,2} ' (' tbl_typ{iR,3} ')'];
        elseif strcmpi(tbl_typ{iR,1},'mean/std')
            tbl_row{iR,1} = [ tbl_typ{iR,2} ' (' tbl_typ{iR,1} ')'];
        end
    end

    % Out
    tbl_out = [ tbl_lbl; num_sbj; tbl_row tbl_out ];

    cell2csv([ new_dta_dir '/' '00_reports'  '/' typ_nme{iDT} '/' typ_nme{iDT}  '_' 'covariates_table' '_' dte_str '.csv'],tbl_out)

end


