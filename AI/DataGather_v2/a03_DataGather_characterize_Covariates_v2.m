
for iDT = 1:numel(dta_typ_nme)

    red_cap_fld_nme = fieldnames(red_cap_fle);

    cde_bok_fld_nme = fieldnames(cde_bok);

    key_fle = mmil_readtext([ dta_dir '/' 'original_covariates' '/' 'key_fle_v2.csv']);

    ejk_chk_dir([ new_dta_dir '/' '02_covariates' '/' typ_nme{iDT} '/'])

    for iD = 1:numel(bid_nme)

        if any(strcmpi(red_cap_fld_nme,bid_nme{iD})) && ...
            exist([ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' bid_nme{iD} '_' typ_nme{iDT} '_' 'file_characterization_' dte_str '.csv'],'file') && ...
            ~isempty(red_cap_fle.(bid_nme{iD}))

            key_row = strcmpi(key_fle(:,1),bid_nme{iD});

            % Load Redcap
            red_cap_dir = [ dta_dir '/' 'original_covariates' '/' bid_nme{iD} '/' ];

            fcfg = [];
            fcfg.dta_loc = [red_cap_dir '/' red_cap_fle.(bid_nme{iD})];
            [ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

            if strcmpi(bid_nme{iD},'ENIGMA_dataset'); a03_DataGather_format_ENIGMA_redcap_v2; end
            if strcmpi(bid_nme{iD},'ECP_dataset_bids'); a03_DataGather_format_ECP_redcap_v2; end
            if strcmpi(bid_nme{iD},'IDEAS_dataset'); a03_DataGather_format_IDEAS_redcap_v2; end

            col_use = nan(1,size(key_fle,2)-1);
            for iK = 2:size(key_fle,2)
                if ~isempty(key_fle{key_row,iK}) && ~strcmpi(key_fle{key_row,iK},'empty_string') && ~strcmpi(key_fle{key_row,iK},'empty_numeric')
                    col_use(iK-1) = find(strcmpi(red_cap_col,key_fle{key_row,iK}));
                elseif strcmpi(key_fle{key_row,iK},'empty_string')
                    red_cap_dta(:,end+1) = repmat({''},size(red_cap_dta,1),1);
                    red_cap_col(:,end+1) = key_fle(key_row,iK);
                    col_use(iK-1) = size(red_cap_dta,2);
                elseif strcmpi(key_fle{key_row,iK},'empty_numeric')
                    red_cap_dta(:,end+1) = repmat({NaN},size(red_cap_dta,1),1);
                    red_cap_col(:,end+1) = key_fle(key_row,iK);
                    col_use(iK-1) = size(red_cap_dta,2);
                end
            end

            % Load Imaging Data
            img_fle_pth = [ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' bid_nme{iD} '_' typ_nme{iDT} '_' 'file_characterization_' dte_str '.csv'];
            
            opts = detectImportOptions(img_fle_pth, 'Delimiter', ',', 'FileType', 'text');
            opts.VariableNamesLine = 0;      % Treat headers as data to ensure everything is read
            opts.DataLines = [1 Inf];
            opts = setvartype(opts, 'char'); % Strictly force all columns to be 'char'
            
            T = readtable(img_fle_pth, opts);
            ses_dta = table2cell(T);         % Convert back to cell array structure

            cov_dta_out = cell(size(ses_dta,1),size(key_fle(:,2:end),2)+1);

            % Match Imaging Data to Covariate Data
            err_cnt = 0;
            for iS = 1:size(ses_dta,1)

                if mod(iS, max(1, floor(size(ses_dta,1)/5))) == 0
                    fprintf('Processing session %d/%d\n', iS, size(ses_dta,1));
                end

                sbj_cov_row = string_find( red_cap_sbj, ses_dta{iS,1} );
                
                if ~isempty(sbj_cov_row)
                    if numel(sbj_cov_row)==1
                        cov_dta_out{iS,1} = red_cap_sbj{sbj_cov_row,1};
                        cov_dta_out(iS,2:end) = red_cap_dta(sbj_cov_row,col_use);
                    else
                        err_cnt = err_cnt+1;
                    end
                end
            end

            % Recode
            if any(strcmpi(cde_bok_fld_nme,bid_nme{iD}))
                cov_nme = fieldnames(cde_bok.(bid_nme{iD}));
                cov_val = cell(1,numel(cov_nme));
                for iCV = 1:numel(cov_nme)
                    cov_val{iCV} = cde_bok.(bid_nme{iD}).(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
                end

                fcfg = [];

                fcfg.dta     = cov_dta_out;
                fcfg.dta_col = key_fle(key_row,:);

                fcfg.rcd_nme = cov_nme;
                fcfg.rcd_val = cov_val;

                fcfg.swt_val = 1;

                cov_dta_out = ejk_recode(fcfg);
            end

            if strcmpi(bid_nme{iD},'ENIGMA_dataset'); cnn_fig_format_ENIGMA_redcap_post; end

            % Save
            cell2csv([ new_dta_dir '/' '02_covariates' '/' typ_nme{iDT} '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ]);

            % Put together missing reports
            cov_dta_hld = [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ];

            has_red_cap = cov_dta_hld(~cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),strcmpi(cov_dta_hld(1,:),'red_cap_nme'));
            mss_red_cap = cov_dta_hld(cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),1);

            ejk_chk_dir([ new_dta_dir '/' '00_reports' '/' typ_nme{iDT} '/' 'missing_redcap' '/'])
            cell2csv([ new_dta_dir '/' '00_reports' '/' typ_nme{iDT} '/' 'missing_redcap' '/' bid_nme{iD} '_' 'missing_redcap' '_' dte_str '.csv'],mss_red_cap);

            [ ~, mss_img, ~] = setxor(red_cap_sbj,has_red_cap);
            ejk_chk_dir([ new_dta_dir '/' '00_reports' '/' typ_nme{iDT} '/' 'missing_imaging' '/'])
            cell2csv([ new_dta_dir '/' '00_reports' '/' typ_nme{iDT} '/' 'missing_imaging' '/' bid_nme{iD} '_' 'missing_imaging' '_' dte_str '.csv'],[ {'sbj_nme'} red_cap_col ; red_cap_sbj(mss_img) red_cap_dta(mss_img,:)])

        end
    end

end