
for iDT = 1:numel(dta_typ_nme)

    t1w_num_out = mmil_readtext([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' 'count' '_' dte_str '.csv']);
    t1w_num_out = [ t1w_num_out ; num2cell(zeros(3,size(t1w_num_out,2)))];

    red_cap_fld_nme = fieldnames(drv_fld{iDT});
    cde_bok_fld_nme = fieldnames(cde_bok);

    key_fle = mmil_readtext([ dta_dir '/' 'original_covariates' '/' 'key_fle.csv']);

    ejk_chk_dir([ new_dta_dir '/' '03_covariates' '/' dta_typ_nme{iDT} '/'])

    for iD = 1:numel(bid_nme)

        if any(strcmpi(red_cap_fld_nme,bid_nme{iD})) && exist([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' 'count' '_' dte_str '.csv'],'file')

            key_row = strcmpi(key_fle(:,1),bid_nme{iD});

            % Load Redcap
            red_cap_dir = [ dta_dir '/' 'original_covariates' '/' bid_nme{iD} '/' ];

            fcfg = [];
            fcfg.dta_loc = [red_cap_dir '/' red_cap_fle.(bid_nme{iD})];
            [ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

            if strcmpi(bid_nme{iD},'enigma_conglom'); a02_DataGather_format_ENIGMA_redcap; end
            if strcmpi(bid_nme{iD},'ECP_dataset_bids'); a02_DataGather_format_ECP_redcap; end

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
            ses_dta = mmil_readtext([ new_dta_dir '/' '02_subjects' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv']);

            cov_dta_out = cell(size(ses_dta,1),size(key_fle(:,2:end),2)+1);

            % Match Imaging Data to Covariate Data
            err_cnt = 0;
            for iS = 1:size(ses_dta,1)

                img_nme_col = strcmpi(red_cap_col,key_fle{key_row,strcmpi(key_fle(1,:),'Imaging_Name')});
                sbj_cov_row = string_find( red_cap_dta(:,img_nme_col), ses_dta{iS,1}(5:end) );
                if strcmpi(bid_nme{iD},'enigma_conglom') && isempty(sbj_cov_row); sbj_cov_row = string_find( red_cap_sbj, ses_dta{iS,1}); end

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

            if strcmpi(bid_nme{iD},'enigma_conglom'); cnn_fig_format_ENIGMA_redcap_post; end

            % Save
            cell2csv([ new_dta_dir '/' '03_covariates' '/' dta_typ_nme{iDT} '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ]);

            % Put together missing reports
            cov_dta_hld = [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ];

            has_red_cap = cov_dta_hld(~cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),strcmpi(cov_dta_hld(1,:),'red_cap_nme'));
            mss_red_cap = cov_dta_hld(cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),1);

            ejk_chk_dir([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' 'missing_redcap' '/'])
            cell2csv([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' 'missing_redcap' '/' bid_nme{iD} '_' 'missing_redcap' '_' dte_str '.csv'],mss_red_cap);

            [ ~, mss_img, ~] = setxor(red_cap_sbj,has_red_cap);
            ejk_chk_dir([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' 'missing_imaging' '/'])
            cell2csv([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' 'missing_imaging' '/' bid_nme{iD} '_' 'missing_imaging' '_' dte_str '.csv'],[ {'sbj_nme'} red_cap_col ; red_cap_sbj(mss_img) red_cap_dta(mss_img,:)])

            % Tally
            t1w_num_out{5,iD+1} = sum(~cellfun(@isempty,cov_dta_out(:,2)));
            t1w_num_out{6,iD+1} = sum(~cellfun(@isempty,cov_dta_out(:,2)) & ~cellfun(@isempty,ses_dta(:,end)));
            t1w_num_out{7,iD+1} = sum(~cellfun(@isempty,cov_dta_out(:,2)) & ~cellfun(@isempty,cov_dta_out(:,end)));

        end
    end

    t1w_num_out{5,1} = 'Redcap Match';
    t1w_num_out{5,end} = sum(cell2mat(t1w_num_out(5,2:end-1)));
    t1w_num_out{6,1} = 'CAT12 Match';
    t1w_num_out{6,end} = sum(cell2mat(t1w_num_out(6,2:end-1)));
    t1w_num_out{7,1} = 'Diagnosis Known';
    t1w_num_out{6,end} = sum(cell2mat(t1w_num_out(6,2:end-1)));
    cell2csv([  new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' 'count' '_' dte_str '.csv'],t1w_num_out)

end