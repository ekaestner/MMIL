t1w_num_out = mmil_readtext([ new_dta_dir '/' 'T1w' '_' 'CAT12' '_' 'count' '_' dte_str '.csv']);
t1w_num_out = [ t1w_num_out ; num2cell(zeros(3,size(t1w_num_out,2)))];

red_cap_fld_nme = fieldnames(drv_fld);
cde_bok_fld_nme = fieldnames(cde_bok);

key_fle = mmil_readtext([ prj_dir '/' 'new_data' '/' 'covariates' '/' 'key_fle.csv']);

for iD = 1:numel(bid_nme)

    if any(strcmpi(red_cap_fld_nme,bid_nme{iD})) && exist([ new_dta_dir '/' 'T1w' '_' 'CAT12' '_' 'count' '_' dte_str '.csv'],'file')
    
        ejk_chk_dir([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' 'reports' '/' ])

        key_row = strcmpi(key_fle(:,1),bid_nme{iD});

        % Load Redcap
        red_cap_dir = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'covariates' '/' ];

        fcfg = [];
        fcfg.dta_loc = [red_cap_dir '/' red_cap_fle.(bid_nme{iD})];
        [ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

        if strcmpi(bid_nme{iD},'enigma_conglom'); cnn_fig_format_ENIGMA_redcap; end

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
        ses_dta = mmil_readtext([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv']);

        cov_dta_out = cell(size(ses_dta,1),size(key_fle(:,2:end),2)+1);

        % Match Imaging Data to Covariate Data
        err_cnt = 0;
        for iS = 1:size(ses_dta,1)

            img_nme_col = strcmpi(red_cap_col,key_fle{key_row,strcmpi(key_fle(1,:),'Imaging_Name')});
            sbj_cov_row = string_find( red_cap_dta(:,img_nme_col), ses_dta{iS,1}(5:end) );

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
        cell2csv([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ]);
        
        % Put together missing reports
        cov_dta_hld = [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ];
        
        has_red_cap = cov_dta_hld(~cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),strcmpi(cov_dta_hld(1,:),'red_cap_nme'));
        mss_red_cap = cov_dta_hld(cellfun(@isempty,cov_dta_hld(:,strcmpi(cov_dta_hld(1,:),'Imaging_Name'))),1);

        cell2csv([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/'  'reports' '/' bid_nme{iD} '_' 'missing_Redcap' '_' dte_str '.csv'],mss_red_cap);

        [ ~, mss_img, ~] = setxor(red_cap_sbj,has_red_cap);
        cell2csv([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/'  'reports' '/' bid_nme{iD} '_' 'missing_Imaging' '_' dte_str '.csv'],[ {'sbj_nme'} red_cap_col ; red_cap_sbj(mss_img) red_cap_dta(mss_img,:)])
        
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
cell2csv([ new_dta_dir '/' 'T1w' '_' 'Redcap' '_' 'count' '_' dte_str '.csv'],t1w_num_out)
