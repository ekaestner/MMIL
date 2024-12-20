onb_dir = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding';

% dte_str = '2024_10_08';

%% Load
ins_dta = mmil_readtext([ onb_dir '/' 'institutions.csv' ]);
% ins_dta = mmil_spec_char(ins_dta,{';'},{','});

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{1} '/' 'covariates' '/' 'ENIGMAEpilepsy-Subjidcheck_DATA_LABELS_2024-10-08_1139.csv'];
[ red_sbj_dta, red_sbj_sbj, red_sbj_col ] = ejk_dta_frm(fcfg);
cnn_fig_format_ENIGMA_redcap_subjid_check

%% Site Specific Reports
for iI = 1:size(ins_dta,1)
    if ~isempty(ins_dta{iI,1}) || ~isempty(ins_dta{iI,2})

        red_sbj_dta_ins = red_sbj_dta(strcmpi(red_sbj_dta(:,strcmpi(red_sbj_col,'Institution Name')),ins_dta{iI,2}),:);
        red_sbj_sbj_ins = red_sbj_sbj(strcmpi(red_sbj_dta(:,strcmpi(red_sbj_col,'Institution Name')),ins_dta{iI,2}),:);

        % %%%%
        img_fle = dir([ onb_dir '/' ins_dta{iI} '/' ]);
        img_fle = { img_fle(:).name};
        img_fle = img_fle(string_find(img_fle,'sub-'));

        % %%%%
        ins_out = cell(numel(img_fle),6);

        red_sbj_rmv = [];
        for iS = 1:numel(img_fle)

            % %%%%%%
            mtc_ind = strcmpi(red_sbj_dta_ins(:,strcmpi(red_sbj_col,'BIDS Subject ID')),img_fle{iS});

            if any(mtc_ind) && (sum(mtc_ind)==1)

                mtc_ind_fnd = find(mtc_ind);

                red_sbj_rmv = [ red_sbj_rmv ; mtc_ind_fnd];

                ins_out{iS,1} = red_sbj_sbj_ins{mtc_ind_fnd};
                ins_out{iS,2} = img_fle{iS};
                ins_out(iS,3:5) = red_sbj_dta_ins(mtc_ind_fnd,3:5);
                ins_out{iS,6} = 'Match';
            else
                ins_out{iS,2} = img_fle{iS};
                ins_out{iS,6} = 'BIDS only';
            end
        end

        red_sbj_dta_ins(red_sbj_rmv,:) = [];
        red_sbj_sbj_ins(red_sbj_rmv,:) = [];

        ins_sve = [ {'redcap_sbj_nme'} {'bids_sbj_nme'} red_sbj_col(3:5) {'Class'} ; ins_out ; red_sbj_sbj_ins cell(numel(red_sbj_sbj_ins),1) red_sbj_dta_ins(:,3:5) repmat({'Redcap Only'},numel(red_sbj_sbj_ins),1)];

        if ~isempty(ins_dta{iI,1})
            ejk_chk_dir([ onb_dir '/' ins_dta{iI} '/' 'covariates' '/' 'original' '/'])
            cell2csv([ onb_dir '/' ins_dta{iI} '/' 'covariates' '/' 'original' '/' ins_dta{iI} '_' 'matching' '_' dte_str '.csv'],ins_sve)
        elseif isempty(ins_dta{iI,1})
            cell2csv([ onb_dir '/' 'reports' '/' 'redcap_hold' '/' num2str(iI) '_' 'matching' '_' dte_str '.csv'],ins_sve)
        end

    end
end

%% Total Report
rpt_out = cell(size(ins_dta,1),5);
for iI = 1:size(ins_dta,1)
    
    rpt_out(iI,1:2) = ins_dta(iI,1:2);
    
    if ~isempty(ins_dta{iI,1}) || ~isempty(ins_dta{iI,2})

        if ~isempty(ins_dta{iI,1})
            fcfg = [];
            fcfg.dta_loc = [ onb_dir '/' ins_dta{iI} '/' 'covariates' '/' 'original' '/' ins_dta{iI} '_' 'matching' '_' dte_str '.csv'];
            [ rpt_dta, rpt_sbj, rpt_col ] = ejk_dta_frm(fcfg);
        elseif isempty(ins_dta{iI,1})
            fcfg = [];
            fcfg.dta_loc = [ onb_dir '/' 'reports' '/' 'redcap_hold' '/' num2str(iI) '_' 'matching' '_' dte_str '.csv'];
            [ rpt_dta, rpt_sbj, rpt_col ] = ejk_dta_frm(fcfg);
        end

        rpt_out{iI,3} = sum(strcmpi(rpt_dta(:,strcmpi(rpt_col,'class')),'Match'));
        rpt_out{iI,4} = sum(strcmpi(rpt_dta(:,strcmpi(rpt_col,'class')),'BIDS only'));
        rpt_out{iI,5} = sum(strcmpi(rpt_dta(:,strcmpi(rpt_col,'class')),'Redcap Only'));

    end
end

rpt_sve = [ {'Site'} {'Folder'} {'Match'} {'BIDS only'} {'Redcap Only'} ; rpt_out ; {'Total'} {''} {sum(cell2mat(rpt_out(:,3)))} {sum(cell2mat(rpt_out(:,4)))} {sum(cell2mat(rpt_out(:,5)))}];

cell2csv([ onb_dir '/' 'reports' '/' 'matching_report_' dte_str '.csv'],rpt_sve);


