

t1w_num_out = mmil_readtext([ new_dta_dir '/' 'T1w' '_' 'count' '_' dte_str '.csv']);
t1w_num_out = [ t1w_num_out ; num2cell(zeros(1,size(t1w_num_out,2)))];

c12_fld_nme = fieldnames(drv_fld);

ejk_chk_dir(ult_t1w_c12_dir);

for iD = 1:numel(bid_nme)

    fprintf("CAT12: %s\n\n",bid_nme{iD})
    ejk_chk_dir([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' ]);
    if any(strcmpi(c12_fld_nme,bid_nme{iD}))
        bid_dir_use = [ bid_dir '/' bid_nme{iD} '/' ];
        drv_dir_use = [ bid_dir_use '/' 'derivatives' '/'];

        drv_nme = drv_fld.(bid_nme{iD});

        ses_dta = mmil_readtext([ new_dta_dir '/' bid_nme{iD} '/' bid_nme{iD} '_' 'sessions' '_' dte_str '.csv']);
        ses_dta_col_num = size(ses_dta,2);
        ses_dta_col_fle = ses_dta_col_num + numel(drv_nme) + 1;
        ses_dta_col_err = ses_dta_col_fle + 1;
        ses_dta = [ ses_dta cell(size(ses_dta,1),numel(drv_nme)+2)]; % Adding: state of derivatives * number of derivatives | final file | error codes

        sbj_nme_hld = cell(size(ses_dta,1),2);

        fprintf("CAT12: %s | searching\n",bid_nme{iD})
        for iS = 1:size(ses_dta,1)

            if rem(iS,100)==0
                fprintf("CAT12: %s | searching | subject #%i\n",bid_nme{iD}, iS)
            end

            ses_use_sbj = ses_dta{iS,4};

            for iDV = 1:numel(drv_nme)

                drv_iDV_dir = [ drv_dir_use '/' drv_nme{iDV} '/' ses_dta{iS,1} '/' ses_dta{iS,4} '/' 'anat'];
                drv_iDV_fle = dir(drv_iDV_dir); drv_iDV_fle = {drv_iDV_fle(:).name};
                if isempty(drv_iDV_fle)
                    ses_dta{iS,ses_dta_col_num+iDV} = 'Subject Not Present';
                else
                    drv_iDV_fle = drv_iDV_fle(string_find(drv_iDV_fle,'.nii'));


                    if isempty(string_find(drv_iDV_fle,'mwp1'))
                        ses_dta{iS,ses_dta_col_num+iDV} = 'Missing mwp1 file';
                    else
                        ses_dta{iS,ses_dta_col_num+iDV} = drv_iDV_fle{string_find(drv_iDV_fle,'mwp1')};
                    end

                    if isempty(ses_dta{iS,ses_dta_col_fle}) && ~strcmpi(ses_dta{iS,ses_dta_col_num+iDV},'Missing mwp1 file')
                        ses_dta{iS,ses_dta_col_fle} = [ drv_dir_use '/' drv_nme{iDV} '/' ses_dta{iS,1} '/' ses_dta{iS,4} '/' 'anat' '/' ses_dta{iS,ses_dta_col_num+iDV}];
                    end
                end
            end

            % Copy files over
            if ~isempty(ses_dta{iS,ses_dta_col_fle})
                sbj_nme_hld{iS,1} = ses_dta{iS,1};
                sbj_nme_hld{iS,2} = [ bid_nme{iD} '_' ses_dta{iS,1} '_' ses_dta{iS,4} '_' 'cat12_mwp1.nii'];
                if ~exist([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12' '/' 't1'  '/' sbj_nme_hld{iS,2} ])
                    copyfile( ses_dta{iS,ses_dta_col_fle}, [ ult_t1w_c12_dir '/' sbj_nme_hld{iS,2} ])
                end
            end
        end

        t1w_num_out{4,iD+1} = sum(~cellfun(@isempty,ses_dta(:,ses_dta_col_fle)));
        cell2csv([ new_dta_dir '/' bid_nme{iD} '/' bid_nme{iD} '_' 'sessions_CAT12' '_' dte_str '.csv'], ses_dta);
        
        sbj_nme_hld(cellfun(@isempty,sbj_nme_hld(:,1)),:) = [];       
        cell2csv([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv'], sbj_nme_hld);
    end
end

t1w_num_out{4,1} = 'CAT12';
t1w_num_out{4,end} = sum(cell2mat(t1w_num_out(4,2:end-1)));
cell2csv([ new_dta_dir '/' 'T1w' '_' 'CAT12' '_' 'count' '_' dte_str '.csv'],t1w_num_out);

