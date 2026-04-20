
for iDT = 1:numel(dta_typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '07_imaging_files' '/' dta_typ_nme{iDT}]);
    img_fle_cll = cell(0);
    img_fle_num = 1;

    t1w_num_out = mmil_readtext([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' 'count' '_' dte_str '.csv']);
    t1w_num_out = [ t1w_num_out ; num2cell(zeros(1,size(t1w_num_out,2)))];

    c12_fld_nme = fieldnames(drv_fld{iDT});

    for iD = 1:numel(bid_nme)

        fprintf("%s: %s\n\n",dta_typ_nme{iDT}, bid_nme{iD})
        ejk_chk_dir([ new_dta_dir '/' '02_subjects' '/' dta_typ_nme{iDT} '/']);
        if any(strcmpi(c12_fld_nme,bid_nme{iD}))
            bid_dir_use = [ bid_dir '/' bid_nme{iD} '/' ];
            drv_dir_use = [ bid_dir_use '/' 'derivatives' '/'];

            drv_nme = drv_fld{iDT}.(bid_nme{iD});

            ses_dta = mmil_readtext([ new_dta_dir '/' '01_sessions' '/' dta_typ_nme{iDT} '/' bid_nme{iD} '_' 'sessions' '_' dte_str '.csv']);
            ses_dta_col_num = size(ses_dta,2);
            ses_dta_col_fle = ses_dta_col_num + numel(drv_nme) + 1;
            ses_dta_col_err = ses_dta_col_fle + 1;
            ses_dta = [ ses_dta cell(size(ses_dta,1),numel(drv_nme)+2)]; % Adding: state of derivatives * number of derivatives | final file | error codes

            sbj_nme_hld = cell(size(ses_dta,1),3);

            fprintf("%s: %s | searching\n",dta_typ_nme{iDT},bid_nme{iD})
            for iS = 1:size(ses_dta,1)

                if rem(iS,100)==0
                    fprintf("%s: %s | searching | subject #%i\n",dta_typ_nme{iDT},bid_nme{iD}, iS)
                end

                ses_use_sbj = ses_dta{iS,4};

                for iDV = 1:numel(drv_nme)

                    drv_iDV_dir = [ drv_dir_use '/' drv_nme{iDV} '/' ses_dta{iS,1} '/' ses_dta{iS,4} '/' fld_ptt{iDT} '/'];                    
                    drv_iDV_fle = dir(drv_iDV_dir); drv_iDV_fle = {drv_iDV_fle(:).name};

                    if isempty(drv_iDV_fle)
                        drv_iDV_dir = [ drv_dir_use '/' drv_nme{iDV} '/' strrep(ses_dta{iS,1},'sub-','') '/' ses_dta{iS,4} '/'  fld_ptt{iDT} '/'];
                        drv_iDV_fle = dir(drv_iDV_dir); drv_iDV_fle = {drv_iDV_fle(:).name};
                    end

                    if isempty(drv_iDV_fle)
                        ses_dta{iS,ses_dta_col_num+iDV} = 'Subject Not Present';
                    else
                        drv_iDV_fle = drv_iDV_fle(string_find(drv_iDV_fle,'.nii'));

                        if isempty(string_find(drv_iDV_fle,fle_ptt{iDT}))
                            ses_dta{iS,ses_dta_col_num+iDV} = ['Missing ' fle_ptt{iDT} ' file'];
                        else
                            ses_dta{iS,ses_dta_col_num+iDV} = drv_iDV_fle{string_find(drv_iDV_fle,fle_ptt{iDT})};
                        end

                        if isempty(ses_dta{iS,ses_dta_col_fle}) && ~strcmpi(ses_dta{iS,ses_dta_col_num+iDV},['Missing ' fle_ptt{iDT} ' file'])
                            ses_dta{iS,ses_dta_col_fle} = [ drv_dir_use '/' drv_nme{iDV} '/' ses_dta{iS,1} '/' ses_dta{iS,4} '/' fld_ptt{iDT} '/' ses_dta{iS,ses_dta_col_num+iDV}];
                            if ~exist(ses_dta{iS,ses_dta_col_fle})
                                ses_dta{iS,ses_dta_col_fle} = [ drv_dir_use '/' drv_nme{iDV} '/' strrep(ses_dta{iS,1},'sub-','') '/' ses_dta{iS,4} '/' fld_ptt{iDT} '/' ses_dta{iS,ses_dta_col_num+iDV}];
                                if ~exist(ses_dta{iS,ses_dta_col_fle})
                                    error('Check paths')
                                end
                            end
                        end
                    end
                end

                % Copy files over
                if ~isempty(ses_dta{iS,ses_dta_col_fle})
                    sbj_nme_hld{iS,1} = ses_dta{iS,1};
                    sbj_nme_hld{iS,2} = [ bid_nme{iD} '_' ses_dta{iS,1} '_' ses_dta{iS,4} '_' fle_ptt{iDT}];
                    sbj_nme_hld{iS,3} = ses_dta{iS,ses_dta_col_fle};
                    if ~exist([ new_dta_dir '/' '07_imaging_files' '/' dta_typ_nme{iDT} '/' sbj_nme_hld{iS,2}])
                        new_fle_loc = [ new_dta_dir '/' '07_imaging_files' '/' dta_typ_nme{iDT} '/' sbj_nme_hld{iS,1} '.nii'];
                        copyfile( ses_dta{iS,ses_dta_col_fle}, new_fle_loc)
                        img_fle_cll{img_fle_num} = sprintf('copyfile( %s, %s)',ses_dta{iS,ses_dta_col_fle},new_fle_loc);
                        img_fle_num = img_fle_num + 1;
                    end
                end
            end

            t1w_num_out{4,iD+1} = sum(~cellfun(@isempty,ses_dta(:,ses_dta_col_fle)));
            cell2csv([ new_dta_dir '/' '01_sessions' '/' dta_typ_nme{iDT} '/' bid_nme{iD} '_' 'sessions_' dta_typ_nme{iDT} '_' dte_str '.csv'], ses_dta);

            sbj_nme_hld(cellfun(@isempty,sbj_nme_hld(:,1)),:) = [];
            cell2csv([ new_dta_dir '/' '02_subjects' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv'], sbj_nme_hld);
        end
    end

    t1w_num_out{4,1} = dta_typ_nme{iDT};
    t1w_num_out{4,end} = sum(cell2mat(t1w_num_out(4,2:end-1)));
    cell2csv([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' dta_typ_nme{iDT} '_' 'count' '_' dte_str '.csv'],t1w_num_out);

    cell2csv([ new_dta_dir '/' '00_reports' '/' dta_typ_nme{iDT} '/' 'copy_file_log.csv'],img_fle_cll)

end
