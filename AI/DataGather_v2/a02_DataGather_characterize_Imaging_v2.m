

for iDT = 1:numel(typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' ]);
    img_fle_cll = cell(0);
    img_fle_num = 1;

    tot_out = cell(0);

    for iD = 1:numel(bid_nme)

        fprintf("%s : %s\n\n", typ_nme{iDT}, bid_nme{iD})

        % Load all raw scan lists
        raw_scn_all = cell(numel(raw_fle_nme{iDT}), 1);
        for iR = 1:numel(raw_fle_nme{iDT})
            try
                raw_scn_all{iR} = mmil_readtext([ new_dta_dir '/' '00_reports' '/' 'file_check' '/' typ_nme{iDT} '/' bid_nme{iD} '_' raw_fle_nme{iDT}{iR} '.txt']);
            catch
                raw_scn_all{iR} = cell(0);
            end
        end

        % Load all processed scan lists
        prc_scn_all = cell(numel(raw_fle_nme{iDT}), 1);
        for iR = 1:numel(raw_fle_nme{iDT})
            prc_scn_all{iR} = cell(numel(dta_typ_nme{iDT}{iR}), 1);
            for iP = 1:numel(dta_typ_nme{iDT}{iR})
                try
                    prc_scn_all{iR}{iP} = mmil_readtext([ new_dta_dir '/' '00_reports' '/' 'file_check' '/' typ_nme{iDT} '/' bid_nme{iD} '_' dta_typ_nme{iDT}{iR}{iP} '.txt']);
                catch
                    prc_scn_all{iR}{iP} = cell(0);
                end
            end
        end

        fprintf("%s : %s | searching\n", typ_nme{iDT}, bid_nme{iD})

        if isempty(raw_scn_all{1})
            fprintf("%s : %s | no raw files found, skipping\n", typ_nme{iDT}, bid_nme{iD})
            continue;
        end

        % Extract subjects/sessions from first raw file
        raw_pattern = 'sub-([^/_]+).*?ses-([^/_]+)';
        tokens = regexp(raw_scn_all{1}, raw_pattern, 'tokens');
        results = vertcat(tokens{:});
        sbj_dta = cellfun(@(x) x{1}, results, 'UniformOutput', false);
        ses_dta = cellfun(@(x) x{2}, results, 'UniformOutput', false);

        % Columns: sub, ses, [raw1, prc1_1, ..., raw2, prc2_1, prc2_2, ...]
        n_cols = 2;
        for iR = 1:numel(raw_fle_nme{iDT})
            n_cols = n_cols + 1 + numel(dta_typ_nme{iDT}{iR});
        end

        scn_dta_out = cell(size(sbj_dta,1), n_cols);
        for iS = 1:size(sbj_dta,1)

            if rem(iS,100)==0
                fprintf("%s : %s | searching | subject #%i\n", typ_nme{iDT}, bid_nme{iD}, iS)
            end

            scn_dta_out{iS,1} = [ 'sub-' sbj_dta{iS,1}];
            scn_dta_out{iS,2} = [ 'ses-' ses_dta{iS,1}];

            col = 3;
            for iR = 1:numel(raw_fle_nme{iDT})

                % Raw file column
                if iR == 1
                    scn_dta_out{iS,col} = raw_scn_all{1}{iS,1};
                elseif ~isempty(raw_scn_all{iR})
                    idx = contains(raw_scn_all{iR}, sbj_dta{iS,1}) & ...
                          contains(raw_scn_all{iR}, ses_dta{iS,1});
                    if any(idx)
                        scn_dta_out{iS,col} = raw_scn_all{iR}{find(idx,1)};
                    else
                        scn_dta_out{iS,col} = 'N/A';
                    end
                else
                    scn_dta_out{iS,col} = 'N/A';
                end
                col = col + 1;

                % Processed file columns for this raw
                for iP = 1:numel(dta_typ_nme{iDT}{iR})
                    if ~isempty(prc_scn_all{iR}{iP})
                        idx = contains(prc_scn_all{iR}{iP}, sbj_dta{iS,1}) & ...
                              contains(prc_scn_all{iR}{iP}, ses_dta{iS,1}) ;
                        if any(idx)
                            scn_dta_out{iS,col} = prc_scn_all{iR}{iP}{find(idx,1)};
                        else
                            scn_dta_out{iS,col} = 'N/A';
                        end
                    else
                        scn_dta_out{iS,col} = 'N/A';
                    end
                    col = col + 1;
                end

            end

        end

        cell2csv([ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' bid_nme{iD} '_' typ_nme{iDT} '_' 'file_characterization_' dte_str '.csv'], scn_dta_out);

        tot_out = [tot_out; scn_dta_out];

    end

    cell2csv([ new_dta_dir '/' '01_imaging_files' '/' typ_nme{iDT} '/' '00' '_' typ_nme{iDT} '_' 'file_characterization_' dte_str '.csv'], tot_out);

end