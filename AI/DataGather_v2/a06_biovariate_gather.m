%% Setup
ejk_chk_dir([ new_dta_dir '/' '07_biovariates' '/'])

%% Find all PicklePrep files
pkl_dir = [ new_dta_dir '/' '05_PicklePreps' '/'];
pkl_fls = dir([ pkl_dir '*.csv']);

%% Loop over each PicklePrep file
for iPP = 1:numel(pkl_fls)

    pkl_nme = pkl_fls(iPP).name(1:end-4);

    ejk_chk_dir([ new_dta_dir '/' '07_biovariates' '/' pkl_nme '/'])

    %% Load PicklePrep
    pkl_fle_pth = [ pkl_dir pkl_fls(iPP).name ];

    opts = detectImportOptions(pkl_fle_pth, 'Delimiter', ',', 'FileType', 'text');
    opts.VariableNamesLine = 0;
    opts.DataLines         = [1 Inf];
    opts = setvartype(opts, 'char');
    T       = readtable(pkl_fle_pth, opts);
    pkl_raw = table2cell(T);
    pkl_hdr = pkl_raw(1,:);
    pkl_dta = pkl_raw(2:end,:);

    % Identify key columns by name
    sbj_col_ind = find(strcmpi(pkl_hdr, 'sbj_nme'));
    ses_col_ind = find(strcmpi(pkl_hdr, 'ses_nme'));
    dts_col_ind = find(strcmpi(pkl_hdr, 'dta_set_str'));

    pkl_dts_str = pkl_dta(:, dts_col_ind);
    unq_dts     = unique(pkl_dts_str);

    %% Loop over derivative folders
    for iF = 1:numel(drv_fld)

        %% Loop over measure files within this derivative folder
        for iM = 1:numel(drv_mes{iF})

            mes_nme  = drv_mes{iF}{iM};
            bio_col  = {};
            bio_dta  = repmat({'NaN'}, size(pkl_dta,1), 0);
            col_init = false;
            n_found  = 0;
            n_miss   = 0;

            out_fle = [ new_dta_dir '/' '07_biovariates' '/' pkl_nme '/' mes_nme ];
            fprintf('Starting: %s\n', out_fle)

            %% Loop over unique datasets present in this PicklePrep
            for iD = 1:numel(unq_dts)

                dts_nme = unq_dts{iD};

                % Skip datasets with no drv_dte entry
                if ~isfield(drv_dte, dts_nme)
                    fprintf('  [skip] no drv_dte for %s\n', dts_nme);
                    continue
                end
                if iF > numel(drv_dte.(dts_nme))
                    fprintf('  [skip] drv_dte index %d missing for %s\n', iF, dts_nme);
                    continue
                end

                % Construct path to biovariate file
                bio_fle_pth = [ bid_dir dts_nme '/derivatives/' drv_fld{iF} ...
                                '/Measures/' drv_dte.(dts_nme){iF} '/' mes_nme ];
                if ~exist(bio_fle_pth, 'file')
                    bio_fle_pth = [ bid_dir dts_nme '/derivatives/' drv_fld{iF} ...
                                '/Measure/' drv_dte.(dts_nme){iF} '/' mes_nme ];
                end

                if ~exist(bio_fle_pth, 'file')
                    fprintf('  [skip] file not found: %s\n', bio_fle_pth);
                    continue
                end

                % Load biovariate file
                opts2 = detectImportOptions(bio_fle_pth, 'Delimiter', ',', 'FileType', 'text');
                opts2.VariableNamesLine = 0;
                opts2.DataLines         = [1 Inf];
                opts2 = setvartype(opts2, 'char');
                T2      = readtable(bio_fle_pth, opts2);
                bio_raw = table2cell(T2);
                bio_hdr = bio_raw(1,:);
                bio_tbl = bio_raw(2:end,:);

                % Identify subject/session index columns in biovariate file
                bio_sbj_col = find(strcmpi(bio_hdr, 'subject'));
                bio_ses_col = find(strcmpi(bio_hdr, 'session'));
                bio_dat_col = setdiff(1:numel(bio_hdr), [bio_sbj_col bio_ses_col]);

                % On first successful load, initialize output columns;
                % on subsequent loads, restrict bio_dat_col to match the
                % reference column set (bio_col) established on first load.
                if ~col_init
                    bio_col  = bio_hdr(bio_dat_col);
                    bio_dta  = repmat({'NaN'}, size(pkl_dta,1), numel(bio_col));
                    col_init = true;
                else
                    % Map each reference column to its index in this file's header
                    bio_dat_col = cellfun(@(c) find(strcmpi(bio_hdr, c), 1), bio_col, ...
                                         'UniformOutput', true);
                end

                % Rows in the pickle that belong to this dataset
                pkl_rows = find(strcmpi(pkl_dts_str, dts_nme));

                % Match each pickle row to the corresponding biovariate row
                for iR = 1:numel(pkl_rows)

                    rIdx = pkl_rows(iR);

                    % Strip BIDS prefixes to get bare IDs used in biovariate files
                    sbj_id = strrep(pkl_dta{rIdx, sbj_col_ind}, 'sub-', '');
                    ses_id = strrep(pkl_dta{rIdx, ses_col_ind}, 'ses-', '');

                    % Find matching row(s)
                    if ~isempty(bio_ses_col)
                        bio_match = find(strcmpi(bio_tbl(:, bio_sbj_col), sbj_id) & ...
                                         strcmpi(bio_tbl(:, bio_ses_col), ses_id));
                    else
                        bio_match = find(strcmpi(bio_tbl(:, bio_sbj_col), sbj_id));
                    end

                    if numel(bio_match) == 1
                        bio_dta(rIdx, :) = bio_tbl(bio_match, bio_dat_col);
                        n_found = n_found + 1;
                    elseif numel(bio_match) > 1
                        fprintf('  [warn] multiple matches for %s %s in %s \x2014 using first\n', ...
                                sbj_id, ses_id, dts_nme);
                        bio_dta(rIdx, :) = bio_tbl(bio_match(1), bio_dat_col);
                        n_found = n_found + 1;
                    else
                        n_miss = n_miss + 1;
                    end
                    % numel == 0: row stays 'NaN'

                end % iR

            end % iD

            %% Save output
            if col_init
                out_hdr = [ {'sbj_nme'} bio_col ];
                out_dta = [ pkl_dta(:, sbj_col_ind) bio_dta ];
                out_csv = [ out_hdr ; out_dta ];

                fprintf('  rows found: %d | rows missing: %d\n', n_found, n_miss)
                cell2csv(out_fle, out_csv)
                fprintf('Saved: %s\n', out_fle)
            else
                fprintf('  [skip] no biovariate data found for %s / %s / %s\n', ...
                        pkl_nme, drv_fld{iF}, mes_nme);
            end

        end % iM
    end % iF

end % iPP