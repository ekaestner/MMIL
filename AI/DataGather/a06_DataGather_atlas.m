
for iDT = 1:numel(dta_typ_nme)

    new_atl_dir = [ dta_dir '/' dte_str '/' '04_atlases' '/' dta_typ_nme{iDT} '/'];
    ejk_chk_dir(new_atl_dir)

    %% Citations
    ejk_chk_dir([new_atl_dir '/' 'citations' '/' ])
    lbl_nme = dir(org_atl_dir); lbl_nme = {lbl_nme(:).name}; lbl_nme = lbl_nme(string_find(lbl_nme,'.txt'));
    for iL = 1:numel(lbl_nme)
        copyfile( [org_atl_dir '/' lbl_nme{iL}], [new_atl_dir '/' 'citations' '/' lbl_nme{iL}] );
    end

    %% labels
    ejk_chk_dir([new_atl_dir '/' 'labels' '/' ])
    lbl_nme = dir(org_atl_dir); lbl_nme = {lbl_nme(:).name}; lbl_nme = lbl_nme(string_find(lbl_nme,'.csv'));
    for iL = 1:numel(lbl_nme)
        copyfile( [org_atl_dir '/' lbl_nme{iL}], [new_atl_dir '/' 'labels' '/' lbl_nme{iL}] );
    end

    %% Atlases
    ejk_chk_dir([new_atl_dir '/' 'atlases' '/' ])
    ejk_chk_dir([new_atl_dir '/' 'atlases' '/' 'orig' '/'])

    atl_nme = dir(org_atl_dir); atl_nme = {atl_nme(:).name}; atl_nme = atl_nme(string_find(atl_nme,'.nii'));

    for iA = 1:numel(atl_nme); copyfile( [org_atl_dir '/' atl_nme{iA}], [new_atl_dir '/' 'atlases' '/' 'orig' '/' atl_nme{iA}] ); end

    for iRS = 1:numel(rsm_dim)
        ejk_chk_dir([new_atl_dir '/' 'atlases' '/'  num2str(rsm_dim(iRS)) '/'])
        atl_rpt = cell(numel(atl_nme),3);

        for iA = 1:numel(atl_nme)

            % Copy original file over
            nii_org = niftiread([new_atl_dir '/' 'atlases' '/' 'orig' '/' atl_nme{iA}]);

            %
            scl_num     = size(nii_org,1) / rsm_dim(iRS);
            scl_num_str = num2str(scl_num);

            % Rescale
            scl_cmd = 'tcsh -c';
            scl_cmd = sprintf('%s "%s',scl_cmd,src_cmd);
            scl_cmd = sprintf('%s; %s"',scl_cmd,[ thr_dim_rsm_cmd{1} scl_num_str ' ' scl_num_str ' ' scl_num_str ' ' thr_dim_rsm_cmd{2} [new_atl_dir '/' 'atlases' '/' 'orig' '/' atl_nme{iA}] ' ' thr_dim_rsm_cmd{3} [new_atl_dir '/' 'atlases' '/' num2str(rsm_dim(iRS)) '/' atl_nme{iA}]]);
            system(scl_cmd);

            %
            nii_scl = niftiread([new_atl_dir '/' 'atlases' '/' num2str(rsm_dim(iRS)) '/' atl_nme{iA}]);
            wrt_nii = 0;

            % Pad if necessary
            if ~isempty(pad_dim{iRS})
                nii_scl = padarray(nii_scl, [pad_dim{iRS}(1, 1), pad_dim{iRS}(1, 2), pad_dim{iRS}(1, 3)], 0, 'pre');
                nii_scl = padarray(nii_scl, [pad_dim{iRS}(2, 1), pad_dim{iRS}(2, 2), pad_dim{iRS}(2, 3)], 0, 'post');
                wrt_nii = 1;
            end

            % Trim if necessary
            if size(nii_scl,2)==138
                nii_scl = nii_scl(:,2:138,:);
                wrt_nii = 1;
            end

            % Write if necessary
            if wrt_nii
                niftiwrite(nii_scl,[new_atl_dir '/' 'atlases' '/' num2str(rsm_dim(iRS)) '/' atl_nme{iA}]);
            end

            %
            atl_rpt{iA,1} = atl_nme{iA};
            atl_rpt{iA,2} = size(nii_org);
            atl_rpt{iA,3} = size(nii_scl);

        end

        cell2csv([new_atl_dir '/' 'atlases' '/' num2str(rsm_dim(iRS)) '/' 'report' '_' dte_str '.csv'], [ {'Atlas'} {'original size'} {'Modified size'} ; atl_rpt])
    end
end
