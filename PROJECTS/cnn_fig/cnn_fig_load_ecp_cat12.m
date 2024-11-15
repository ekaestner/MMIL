
%% Identify CAT12 Runs (Donatello run)
c12_sbj_lst = dir(ecp_c12_drv_dir); c12_sbj_lst = { c12_sbj_lst(:).name }; c12_sbj_lst = c12_sbj_lst(string_find(c12_sbj_lst,'sub-'));
err_lst = {};

c12_dta_out = [ c12_sbj_lst cell(numel(c12_sbj_lst),2) ];

for iS = 1:numel(c12_sbj_lst)

    % Find possible session dirs
    ses_dir = dir([ ecp_c12_drv_dir '/' c12_sbj_lst{iS} ]); ses_dir = {ses_dir(:).name};
    ses_dir = ses_dir(string_find(ses_dir,'ses'));

    % Identify possible files
    mwp_pth = [ ecp_c12_drv_dir '/' c12_sbj_lst{iS} '/' ses_dir{1} '/' 'anat' '/'];
    mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
    try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end

    if ~isempty(mwp_fle)
        c12_dta_out{iS,2} = mwp_fle;
        c12_dta_out{iS,3} = [ mwp_pth '/' mwp_fle];
    else
        c12_dta_out{iS,2} = 'MISSING_mwp_fle';
    end

end

%% Check the outputs %%%
clear sum_tbl
sum_tbl{1,2} = numel(string_find(c12_dta_out(:,2),'.nii'));                                           sum_tbl{1,1} = 'derivatives/cat12 has mwp file';
sum_tbl{2,2} = sum(strcmpi(c12_dta_out(:,2),'MISSING_mwp_fle') & cellfun(@isempty,c12_dta_out(:,3))); sum_tbl{2,1} = 'derivatives/cat12 missing mwp file';

sum_tbl

%% Save Out %%%
cell2csv([ dta_set_ifo '/' 'ECP_CAT12_wmp1_file_locations' '_' dte_str '.csv'],[{'subject'} {'derivatives/cat12'} {'file_location'} ; c12_dta_out])
cell2csv([ dta_set_ifo '/' 'ECP_CAT12_wmp1_file_totals' '_' dte_str '.csv'],sum_tbl);

%% Copy Files over %%%
for iS = 1:size(c12_dta_out,1)
    if ~isempty(c12_dta_out{iS,3}) && ~exist([ t1w_c12_ecp_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii'],'file') 
        copyfile(c12_dta_out{iS,3}, [ t1w_c12_ecp_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii']);
    end
end