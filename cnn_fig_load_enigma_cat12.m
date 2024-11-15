
%% Identify CAT12 Runs (Donatello run)
c12_ver_one_sbj_lst = dir(c12_drv_ver_one_dir); c12_ver_one_sbj_lst = { c12_ver_one_sbj_lst(:).name }; c12_ver_one_sbj_lst = c12_ver_one_sbj_lst(string_find(c12_ver_one_sbj_lst,'sub-'));
c12_ver_two_sbj_lst = dir(c12_drv_ver_two_dir); c12_ver_two_sbj_lst = { c12_ver_two_sbj_lst(:).name }; c12_ver_two_sbj_lst = c12_ver_two_sbj_lst(string_find(c12_ver_two_sbj_lst,'sub-'));
err_lst = {};

tot_sbj_lst = unique( [ c12_ver_one_sbj_lst' ; c12_ver_two_sbj_lst' ]);

c12_dta_out = [ tot_sbj_lst cell(numel(tot_sbj_lst),3) ];

for iS = 1:numel(tot_sbj_lst)

    % c12_ver_one_sbj_lst %%%
    ses_dir = dir([ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} ]); ses_dir = {ses_dir(:).name};
    if any(string_find(ses_dir,'ses'))
        ses_dir = ses_dir{string_find(ses_dir,'ses')};
        mwp_pth = [ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} '/' ses_dir '/' 'anat' '/'];
        mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
        try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        
    else
        mwp_pth = [ c12_drv_ver_one_dir '/' tot_sbj_lst{iS} '/' 'anat' '/'];
        mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
        try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
    end
    if ~isempty(mwp_fle)
        c12_dta_out{iS,2} = mwp_fle;
        c12_dta_out{iS,4} = [ mwp_pth '/' mwp_fle];
    else
        c12_dta_out{iS,2} = 'MISSING_mwp_fle';
    end

    % c12_ver_two_sbj_lst %%%
    if any(strcmpi(c12_ver_two_sbj_lst,tot_sbj_lst{iS}))
        ses_dir = dir([ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} ]); ses_dir = {ses_dir(:).name};
        if any(string_find(ses_dir,'ses'))
            ses_dir = ses_dir{string_find(ses_dir,'ses')};
            mwp_pth = [ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} '/' ses_dir '/' 'anat' '/'];
            mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
            try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        else
            mwp_pth = [ c12_drv_ver_two_dir '/' tot_sbj_lst{iS} '/' 'anat' '/'];
            mwp_fle = dir(mwp_pth); mwp_fle = {mwp_fle(:).name}; mwp_fle = mwp_fle(string_find(mwp_fle,'.nii'));
            try mwp_fle = mwp_fle{string_find(mwp_fle,'mwp1')}; catch mwp_fle = cell(0); end
        end
        if ~isempty(mwp_fle)
            c12_dta_out{iS,3} = mwp_fle;
            if isempty(c12_dta_out{iS,4})
                c12_dta_out{iS,4} = [ mwp_pth '/' mwp_fle];
            else
                c12_dta_out{iS,5} = 'Double Files';
            end
        else
            c12_dta_out{iS,3} = 'MISSING_mwp_fle';
        end
    end
end

%% Check the outputs %%%
clear sum_tbl
sum_tbl{1,2} = numel(string_find(c12_dta_out(:,2),'.nii')); sum_tbl{1,1} = 'derivatives/cat12 has mwp file';
sum_tbl{2,2} = numel(string_find(c12_dta_out(:,3),'.nii')); sum_tbl{2,1} = 'derivatives/cat12_2 has mwp file';

sum_tbl{3,2} = sum(strcmpi(c12_dta_out(:,2),'MISSING_mwp_fle') & cellfun(@isempty,c12_dta_out(:,3))); sum_tbl{3,1} = 'derivatives/cat12 missing mwp file';
sum_tbl{4,2} = sum(strcmpi(c12_dta_out(:,3),'MISSING_mwp_fle')); sum_tbl{4,1} = 'derivatives/cat12_2 missing mwp file';

sum_tbl{5,2} = sum(strcmpi(c12_dta_out(:,5),'Double Files')); sum_tbl{5,1} = 'cat12 & cat12_2 both have file';

sum_tbl

%% Save Out %%%
cell2csv([ prj_dir '/' 'CAT12_wmp1_file_locations' '_' dte_str '.csv'],[{'subject'} {'derivatives/cat12'} {'derivatives/cat12_2'} {'file_location'} {'notes'}; c12_dta_out])
cell2csv([ prj_dir '/' 'CAT12_wmp1_file_totals' '_' dte_str '.csv'],sum_tbl);

%% Copy Files over %%%
for iS = 1:size(c12_dta_out,1)
    if ~isempty(c12_dta_out{iS,4}) && ~exist([ t1w_c12_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii'],'file') 
        copyfile(c12_dta_out{iS,4}, [ t1w_c12_dir '/' c12_dta_out{iS,1} '_' 'cat12' '_' 'mwp1.nii']);
    end
end
