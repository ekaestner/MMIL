clear; clc;

eng_dir = '/space/mcdonald-syn01/1/BIDS/enigma_conglom/';

mca_drv_dir = [ eng_dir '/' 'derivatives' '/' 'micapipe_v0.2.0' '/' ];

%% fsaveragepro
eng_sbj_lst = dir(mca_drv_dir); eng_sbj_lst = { eng_sbj_lst(:).name }; eng_sbj_lst = eng_sbj_lst(string_find(eng_sbj_lst,'sub-'));

ssn_out = cell(numel(eng_sbj_lst),2);
ssn_out(:,1) = eng_sbj_lst;
for iS = 1:numel(eng_sbj_lst)
    ses_fld = dir([ mca_drv_dir '/' eng_sbj_lst{iS} '/']); ses_fld = {ses_fld(3:end).name};
    if numel(ses_fld)==1
        ssn_out{iS,2} = ses_fld{1};
    elseif any(strcmpi(ses_fld,'anat'))
        ssn_out{iS,2} = '';
    elseif any(strcmpi(ses_fld,'ses-01'))
        ssn_out{iS,2} = ses_fld{strcmpi(ses_fld,'ses-01')};
    elseif any(strcmpi(ses_fld,'ses-BL'))
        ssn_out{iS,2} = ses_fld{strcmpi(ses_fld,'ses-BL')};
    elseif isempty(ses_fld)
    
    else
        error('stop here')
    end
end

cell2csv([ '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/' '/' 'enigma_session_2024_05_10.csv'],ssn_out)



