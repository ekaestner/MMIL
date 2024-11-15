

%% Redcap Codes
% Variables of interest
cov_var = { 'site' 'dx' 'sex' 'age' 'handedness' 'localization' 'lateralization' 'mts'};
sbj_nme_var = 'subjid';
sbj_nme_org = 'old_subjid';

% Put together codebook
cde_bok.dx = { 0 'HC' ; ...
               1 'EPD' ; ...
               9 'other' };

cde_bok.sex = { 1 'Male' ; ...
                2 'Female' ; ...
                3 'other' ; ...
                4 'unknown'};

cde_bok.handedness = { 1 'Right' ; ...
                       2 'Left' ; ...
                       3 'Ambi' ; ...
                       4 'unknown'};

cde_bok.localization = { 1 'temporal' ; ...
                         2 'temporal_plus' ; ...
                         3 'ex-tle-frontal' ; ...
                         4 'unknown'};

cde_bok.lateralization = { 1 'left' ; ...
                           2 'right' ; ...
                           3 'bilateral' ; ...
                           4 'multifocal' ; ...
                           9 'unknown' };

cde_bok.mts = { 0 'noHS' ; ...
                1 'HS_mri' ; ...
                2 'HS_pathology' ; ...
                9 'unknown' };

%% Load ENIGMA Redcap
fcfg = [];
fcfg.dta_loc = [red_cap_dir '/' red_cap_fle];
[ red_cap_dta, red_cap_sbj, red_cap_col ] = ejk_dta_frm(fcfg);

%% Load Imaging Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' 'CAT12_wmp1_file_locations' '_' dte_str '.csv'];
[ img_dta, img_sbj, img_col ] = ejk_dta_frm(fcfg);

%% Load and Format ENIGMA Redcap


reg_dta_row = [ 3:9 14:19];
img_dta_row = [ 1:2 10:13 ];

cmb_sbj = unique(red_cap_sbj);
cmb_col = red_cap_col;
cmb_dta = cell( numel(cmb_sbj), numel(red_cap_col) );
for iS = 1:numel(cmb_sbj)

    sbj_row = find(strcmpi(red_cap_sbj,cmb_sbj{iS}));
    if numel(sbj_row)>2; error('rows greater than 2'); end
    reg_row = sbj_row( ~strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );
    img_row = sbj_row( strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );

    cmb_dta(iS,reg_dta_row) = red_cap_dta(reg_row,reg_dta_row);
    if ~isempty(img_row); cmb_dta(iS,img_dta_row) = red_cap_dta(img_row,img_dta_row); end

    % fix ses BS
    ses_pos = string_find( cmb_dta(iS,strcmpi(red_cap_col,'subjid')), 'ses' );
    if ~isempty(ses_pos)
        ses_pos = strfind( cmb_dta{iS,strcmpi(red_cap_col,'subjid')}, 'ses' );
        cmb_dta{iS,strcmpi(red_cap_col,'subjid')}(ses_pos-1:end) = [];
    end

end

%% Put together filenames & characteristics
cov_img_dta = cell(size(cmb_dta,1),4+numel(cov_var));
cov_img_col = [ {'sbj_nme'} {'img_nme'} {'t1w_c12_fle'} {'t1w_mca_fle'} cov_var ];

t1w_c12_fle = dir(t1w_c12_dir); t1w_c12_fle = {t1w_c12_fle(:).name}; t1w_c12_fle = t1w_c12_fle(string_find(t1w_c12_fle,".nii"));
t1w_mca_fle = dir(t1w_mca_dir); t1w_mca_fle = {t1w_mca_fle(:).name}; t1w_mca_fle = t1w_mca_fle(string_find(t1w_mca_fle,".nii"));

[~, cov_col ] = ismember(cov_var,cmb_col);

for iS = 1:size(cmb_dta,1)

    cov_img_dta{iS,1} = cmb_dta{iS,strcmpi(cmb_col,sbj_nme_org)};
    cov_img_dta{iS,2} = cmb_dta{iS,strcmpi(cmb_col,sbj_nme_var)};
    try cov_img_dta{iS,3} = t1w_c12_fle{string_find( t1w_c12_fle, cov_img_dta{iS,2} )}; catch cov_img_dta{iS,3} = ''; end
    try cov_img_dta{iS,4} = t1w_mca_fle{string_find( t1w_mca_fle, cov_img_dta{iS,2} )}; catch cov_img_dta{iS,4} = ''; end
    cov_img_dta(iS,5:end) = cmb_dta(iS,cov_col);

end

% Find Missing
mss_fle = sort(setxor(t1w_c12_fle, unique(cov_img_dta_rcd(:,3)) ));
cell2csv( [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'missing_enigma_imaging_and_covariates' '_' dte_str '.csv'],mss_fle);

% Fix covariates
cov_nme = fieldnames(cde_bok);
cov_val = cell(1,numel(cov_nme));
for iCV = 1:numel(cov_nme)
    cov_val{iCV} = cde_bok.(cov_nme{iCV}); %cov_val{strcmpi(cov_img_col,cov_nme{iCV})} = cde_bok.(cov_nme{iCV});
end

fcfg = [];

fcfg.dta     = cov_img_dta;
fcfg.dta_col = cov_img_col;

fcfg.rcd_nme = cov_nme;
fcfg.rcd_val = cov_val;

fcfg.swt_val = 1;

cov_img_dta_rcd = ejk_recode(fcfg);

%
[ sum(~strcmpi(cov_img_dta_rcd(:,3),'')) numel(t1w_c12_fle); 
  sum(~strcmpi(cov_img_dta_rcd(:,4),'')) numel(t1w_mca_fle) ]

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'dx'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'localization'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'lateralization'))))

tabulate(cov_img_dta_rcd(:,find(strcmpi(cov_img_col,'mts'))))

cell2csv( [ prj_dir '/' 'new_data' '/' 'enigma_conglom' '/' 'enigma_imaging_and_covariates' '_' dte_str '.csv'],[ {'sbj_nme'} cov_img_col ; cmb_sbj cov_img_dta ])]

