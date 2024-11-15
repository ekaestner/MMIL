dte_str = '2024_10_29';

%%
cov_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_eng/data/covariates/';
cov_fle = 'cat12_t1_covariates_2024_10_08.csv';
loc_fle = 'cat12_t1_files_2024_10_08.csv';

old_dta_loc = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig//new_data//_cat12_t1w_data_repo//';
new_dta_loc = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_eng/data/numpy_3D_minmax/';
new_dta_pkl = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_eng/data/covariates/';

old_dta_end = '.nii';
new_dta_end = '.npy';

%%
cov_col_nme_cph = { 'cov_sbj'        'ID' ; ...
                    'Age'            'age' ; ...
                    'Sex'            'sex' ; ...
                    'Diagnosis'      'diagnosis' ; ...
                    'Lateralization' 'side' ; ...
                    'Site'           'site' ; ...
                    'pre_dummy'      'pre' ; ...
                    'post_dummy'     'pos' ; ...
                    'file_loc'       'smriPath' };


%%
fcfg = [];
fcfg.dta_loc = [ cov_dir '/' cov_fle];
fcfg.dta_col = 2;
[ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ cov_dir '/' loc_fle];
fcfg.dta_col = 2;
[ loc_dta, loc_sbj, loc_col] = ejk_dta_frm( fcfg );

cov_col = [ cov_col {'file_loc'} {'pre_dummy'} {'post_dummy'} {'cov_sbj'} ];
cov_dta = [ cov_dta loc_sbj num2cell(nan(numel(cov_sbj),2)) cov_sbj];

cov_col_ind = nan(1,size(cov_col_nme_cph,1));
for iI = 1:numel(cov_col_ind)
    cov_col_ind(iI) = find(strcmpi(cov_col,cov_col_nme_cph{iI,1}));
end

%%
% Keep only HC & EPD
kep_ind = strcmpi(cov_dta(:,strcmpi(cov_col,'localization')),'temporal') + strcmpi(cov_dta(:,strcmpi(cov_col,'diagnosis')),'HC');
cov_dta = cov_dta(logical(kep_ind),:);

% Update file location for new project
cov_dta(:,strcmpi(cov_col,'file_loc')) = strrep(cov_dta(:,strcmpi(cov_col,'file_loc')),old_dta_loc,new_dta_loc);
cov_dta(:,strcmpi(cov_col,'file_loc')) = strrep(cov_dta(:,strcmpi(cov_col,'file_loc')),old_dta_end,new_dta_end);

% Recode to match Reihaneh
cde_nme = { 'Diagnosis' 'Lateralization' 'Sex'};

cde_bok{1} = { 'HC'  0 ; ...
               'EPD' 1 };

cde_bok{2} = { 'left'  1 ; ...
               'right' 2 };

cde_bok{3} = { 'Male'   0 ; ...
               'Female' 1 };

fcfg = [];

fcfg.dta     = cov_dta;
fcfg.dta_col = cov_col;

fcfg.rcd_nme = cde_nme;
fcfg.rcd_val = cde_bok;

fcfg.swt_val = 1;

cov_dta_rcd = ejk_recode(fcfg);

%%
rei_pkl_out = cov_dta_rcd( :, cov_col_ind);

cell2csv([ new_dta_pkl '/' 'pickle_prep_' dte_str '.csv'], [ cov_col_nme_cph(:,2)' ; rei_pkl_out])








