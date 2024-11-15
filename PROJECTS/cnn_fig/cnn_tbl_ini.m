ejk_chk_dir([ ult_t1w_c12_dir '/' 'covariates' '/' ]);

%% Load Data & Combine
tot_dta = cell(0);
tot_sbj = cell(0);
grp_idn = nan(0);
grp_idn_nme = cell(0);
for iD = 1:numel(bid_nme)
    if exist([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'])

        fcfg = [];
        fcfg.dta_loc = [ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'covariates' '_' dte_str '.csv'];
        fcfg.dta_col = 2;
        [ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

        tot_dta = [ tot_dta ; cov_dta];
        tot_sbj = [ tot_sbj ; cov_sbj];
        grp_idn = [ grp_idn ; ones(numel(cov_sbj),1)*iD];
        grp_idn_nme = [ grp_idn_nme ; repmat(bid_nme(iD),numel(cov_sbj),1)];
    end
end

%cell2csv([ prj_dir '/' 'new_data' '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} cov_col {'dataset'};  tot_sbj tot_dta grp_idn_nme])
cell2csv([ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'], [ {'sbj_nme'} cov_col {'dataset'};  tot_sbj tot_dta grp_idn_nme])


%% Combine session
tot_ses = cell(0);
for iD = 1:numel(bid_nme)
    if exist([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv'])
        ses_dta = mmil_readtext([ prj_dir '/' 'new_data' '/' bid_nme{iD} '/' 'cat12_t1w'  '/' bid_nme{iD} '_' 'subjects' '_' dte_str '.csv']);

        tot_ses = [ tot_ses ; strcat([ ult_t1w_c12_dir '/'], ses_dta(:,2)) ];
    end
end
%cell2csv([ prj_dir '/' 'new_data' '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'files' '_' dte_str '.csv'], tot_ses);
cell2csv([ ult_t1w_c12_dir '/' 'covariates'  '/' 'cat12' '_' 't1'  '_' 'files' '_' dte_str '.csv'], tot_ses);

%% Make groups
for iD = 1:numel(bid_nme)
    if any(grp_idn==iD)
        grp.tot_tbl.(bid_nme{iD}) = find(grp_idn==iD);
    end
end
grp.tot_tbl.total = [1:numel(grp_idn)]';
fld_nme = fieldnames(grp.tot_tbl);

%% Initialize
tbl_typ = { 'count'    'Diagnosis'               'HC/EPD/other'            1 ; ...
            'count'    'Sex'                     'Male/Female/other'       1 ; ...
            'mean/std' 'Age'                     ''                          1 ; ... 
            'count'    'Handedness'              'Right/Left/Ambi'         1 ; ...
            'count'    'Focal'                   'focal/generalized/mixed' 1; ...
            'count'    'Localization'            'temporal/temporal_plus/ex-tle-frontal/unknown'   1 ; ...
            'count'    'Lateralization'          'left/right/bilateral/multifocal/unknown'         1 ; ...
            'count'    'MTS'                     'noHS/HS_mri/HS_pathology/unknown'                1 ; ...           
            };


%% Make Table
% Create Table
fcfg = [];

for iR = 1:size(tbl_typ)
    if strcmpi(tbl_typ{iR,1},'mean/std')
        for iN = 1:numel(fld_nme)
            fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
        end
    elseif  strcmpi(tbl_typ{iR,1},'count')
        for iN = 1:numel(fld_nme)
            fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
        end
    end
end

fcfg.dta = {[ cov_col ; tot_dta ]};
fcfg.grp = grp.tot_tbl;
tbl_out  = ejk_create_table( fcfg );

% Add-ons
num_sbj{1} = 'N';
tbl_lbl{1} = '';
for iN = 1:numel(fld_nme)
    num_sbj{iN+1} = numel(grp.tot_tbl.(fld_nme{iN}));
    tbl_lbl{iN+1} = fld_nme{iN};
end
            
%% Save it
tbl_row = cell(size(tbl_typ,1),1);
for iR = 1:size(tbl_typ,1)
    if strcmpi(tbl_typ{iR,1},'count')
        tbl_row{iR,1} = [ tbl_typ{iR,2} ' (' tbl_typ{iR,3} ')'];
    elseif strcmpi(tbl_typ{iR,1},'mean/std')
        tbl_row{iR,1} = [ tbl_typ{iR,2} ' (' tbl_typ{iR,1} ')'];
    end
end

% Out
tbl_out = [ tbl_lbl; num_sbj; tbl_row tbl_out ];

%cell2csv([ prj_dir '/' 'new_data' '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates_table' '_' dte_str '.csv'],tbl_out)
cell2csv([ ult_t1w_c12_dir '/' 'covariates'  '/' 'cat12' '_' 't1'  '_' 'covariates_table' '_' dte_str '.csv'],tbl_out)




