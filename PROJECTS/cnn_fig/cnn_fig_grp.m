%%
fcfg = [];
fcfg.dta_loc = [ ult_t1w_c12_dir '/' 'covariates' '/' 'cat12' '_' 't1'  '_' 'covariates' '_' dte_str '.csv'];
fcfg.dta_col = 2;
[ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

%%
ste_var = 'Site';           ste_var_col = strcmpi(cov_col,ste_var); 
con_var = 'Diagnosis';      con_var_col = strcmpi(cov_col,con_var); 
lat_var = 'Lateralization'; lat_var_col = strcmpi(cov_col,lat_var); 
dta_set_var = 'dataset';    dta_set_col = strcmpi(cov_col,dta_set_var); 
ste_var     = 'Site';           ste_var_col = strcmpi(cov_col,ste_var); 

%%
epd_sbj = strcmpi(cov_dta(:,con_var_col),'EPD');
con_sbj = strcmpi(cov_dta(:,con_var_col),'HC');

lft_sbj = strcmpi(cov_dta(:,lat_var_col),'left');
rgh_sbj = strcmpi(cov_dta(:,lat_var_col),'right');

dta_set_nme = unique(cov_dta(:,dta_set_col));

ste_nme = unique(cov_dta(:,ste_var_col)); ste_nme(1) = [];

%% Total
% Diagnosis
grp.total.all.diagnosis.EPD = find(epd_sbj);
grp.total.all.diagnosis.HC  = find(con_sbj); 

% Lateralization
grp.total.all.lateralization.left  = find(lft_sbj);
grp.total.all.lateralization.right = find(rgh_sbj);
grp.total.all.lateralization.HC    = find(con_sbj); 

%% Dataset
for iD = 1:numel(dta_set_nme)

    grp.dataset.(dta_set_nme{iD}).diagnosis.EPD = find(epd_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
    grp.dataset.(dta_set_nme{iD}).diagnosis.HC  = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

    % Lateralization
    grp.dataset.(dta_set_nme{iD}).lateralization.left  = find(lft_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
    grp.dataset.(dta_set_nme{iD}).lateralization.right = find(rgh_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));
    grp.dataset.(dta_set_nme{iD}).lateralization.HC    = find(con_sbj & strcmpi(cov_dta(:,dta_set_col),dta_set_nme{iD}));

end

%% Site
for iS = 1:numel(ste_nme)
    try
    grp.site.(ste_nme{iS}).diagnosis.EPD = find(epd_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
    grp.site.(ste_nme{iS}).diagnosis.HC  = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));

    % Lateralization
    grp.site.(ste_nme{iS}).lateralization.left  = find(lft_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
    grp.site.(ste_nme{iS}).lateralization.right = find(rgh_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
    grp.site.(ste_nme{iS}).lateralization.HC    = find(con_sbj & strcmpi(cov_dta(:,ste_var_col),ste_nme{iS}));
    catch; end
end

%% Save
save([ prj_dir '/' 'misc' '/' 'group.mat' ], 'grp');

