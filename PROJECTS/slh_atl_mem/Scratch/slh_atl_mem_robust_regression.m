
% Final Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv' ];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] = ejk_dta_frm( fcfg );

fnl_dta(cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_pre'))) & cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_chg'))),strcmpi(fnl_dta_col,'lm2_rci')) = {NaN};

%% Table Build
tbl_typ = { 'mean/std' 'lm2_pre' '' 1; ...
            'mean/std' 'lm2_chg' '' 1 ; ...
            'mean/std' 'lm2_rci' '' 1};

load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])        
        
%% Robust Correlations
grp_nme     = { 'surgery' };
grp_nme_sub = { {'pst_cog' 'pst_cog_dti'} };

cog_var = tbl_typ(:,2);
cog_ind = find(ismember(fnl_dta_col,cog_var));

cor_ind = 1:size(fnl_dta_col,2); cor_ind = setxor(cor_ind,cog_ind); cor_ind = cor_ind(cellfun(@isnumeric,fnl_dta(1,cor_ind))); cor_ind = [ find(strcmpi(fnl_dta_col,'lm2_pre')) cor_ind ];
cor_var = fnl_dta_col(cor_ind);

for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        for iFN = 1:numel(fld_nme)
            grp_ind = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iFN});
            
            fcfg = [];
            
            fcfg.sbj_nme = fnl_dta_sbj( grp_ind, 1);
            
            fcfg.dta_two = cell2mat(fnl_dta( grp_ind, cog_ind));
            fcfg.lbl_two = strcat('cog_',fnl_dta_col(cog_ind));
            
            fcfg.dta_one = cell2mat(fnl_dta( grp_ind, cor_ind));
            fcfg.lbl_one = fnl_dta_col(cor_ind);
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.10;
            
            fcfg.force_plot = 0;
            
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' 'robust' '/' fld_nme{iFN}  '/' ];
            
            ejk_cross_cor_robust( fcfg );
        end
        
    end
end

%% Make Table


