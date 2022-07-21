%% Load Data
out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Cognitive' '/' 'Correlations' '/'];
    ejk_chk_dir(out_dir);

% Groups
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'])

% Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

% Cognitive
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%
inp_fle     = { mri_dev_fle dti_dev_fle };
inp_pre     = { mri_dev_pre dti_dev_pre };
inp_suf     = { mri_dev_suf dti_dev_suf };
inp_mse     = { mri_dev_mse dti_dev_mse };
inp_roi     = { mri_dev_roi dti_dev_roi };
inp_lat     = { mri_dev_lat dti_dev_lat };
inp_icv     = { mri_dev_icv dti_dev_icv };

%% Set-up correlations
cor_typ = {'pre' 'pst'};

cog_col.pre = string_find(cog_dta_col,'_pre');
cog_col.pst = [ string_find(cog_dta_col,'_chg') string_find(cog_dta_col,'_pct')];

grp_use.nme = { 'pre_cog' 'surgery' };
grp_use.pre = fieldnames(grp.pre_cog);
grp_use.pst = fieldnames(grp.surgery);

%% Run Cognitive/Neurobiological Variable Correlations
% Load Neurobiological
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                sve_fld = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '_' 'laterality' ];
            else
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                sve_fld = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '_' 'laterality' ];
            end
            
            % Load Neurobio
            if ~strcmpi(inp_icv{iF}{iM},'')
                fcfg = [];
                fcfg.dta_loc = [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '_ComBat.csv'];
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );
                out_fld = [];
            else
                fcfg = [];
                fcfg.dta_loc = [roi_fle_nme(1:end-4) '_ComBat.csv'];
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );
                out_fld = [];
            end
            
            % Load Neurobio Lateralization
            fcfg = [];
            fcfg.dta_loc = [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '.csv'];
            fcfg.dta_col = 5;
            [ neu_dta_lat, neu_dta_sbj_lat, neu_dta_col_lat] = ejk_dta_frm( fcfg );
            
            
            % Normal Correlation
            for iCT = 1:numel(cor_typ)
                for iG = 1:numel(grp_use.(cor_typ{iCT}))
                    
                    grp_ind = grp.(grp_use.nme{iCT}).(grp_use.(cor_typ{iCT}){iG});
                    
                    fcfg = [];
                    
                    fcfg.sbj_nme = cln_dta_sbj( grp_ind, 1);
                    
                    fcfg.dta_two = cell2mat(cog_dta( grp_ind, cog_col.(cor_typ{iCT})));
                    fcfg.lbl_two = cog_dta_col(cog_col.(cor_typ{iCT}));
                    
                    fcfg.cor_typ = 'spearman';
                    
                    fcfg.dta_one = cell2mat( neu_dta( grp_ind, :));
                    fcfg.lbl_one = neu_dta_col;
                    
                    fcfg.pvl_cut = 0.05;
                    fcfg.pvl_lib = 0.10;
                    
                    fcfg.force_plot = 0;
                    
                    fcfg.out_dir = [ out_dir '/' grp_use.nme{iCT} '_' grp_use.(cor_typ{iCT}){iG} '/' sve_fld];
                    
                    ejk_cross_cor( fcfg );
                    
                end
            end
            
            % Laterality Correlation
            for iCT = 1:numel(cor_typ)
                for iG = 1:numel(grp_use.(cor_typ{iCT}))
                    
                    grp_ind = grp.(grp_use.nme{iCT}).(grp_use.(cor_typ{iCT}){iG});
                    
                    fcfg = [];
                    
                    fcfg.sbj_nme = cln_dta_sbj( grp_ind, 1);
                    
                    fcfg.dta_two = cell2mat(cog_dta( grp_ind, cog_col.(cor_typ{iCT})));
                    fcfg.lbl_two = cog_dta_col(cog_col.(cor_typ{iCT}));
                    
                    fcfg.cor_typ = 'spearman';
                    
                    fcfg.dta_one = cell2mat( neu_dta_lat( grp_ind, :));
                    fcfg.lbl_one = neu_dta_col_lat;
                    
                    fcfg.pvl_cut = 0.05;
                    fcfg.pvl_lib = 0.10;
                    
                    fcfg.force_plot = 0;
                    
                    fcfg.out_dir = [ out_dir '/' grp_use.nme{iCT} '_' grp_use.(cor_typ{iCT}){iG} '/' sve_fld_lat ];
                    
                    ejk_cross_cor( fcfg );
                    
                end
            end
            
        end
    end
end