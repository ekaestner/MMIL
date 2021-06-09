load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

out_put = [ prj_dir '/' prj_nme '/' 'SpecificCor'];

run_grp = { 'tle_controls_pre_3T_allSurg_all'   ...
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_all' ...
            'tle_post_3T_ATLonly_left' ...
            'tle_post_3T_ATLonly_right' };

%% 
cln_dta = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ]);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
cog_dta_col = cog_dta(1,2:end);
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

%% Cognitive
naming_correlation_cognitive

%% Clinical
naming_correlation_clinical

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% DTI
for iG = 1:numel(run_grp)
    for iN = 1:numel(dti_mse)
        for iR = 1:numel(dti_roi_use{iN})
            
            if ~(dti_roi_use{iN}(iR)==0)
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
            else
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_QC.csv'];
            end
            
            dti_dta = mmil_readtext(dti_dta_nme);
            dti_dta_col = ejk_fix_column_names(dti_dta(1,5:end));
            dti_dta_sbj = dti_dta(2:end,1);
            dti_dta     = dti_dta(2:end,5:end);

            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(dti_dta( grp.(run_grp{iG}), :));
            fcfg.lbl_one = dti_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'Raw' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
                        
        end
    end
end

%% MRI
for iG = 1:numel(run_grp)
    for iN = 1:numel(mri_mse)
        for iR = 1:numel(mri_roi_use{iN})
            
            if ~(mri_roi_use{iN}(iR)==0) && ~(mri_roi_use{iN}(iR)==-1)
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
            else
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_QC.csv'];
            end
            
            mri_dta = mmil_readtext(mri_dta_nme);
            mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
            mri_dta_sbj = mri_dta(2:end,1);
            mri_dta     = mri_dta(2:end,5:end);

            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(mri_dta( grp.(run_grp{iG}), :));
            fcfg.lbl_one = mri_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'Raw' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
                        
        end
    end
end

