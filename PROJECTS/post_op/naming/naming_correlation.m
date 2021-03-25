load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

out_put = [ prj_dir '/' prj_nme '/' 'TotalCor'];

run_grp = fieldnames(grp);

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% DTI
for iG = 1:numel(run_grp)
    for iN = 1:numel(dti_mse)
        for iR = 1:numel(dti_roi_use{iN})
            
            if ~(dti_roi_use{iN}(iR)==0)
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
                dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_LI_QC.csv'];
            else
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_QC.csv'];
                dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_LI_QC.csv'];
            end
            
            dti_dta = mmil_readtext(dti_dta_nme);
            dti_dta_col = ejk_fix_column_names(dti_dta(1,5:end));
            dti_dta_sbj = dti_dta(2:end,1);
            dti_dta     = dti_dta(2:end,5:end);
            dti_dta_lat = mmil_readtext(dti_dta_lat_nme);
            dti_dta_lat_col = dti_dta_lat(1,5:end);
            dti_dta_lat_sbj = dti_dta_lat(2:end,1);
            dti_dta_lat     = dti_dta_lat(2:end,5:end);
            
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
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(dti_dta_lat( grp.(run_grp{iG}), :));
            fcfg.lbl_one = dti_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'LI' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% MRI
for iG = 1:numel(run_grp)
    for iN = 1%:numel(mri_mse)
        for iR = 1:numel(mri_roi_use{iN})
            
            if ~(mri_roi_use{iN}(iR)==0)
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
                mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_LI_QC.csv'];
            else
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_QC.csv'];
                mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_LI_QC.csv'];
            end
            
            mri_dta = mmil_readtext(mri_dta_nme);
            mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
            mri_dta_sbj = mri_dta(2:end,1);
            mri_dta     = mri_dta(2:end,5:end);
            mri_dta_lat = mmil_readtext(mri_dta_lat_nme);
            mri_dta_lat_col = mri_dta_lat(1,5:end);
            mri_dta_lat_sbj = mri_dta_lat(2:end,1);
            mri_dta_lat     = mri_dta_lat(2:end,5:end);
            
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
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(mri_dta_lat( grp.(run_grp{iG}), :));
            fcfg.lbl_one = mri_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'LI' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% fMRI
for iG = 1:numel(run_grp)
    for iN = 1:numel(fmr_mse)
        for iR = 1:numel(fmr_roi_use{iN})
            
            if ~(fmr_roi_use{iN}(iR)==0)
                fmr_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_' mmil_spec_char(roi_nme{fmr_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
                fmr_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_' mmil_spec_char(roi_nme{fmr_roi_use{iN}(iR)},{'.'}) '_LI_QC.csv'];
            else
                fmr_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_QC.csv'];
                fmr_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_LI_QC.csv'];
            end
            
            fmr_dta = mmil_readtext(fmr_dta_nme);
            fmr_dta_col = ejk_fix_column_names(fmr_dta(1,5:end));
            fmr_dta_sbj = fmr_dta(2:end,1);
            fmr_dta     = fmr_dta(2:end,5:end);
            fmr_dta_lat = mmil_readtext(fmr_dta_lat_nme);
            fmr_dta_lat_col = fmr_dta_lat(1,5:end);
            fmr_dta_lat_sbj = fmr_dta_lat(2:end,1);
            fmr_dta_lat     = fmr_dta_lat(2:end,5:end);
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(fmr_dta( grp.(run_grp{iG}), :));
            fcfg.lbl_one = fmr_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'fMRI' '/' fmr_mse{iN} '/' 'Raw' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(fmr_dta_lat( grp.(run_grp{iG}), :));
            fcfg.lbl_one = fmr_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'fMRI' '/' fmr_mse{iN} '/' 'LI' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% rsfMRI
for iG = 1:numel(run_grp)
    for iN = 1:numel(rsf_mse)
        for iR = 1:numel(rsf_roi_use{iN})
            
            if ~(rsf_roi_use{iN}(iR)==0)
                rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
                if iN ~= 3; rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '_LI_QC.csv']; end
            else
                rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_QC.csv'];
                if iN ~= 3; rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_LI_QC.csv']; end
            end
            
            rsf_dta = mmil_readtext(rsf_dta_nme);
            rsf_dta_col = cellfun(@(x) mmil_spec_char(x,{'.'}), ejk_fix_column_names(rsf_dta(1,5:end)), 'uni', 0);
            rsf_dta_sbj = rsf_dta(2:end,1);
            rsf_dta     = rsf_dta(2:end,5:end);
            if iN ~= 3
                rsf_dta_lat = mmil_readtext(rsf_dta_lat_nme);
                rsf_dta_lat_col = rsf_dta_lat(1,5:end);
                rsf_dta_lat_sbj = rsf_dta_lat(2:end,1);
                rsf_dta_lat     = rsf_dta_lat(2:end,5:end);
            end
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
            
            fcfg.dta_one = cell2mat(rsf_dta( grp.(run_grp{iG}), :));
            fcfg.lbl_one = rsf_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'Raw' '/' run_grp{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            if iN ~= 3
                fcfg = [];
                
                fcfg.sbj_nme = cln_dta_sbj( grp.(run_grp{iG}), 1);
                
                fcfg.dta_one = cell2mat(rsf_dta_lat( grp.(run_grp{iG}), :));
                fcfg.lbl_one = rsf_dta_lat_col;
                
                fcfg.cor_typ = 'spearman';
                
                fcfg.dta_two = cell2mat(cog_dta( grp.(run_grp{iG}), cog_col.(run_grp{iG})));
                fcfg.lbl_two = cog_dta_col(cog_col.(run_grp{iG}));
                
                fcfg.pvl_cut = 0.05;
                fcfg.pvl_lib = 0.20;
                
                fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'LI' '/' run_grp{iG} '/'];
                
                ejk_cross_cor( fcfg );
            end
            
        end
    end
end