%% Load Data
out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Cognitive' '/'];

% Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
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
% Cognitive
cor_typ = {'pre_tot' 'pst'};

cog_col.pre_tot = string_find(cog_dta_col,'_pre');
cog_col.pst     = [ string_find(cog_dta_col,'_chg') string_find(cog_dta_col,'_pct')];

% Groups
cmp_typ     = { 'QC_none' 'QC_visual_imaging' 'QC_post_imaging' 'QC_imaging_LowCog' };
grp_nme     = { 'grp.mat' 'grp_vis_qal.mat'   'grp_img_qal.mat' 'grp_img_cog_qal.mat' };
grp_run     = { [1 2]     [1 2]               [1 2]             [2]};

grp_use.nme     = { 'total_pre_cog' 'surgery' };
grp_use.sub_nme = { {'pre_cog_dti'} {'pst_cog_dti'} };

%% Run Cognitive/Neurobiological Variable Correlations
for iGU = 1:numel(cmp_typ)
    
    % Load Groups
    load([ prj_dir '/' prj_nme '/' 'Data' '/' grp_nme{iGU}])
    
    sub_out_dir = [ out_dir '/' 'Correlations' '_' cmp_typ{iGU} '/'];
    ejk_chk_dir(sub_out_dir);
    
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
                for iCT = 1:numel(grp_run{iGU})
                    fld_nme = fieldnames(grp.(grp_use.nme{grp_run{iGU}(iCT)}).(grp_use.sub_nme{grp_run{iGU}(iCT)}{1}));
                    for iG = 1:numel(fld_nme)
                        
                        grp_ind = grp.(grp_use.nme{grp_run{iGU}(iCT)}).(grp_use.sub_nme{grp_run{iGU}(iCT)}{1}).(fld_nme{iG});
                        
                        fcfg = [];
                        
                        fcfg.sbj_nme = cln_dta_sbj( grp_ind, 1);
                        
                        fcfg.dta_two = cell2mat(cog_dta( grp_ind, cog_col.(cor_typ{grp_run{iGU}(iCT)})));
                        fcfg.lbl_two = cog_dta_col(cog_col.(cor_typ{grp_run{iGU}(iCT)}));
                        
                        fcfg.cor_typ = 'spearman';
                        
                        fcfg.dta_one = cell2mat( neu_dta( grp_ind, :));
                        fcfg.lbl_one = neu_dta_col;
                        
                        fcfg.pvl_cut = 0.05;
                        fcfg.pvl_lib = 0.10;
                        
                        fcfg.force_plot = 0;
                        
                        fcfg.out_dir = [ sub_out_dir '/' grp_use.nme{grp_run{iGU}(iCT)} '_' fld_nme{iG} '/' sve_fld];
                        
                        ejk_cross_cor( fcfg );
                        
                    end
                end
                
                % Laterality Correlation
                for iCT = 1:numel(grp_run{iGU})
                    fld_nme = fieldnames(grp.(grp_use.nme{grp_run{iGU}(iCT)}).(grp_use.sub_nme{grp_run{iGU}(iCT)}{1}));
                    for iG = 1:numel(fld_nme)
                        
                        grp_ind = grp.(grp_use.nme{grp_run{iGU}(iCT)}).(grp_use.sub_nme{grp_run{iGU}(iCT)}{1}).(fld_nme{iG});
                        
                        fcfg = [];
                        
                        fcfg.sbj_nme = cln_dta_sbj( grp_ind, 1);
                        
                        fcfg.dta_two = cell2mat(cog_dta( grp_ind, cog_col.(cor_typ{grp_run{iGU}(iCT)})));
                        fcfg.lbl_two = cog_dta_col(cog_col.(cor_typ{grp_run{iGU}(iCT)}));
                        
                        fcfg.cor_typ = 'spearman';
                        
                        fcfg.dta_one = cell2mat( neu_dta_lat( grp_ind, :));
                        fcfg.lbl_one = neu_dta_col_lat;
                        
                        fcfg.pvl_cut = 0.05;
                        fcfg.pvl_lib = 0.10;
                        
                        fcfg.force_plot = 0;
                        
                        fcfg.out_dir = [ sub_out_dir '/' grp_use.nme{grp_run{iGU}(iCT)} '_' fld_nme{iG} '/' sve_fld_lat];
                        
                        ejk_cross_cor( fcfg );
                        
                    end
                end
                
            end
        end
    end
    
end

%% Table it up - LM2
grp_nme = { 'total_pre_cog_tle_hc' 'surgery_ltle_atl' 'surgery_ltle_slah' 'surgery_rtle_atl' 'surgery_rtle_slah' };
tst_nme = { 'lm2.pre'              'lm2.chg'          'lm2.chg'           'lm2.chg'          'lm2.chg' };
grp_run_ind = [ 1                  2                  2                   2                   2]; 

% Loop through neurobio
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            % Identify Neurobio Name
            if inp_roi{iF}{iM}(iR)==0
                sve_fld     = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '_' 'laterality' ];
            else
                sve_fld     = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '_' 'laterality' ];
            end
            
            % bilateral
            ref_tbl = mmil_readtext([ out_dir '/' 'Correlations' '_' cmp_typ{1} '/' '/' grp_nme{2} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
            ovr_hld_tbl = cell(size(ref_tbl,1),1+(numel(grp_nme)*numel(cmp_typ))+(1*numel(cmp_typ)));
            ovr_hld_tbl(:,1)=ref_tbl(:,1);
            for iGU = 1:numel(cmp_typ)
                
                % Load R & P
                sub_out_dir = [ out_dir '/' 'Correlations' '_' cmp_typ{iGU} '/'];
                                                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                for iG = 1:numel(grp_nme)
                    if any(ismember(grp_run{iGU},grp_run_ind(iG)))
                        
                        bil_tbl_hld{1,iG+1} = grp_nme{iG};
                        
                        cor_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
                        cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                        cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                        pvl_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_pvalues.csv' ]);
                        pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                        pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                        num_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_n.csv' ]);
                        num_hld = num_hld(2:end,strcmpi(num_hld(1,:),tst_nme{iG}));
                        num_hld(cellfun(@isstr,num_hld)) = {NaN};
                        
                        bil_tbl_hld(2:end,iG+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                    end
                end
                
                ovr_hld_tbl(:, (((iGU-1)*numel(grp_nme))+iGU)+1:((iGU*numel(grp_nme))+iGU) ) = bil_tbl_hld(:,2:end);
                
            end
            
            cell2csv([ out_dir '/' 'correlations_hld' '/' 'LM2' '/' sve_fld '.csv' ], ovr_hld_tbl)
            
            % laterality index
            ref_tbl = mmil_readtext([ out_dir '/' 'Correlations' '_' cmp_typ{1} '/' '/' grp_nme{2} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
            ovr_hld_tbl = cell(size(ref_tbl,1),1+(numel(grp_nme)*numel(cmp_typ))+(1*numel(cmp_typ)));
            ovr_hld_tbl(:,1)=ref_tbl(:,1);
            for iGU = 1:numel(cmp_typ)
                
                % Load R & P
                sub_out_dir = [ out_dir '/' 'Correlations' '_' cmp_typ{iGU} '/'];
                                                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                for iG = 1:numel(grp_nme)
                    if any(ismember(grp_run{iGU},grp_run_ind(iG)))
                        
                        bil_tbl_hld{1,iG+1} = grp_nme{iG};
                        
                        cor_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
                        cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                        cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                        pvl_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_pvalues.csv' ]);
                        pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                        pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                        num_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_n.csv' ]);
                        num_hld = num_hld(2:end,strcmpi(num_hld(1,:),tst_nme{iG}));
                        num_hld(cellfun(@isstr,num_hld)) = {NaN};
                        
                        bil_tbl_hld(2:end,iG+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                    end
                end
                
                ovr_hld_tbl(:, (((iGU-1)*numel(grp_nme))+iGU)+1:((iGU*numel(grp_nme))+iGU) ) = bil_tbl_hld(:,2:end);
                
            end
            
            cell2csv([ out_dir '/' 'correlations_hld' '/' 'LM2' '/' sve_fld_lat '.csv' ], ovr_hld_tbl)
            
        end
    end
    
    
end

%% Table it up - LM2
grp_nme = { 'total_pre_cog_tle_hc' 'surgery_ltle_atl' 'surgery_ltle_slah' 'surgery_rtle_atl' 'surgery_rtle_slah' };
tst_nme = { 'vp2.pre'              'vp2.chg'          'vp2.chg'           'vp2.chg'          'vp2.chg' };
grp_run_ind = [ 1                  2                  2                   2                   2]; 

% Loop through neurobio
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            % Identify Neurobio Name
            if inp_roi{iF}{iM}(iR)==0
                sve_fld     = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '_' 'laterality' ];
            else
                sve_fld     = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '_' 'laterality' ];
            end
            
            % bilateral
            ref_tbl = mmil_readtext([ out_dir '/' 'Correlations' '_' cmp_typ{1} '/' '/' grp_nme{2} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
            ovr_hld_tbl = cell(size(ref_tbl,1),1+(numel(grp_nme)*numel(cmp_typ))+(1*numel(cmp_typ)));
            ovr_hld_tbl(:,1)=ref_tbl(:,1);
            for iGU = 1:numel(cmp_typ)
                
                % Load R & P
                sub_out_dir = [ out_dir '/' 'Correlations' '_' cmp_typ{iGU} '/'];
                                                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                for iG = 1:numel(grp_nme)
                    if any(ismember(grp_run{iGU},grp_run_ind(iG)))
                        
                        bil_tbl_hld{1,iG+1} = grp_nme{iG};
                        
                        cor_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
                        cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                        cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                        pvl_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_pvalues.csv' ]);
                        pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                        pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                        num_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_n.csv' ]);
                        num_hld = num_hld(2:end,strcmpi(num_hld(1,:),tst_nme{iG}));
                        num_hld(cellfun(@isstr,num_hld)) = {NaN};
                        
                        bil_tbl_hld(2:end,iG+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                    end
                end
                
                ovr_hld_tbl(:, (((iGU-1)*numel(grp_nme))+iGU)+1:((iGU*numel(grp_nme))+iGU) ) = bil_tbl_hld(:,2:end);
                
            end
            
            cell2csv([ out_dir '/' 'correlations_hld' '/' 'vp2' '/' sve_fld '.csv' ], ovr_hld_tbl)
            
            % laterality index
            ref_tbl = mmil_readtext([ out_dir '/' 'Correlations' '_' cmp_typ{1} '/' '/' grp_nme{2} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
            ovr_hld_tbl = cell(size(ref_tbl,1),1+(numel(grp_nme)*numel(cmp_typ))+(1*numel(cmp_typ)));
            ovr_hld_tbl(:,1)=ref_tbl(:,1);
            for iGU = 1:numel(cmp_typ)
                
                % Load R & P
                sub_out_dir = [ out_dir '/' 'Correlations' '_' cmp_typ{iGU} '/'];
                                                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                for iG = 1:numel(grp_nme)
                    if any(ismember(grp_run{iGU},grp_run_ind(iG)))
                        
                        bil_tbl_hld{1,iG+1} = grp_nme{iG};
                        
                        cor_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
                        cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                        cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                        pvl_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_pvalues.csv' ]);
                        pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                        pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                        num_hld = mmil_readtext([ sub_out_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_n.csv' ]);
                        num_hld = num_hld(2:end,strcmpi(num_hld(1,:),tst_nme{iG}));
                        num_hld(cellfun(@isstr,num_hld)) = {NaN};
                        
                        bil_tbl_hld(2:end,iG+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                    end
                end
                
                ovr_hld_tbl(:, (((iGU-1)*numel(grp_nme))+iGU)+1:((iGU*numel(grp_nme))+iGU) ) = bil_tbl_hld(:,2:end);
                
            end
            
            cell2csv([ out_dir '/' 'correlations_hld' '/' 'vp2' '/' sve_fld_lat '.csv' ], ovr_hld_tbl)
            
        end
    end
    
    
end
    
    















