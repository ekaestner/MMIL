%% Get participants
% Load Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = 1:numel(cog_tst_nme)
    sbj_nme_hld = [ sbj_nme_hld ; sbj_cog.sbj_nme(~isnan(sbj_cog.(cog_tst_nme{iT}))) ];
end
sbj_nme = unique(sbj_nme_hld);

clear sbj_cog

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

%%
ejk_chk_dir([prj_dir '/' prj_nme '/' 'Data' '/'])

% Calculate Change Scores
fcfg = [];
fcfg.rci = 1;
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

% Cognitive Save
cog_dta_out(:,1) = pst_cog_dta.sbj_nme;
for iC = 1:numel(cog_tst_nme)/2
        cog_dta_out(:,iC+1) = num2cell(sbj_cog.(cog_tst_nme{iC}));
        cog_dta_out(:,iC+(numel(cog_tst_nme)/2+1)) = num2cell(pst_cog_dta.(cog_tst_nme{iC+(numel(cog_tst_nme)/2)}));
end
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv'], [ 'SubjID' cog_tst_nme ; cog_dta_out ]);

fcfg = [];
fcfg.sbj_nme = cog_dta_out(:,1);
fcfg.dta     = cell2mat(cog_dta_out(:,2:end));
fcfg.dta_lbl = cog_tst_nme;
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'Cognitive' '/' ];
fcfg.out_pre_fix = 'Cognitive';
ejk_qc_roi(fcfg)

%% Load MRI - 238 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(mri_238_mse)
    for iR = 1:numel(mri_238_roi_use{iM})
        
        if mri_238_roi_use{iM}(iR)>=0
            
            fcfg = [];
            fcfg.sbj_nme = sbj_nme;
            fcfg.fle_nme = mri_238_fle;
            fcfg.mes_nme = mri_238_mse{iM};
            fcfg.rcn_nme = rcn_fle;
            if mri_238_roi_use{iM}(iR)==0
                fcfg.roi_nme = [];
            else
                fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{mri_238_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{mri_238_roi_use{iM}(iR)}] };
                fcfg.roi_hms = { 'lh' 'rh' };
            end
            
            mri_238_dta = ejk_extract_mmps_roi( fcfg );
            mri_238_dta(1,:) = cellfun(@(x) strrep(x,[mri_238_mse{iM} '_'],''),mri_238_dta(1,:),'uni',0);
            
        else
            
            dta_inp = mmil_readtext( [prj_dir '/' prj_nme '/' 'Data' '/' mri_238_mse{iM-1} '_238.csv'] );
            
            fcfg = [];
            
            fcfg.dta     = dta_inp;
            fcfg.cor_col = 'IntracranialVolume';
            
            mri_238_dta = ejk_cor_roi( fcfg );
            
        end
        
        % Save
        if mri_238_roi_use{iM}(iR)==0 || mri_238_roi_use{iM}(iR)==-1
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_238_mse{iM} '_238.csv'], mri_238_dta)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_238_mse{iM} '_' mmil_spec_char(roi_nme{mri_238_roi_use{iM}(iR)},{'.'}) '_238.csv'], mri_238_dta)
        end
        
        % QC
        if mri_238_roi_use{iM}(iR)==0 || mri_238_roi_use{iM}(iR)==-1
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI_238' '/' mri_238_mse{iM} '/' ];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI_238' '/' mri_238_mse{iM}  '_' mmil_spec_char(roi_nme{mri_238_roi_use{iM}(iR)},{'.'}) '/' ];
        end
        
        fcfg = [];
        fcfg.sbj_nme = mri_238_dta(2:end,1);
        fcfg.dta     = cell2mat(mri_238_dta(2:end,5:end));
        fcfg.dta_lbl = mri_238_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'MRI';
        ejk_qc_roi(fcfg)
        
    end
end

%% Load DTI - 238 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(dti_238_mse)
    for iR = 1:numel(dti_238_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = dti_238_fle;
        fcfg.mes_nme = dti_238_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if dti_238_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{dti_238_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{dti_238_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        dti_238_dta = ejk_extract_mmps_roi( fcfg );
        dti_238_dta(1,:) = cellfun(@(x) strrep(x,[dti_238_mse{iM} '_'],''),dti_238_dta(1,:),'uni',0);
                
        % Save
        if dti_238_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_238_mse{iM} '_238.csv'], dti_238_dta)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_238_mse{iM} '_' mmil_spec_char(roi_nme{dti_238_roi_use{iM}(iR)},{'.'}) '_238.csv'], dti_238_dta)
        end
        
        % QC
        if dti_238_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI_238' '/' dti_238_mse{iM} '/' ];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI_238' '/' dti_238_mse{iM}  '_' mmil_spec_char(roi_nme{dti_238_roi_use{iM}(iR)},{'.'}) '/' ];
        end

        fcfg = [];
        fcfg.sbj_nme = dti_238_dta(2:end,1);
        fcfg.dta     = cell2mat(dti_238_dta(2:end,5:end));
        fcfg.dta_lbl = dti_238_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'DTI';
        ejk_qc_roi(fcfg)
        
    end
end

%% Load MRI - dev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(mri_dev_mse)
    for iR = 1:numel(mri_dev_roi_use{iM})
        
        if mri_dev_roi_use{iM}(iR)>=0
            
            fcfg = [];
            fcfg.sbj_nme = sbj_nme;
            fcfg.fle_nme = mri_dev_fle;
            fcfg.mes_nme = mri_dev_mse{iM};
            fcfg.rcn_nme = rcn_fle;
            if mri_dev_roi_use{iM}(iR)==0
                fcfg.roi_nme = [];
            else
                fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{mri_dev_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{mri_dev_roi_use{iM}(iR)}] };
                fcfg.roi_hms = { 'lh' 'rh' };
            end
            
            mri_dev_dta = ejk_extract_mmps_roi( fcfg );
            mri_dev_dta(1,:) = cellfun(@(x) strrep(x,[mri_dev_mse{iM} '_'],''),mri_dev_dta(1,:),'uni',0);
            
        else
            
            dta_inp = mmil_readtext( [prj_dir '/' prj_nme '/' 'Data' '/' mri_dev_mse{iM-1} '_dev.csv'] );
            
            fcfg = [];
            
            fcfg.dta     = dta_inp;
            fcfg.cor_col = 'IntracranialVolume';
            
            mri_dev_dta = ejk_cor_roi( fcfg );
            
        end
        
        % Save
        if mri_dev_roi_use{iM}(iR)==0 || mri_dev_roi_use{iM}(iR)==-1
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_dev_mse{iM} '_dev.csv'], mri_dev_dta)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_dev_mse{iM} '_' mmil_spec_char(roi_nme{mri_dev_roi_use{iM}(iR)},{'.'}) '_dev.csv'], mri_dev_dta)
        end
        
        % QC
        if mri_dev_roi_use{iM}(iR)==0 || mri_dev_roi_use{iM}(iR)==-1
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI_dev' '/' mri_dev_mse{iM} '/' ];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI_dev' '/' mri_dev_mse{iM}  '_' mmil_spec_char(roi_nme{mri_dev_roi_use{iM}(iR)},{'.'}) '/' ];
        end
        
        fcfg = [];
        fcfg.sbj_nme = mri_dev_dta(2:end,1);
        fcfg.dta     = cell2mat(mri_dev_dta(2:end,5:end));
        fcfg.dta_lbl = mri_dev_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'MRI';
        ejk_qc_roi(fcfg)
        
    end
end

%% Load DTI - dev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(dti_dev_mse)
    for iR = 1:numel(dti_dev_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = dti_dev_fle;
        fcfg.mes_nme = dti_dev_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if dti_dev_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{dti_dev_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{dti_dev_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        dti_dev_dta = ejk_extract_mmps_roi( fcfg );
        dti_dev_dta(1,:) = cellfun(@(x) strrep(x,[dti_dev_mse{iM} '_'],''),dti_dev_dta(1,:),'uni',0);
                
        % Save
        if dti_dev_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_dev_mse{iM} '_dev.csv'], dti_dev_dta)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_dev_mse{iM} '_' mmil_spec_char(roi_nme{dti_dev_roi_use{iM}(iR)},{'.'}) '_dev.csv'], dti_dev_dta)
        end
        
        % QC
        if dti_dev_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI_dev' '/' dti_dev_mse{iM} '/' ];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI_dev' '/' dti_dev_mse{iM}  '_' mmil_spec_char(roi_nme{dti_dev_roi_use{iM}(iR)},{'.'}) '/' ];
        end

        fcfg = [];
        fcfg.sbj_nme = dti_dev_dta(2:end,1);
        fcfg.dta     = cell2mat(dti_dev_dta(2:end,5:end));
        fcfg.dta_lbl = dti_dev_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'DTI';
        ejk_qc_roi(fcfg)
        
    end
end

%% Save Demographics & Cognitive Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save
cln_dta_nme = { 'SubjID'      'VisitID'           'SideOfSeizureFocus' 'FieldStrength' 'AgeAtSurgery' ...
                'Educ'        'AgeOfSeizureOnset' 'NumAEDs'            'SeizureFreq'   'SurgicalSide' ...
                'SurgeryType' 'Sex'               'Handedness'         'MTS'           'EngelOutcome' };
cln_dta_out = [ sbj_dem.sbj_nme           dti_dev_dta(2:end,2)          sbj_sze.sbj_sde_ons           dti_dev_dta(2:end,3)          num2cell(sbj_srg.srg_age) ...
                num2cell(sbj_dem.sbj_edu) num2cell(sbj_sze.sbj_age_ons) num2cell(sbj_sze.sbj_aed_num) num2cell(sbj_sze.sbj_sze_frq) sbj_srg.srg_sde ...
                sbj_srg.srg_typ           sbj_dem.sbj_sex               sbj_dem.sbj_hnd               sbj_sze.sbj_mts               sbj_srg.eng_out          ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical' '.csv'], [ cln_dta_nme ; cln_dta_out ]);


