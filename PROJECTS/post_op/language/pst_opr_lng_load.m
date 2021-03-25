% Calculate Change Scores
fcfg = [];
fcfg.rci = 1;
[ pst_cog_dta , ~ ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

% Cognitive Save
cog_dta_out(:,1) = pst_cog_dta.sbj_nme;
for iC = 1:numel(cog_tst_nme)/2
        cog_dta_out(:,iC+1) = num2cell(sbj_cog.(cog_tst_nme{iC}));
        cog_dta_out(:,iC+5) = num2cell(pst_cog_dta.(cog_tst_nme{iC+4}));
end
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv'], [ 'SubjID' cog_tst_nme ; cog_dta_out ]);

% Load MRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(mri_mse)
    for iR = 1:numel(mri_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = mri_fle;
        fcfg.mes_nme = mri_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if mri_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{mri_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{mri_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        mri_dta = ejk_extract_mmps_roi( fcfg );
        mri_dta(1,:) = cellfun(@(x) strrep(x,[mri_mse{iM} '_'],''),mri_dta(1,:),'uni',0);
        
        % Asymmetry score
        [ dta_out , dta_lbl ] = ejk_create_laterality_index( mri_dta(2:end,:) , mri_dta(1,:));
        mri_dta_lat = [ mri_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
        
        % Save
        if mri_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iM} '.csv'], mri_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iM} '_LI.csv'], mri_dta_lat)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iM} '_' mmil_spec_char(roi_nme{mri_roi_use{iM}(iR)},{'.'}) '.csv'], mri_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iM} '_' mmil_spec_char(roi_nme{mri_roi_use{iM}(iR)},{'.'}) '_LI.csv'], mri_dta_lat)
        end
        
        % QC
        if mri_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI' '/' mri_mse{iM} '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI' '/' mri_mse{iM} '_LI' '/'];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI' '/' mri_mse{iM}  '_' mmil_spec_char(roi_nme{mri_roi_use{iM}(iR)},{'.'}) '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'MRI' '/' mri_mse{iM}  '_' mmil_spec_char(roi_nme{mri_roi_use{iM}(iR)},{'.'}) '_LI' '/'];
        end

        fcfg = [];
        fcfg.sbj_nme = mri_dta(2:end,1);
        fcfg.dta     = cell2mat(mri_dta(2:end,5:end));
        fcfg.dta_lbl = mri_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'MRI';
        ejk_qc_roi(fcfg)
        
        fcfg = [];
        fcfg.sbj_nme = mri_dta_lat(2:end,1);
        fcfg.dta     = cell2mat(mri_dta_lat(2:end,5:end));
        fcfg.dta_lbl = mri_dta_lat(1,5:end);
        fcfg.out_dir = lat_dir;
        fcfg.out_pre_fix = 'MRI_lat';
        ejk_qc_roi(fcfg)
        
    end
end

% Load DTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(dti_mse)
    for iR = 1:numel(dti_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = dti_fle;
        fcfg.mes_nme = dti_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if dti_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{dti_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{dti_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        dti_dta = ejk_extract_mmps_roi( fcfg );
        dti_dta(1,:) = cellfun(@(x) strrep(x,[dti_mse{iM} '_'],''),dti_dta(1,:),'uni',0);
        
        % Asymmetry score
        [ dta_out , dta_lbl ] = ejk_create_laterality_index( dti_dta(2:end,:) , dti_dta(1,:));
        dti_dta_lat = [ dti_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
        
        % Save
        if dti_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iM} '.csv'], dti_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iM} '_LI.csv'], dti_dta_lat)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iM} '_' mmil_spec_char(roi_nme{dti_roi_use{iM}(iR)},{'.'}) '.csv'], dti_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iM} '_' mmil_spec_char(roi_nme{dti_roi_use{iM}(iR)},{'.'}) '_LI.csv'], dti_dta_lat)
        end
        
        % QC
        if dti_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI' '/' dti_mse{iM} '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI' '/' dti_mse{iM} '_LI' '/'];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI' '/' dti_mse{iM}  '_' mmil_spec_char(roi_nme{dti_roi_use{iM}(iR)},{'.'}) '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'DTI' '/' dti_mse{iM}  '_' mmil_spec_char(roi_nme{dti_roi_use{iM}(iR)},{'.'}) '_LI' '/'];
        end

        fcfg = [];
        fcfg.sbj_nme = dti_dta(2:end,1);
        fcfg.dta     = cell2mat(dti_dta(2:end,5:end));
        fcfg.dta_lbl = dti_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'DTI';
        ejk_qc_roi(fcfg)
        
        fcfg = [];
        fcfg.sbj_nme = dti_dta_lat(2:end,1);
        fcfg.dta     = cell2mat(dti_dta_lat(2:end,5:end));
        fcfg.dta_lbl = dti_dta_lat(1,5:end);
        fcfg.out_dir = lat_dir;
        fcfg.out_pre_fix = 'DTI_lat';
        ejk_qc_roi(fcfg)
        
    end
end

% Load fMRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(fmr_mse)
    for iR = 1:numel(fmr_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = fmr_fle;
        fcfg.mes_nme = fmr_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if fmr_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{fmr_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{fmr_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        fmr_dta = ejk_extract_mmps_roi( fcfg );
        fmr_dta(1,:) = cellfun(@(x) strrep(x,[fmr_mse{iM} '_'],''),fmr_dta(1,:),'uni',0);
        
        % Asymmetry score
        [ dta_out , dta_lbl ] = ejk_create_laterality_index( fmr_dta(2:end,:) , fmr_dta(1,:));
        fmr_dta_lat = [ fmr_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
        
        % Save
        if fmr_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iM} '.csv'], fmr_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iM} '_LI.csv'], fmr_dta_lat)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iM} '_' mmil_spec_char(roi_nme{fmr_roi_use{iM}(iR)},{'.'}) '.csv'], fmr_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iM} '_' mmil_spec_char(roi_nme{fmr_roi_use{iM}(iR)},{'.'}) '_LI.csv'], fmr_dta_lat)
        end
        
        % QC
        if fmr_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'fMRI' '/' fmr_mse{iM} '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'fMRI' '/' fmr_mse{iM} '_LI' '/'];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'fMRI' '/' fmr_mse{iM}  '_' mmil_spec_char(roi_nme{fmr_roi_use{iM}(iR)},{'.'}) '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'fMRI' '/' fmr_mse{iM}  '_' mmil_spec_char(roi_nme{fmr_roi_use{iM}(iR)},{'.'}) '_LI' '/'];
        end

        fcfg = [];
        fcfg.sbj_nme = fmr_dta(2:end,1);
        fcfg.dta     = cell2mat(fmr_dta(2:end,5:end));
        fcfg.dta_lbl = fmr_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'fMRI';
        ejk_qc_roi(fcfg)
        
        fcfg = [];
        fcfg.sbj_nme = fmr_dta_lat(2:end,1);
        fcfg.dta     = cell2mat(fmr_dta_lat(2:end,5:end));
        fcfg.dta_lbl = fmr_dta_lat(1,5:end);
        fcfg.out_dir = lat_dir;
        fcfg.out_pre_fix = 'fMRI_lat';
        ejk_qc_roi(fcfg)
        
    end
end

% Load rs-fMRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iM = 1:numel(rsf_mse)
    for iR = 1:numel(rsf_roi_use{iM})
        
        fcfg = [];
        fcfg.sbj_nme = sbj_nme;
        fcfg.fle_nme = rsf_fle;
        fcfg.mes_nme = rsf_mse{iM};
        fcfg.rcn_nme = rcn_fle;
        if rsf_roi_use{iM}(iR)==0
            fcfg.roi_nme = [];
        else
            fcfg.roi_nme = { [ roi_loc '/' 'lh.' roi_nme{rsf_roi_use{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{rsf_roi_use{iM}(iR)}] };
            fcfg.roi_hms = { 'lh' 'rh' };
        end
        
        rsf_dta = ejk_extract_mmps_roi( fcfg );
        rsf_dta(1,:) = cellfun(@(x) strrep(x,[rsf_mse{iM} '_'],''),rsf_dta(1,:),'uni',0);
        
        % Asymmetry score
        [ dta_out , dta_lbl ] = ejk_create_laterality_index( rsf_dta(2:end,:) , rsf_dta(1,:));
        rsf_dta_lat = [ rsf_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
        
        % Save
        if rsf_roi_use{iM}(iR)==0
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iM} '.csv'], rsf_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iM} '_LI.csv'], rsf_dta_lat)
        else
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iM} '_' mmil_spec_char(roi_nme{rsf_roi_use{iM}(iR)},{'.'}) '.csv'], rsf_dta)
            cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iM} '_' mmil_spec_char(roi_nme{rsf_roi_use{iM}(iR)},{'.'}) '_LI.csv'], rsf_dta_lat)
        end
        
        % QC
        if rsf_roi_use{iM}(iR)==0
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'rsfMRI' '/' rsf_mse{iM} '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'rsfMRI' '/' rsf_mse{iM} '_LI' '/'];
        else
            roi_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'rsfMRI' '/' rsf_mse{iM}  '_' mmil_spec_char(roi_nme{rsf_roi_use{iM}(iR)},{'.'}) '/' ];
            lat_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'rsfMRI' '/' rsf_mse{iM}  '_' mmil_spec_char(roi_nme{rsf_roi_use{iM}(iR)},{'.'}) '_LI' '/'];
        end

        fcfg = [];
        fcfg.sbj_nme = rsf_dta(2:end,1);
        fcfg.dta     = cell2mat(rsf_dta(2:end,5:end));
        fcfg.dta_lbl = rsf_dta(1,5:end);
        fcfg.out_dir = roi_dir;
        fcfg.out_pre_fix = 'rsfMRI';
        ejk_qc_roi(fcfg)
        
        fcfg = [];
        fcfg.sbj_nme = rsf_dta_lat(2:end,1);
        fcfg.dta     = cell2mat(rsf_dta_lat(2:end,5:end));
        fcfg.dta_lbl = rsf_dta_lat(1,5:end);
        fcfg.out_dir = lat_dir;
        fcfg.out_pre_fix = 'rsfMRI_lat';
        ejk_qc_roi(fcfg)
        
    end
end

% Save Demographics & Cognitive Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alena Dominance
aln_dom = mmil_readtext(aln_dom_fle);

dom_add = cell(numel(sbj_dem.sbj_nme),1);
for iS = 1:size(dom_add,1)
    sbj_row = find(strcmpi(aln_dom(:,1),sbj_dem.sbj_nme{iS}));
    if ~isempty(sbj_row)
        dom_add{iS} = aln_dom{sbj_row,3};
    else
        dom_add{iS} = '';
    end
end

% Save
cln_dta_nme = { 'SubjID'      'VisitID'           'SideOfSeizureFocus' 'FieldStrength' 'AgeAtSurgery' ...
                'Educ'        'AgeOfSeizureOnset' 'NumAEDs'            'SeizureFreq'   'SurgicalSide' ...
                'SurgeryType' 'Sex'               'Handedness'         'MTS'           'EngelOutcome' 'LangDominance' };
cln_dta_out = [ sbj_dem.sbj_nme           dti_dta(2:end,2)              sbj_sze.sbj_sde_ons           dti_dta(2:end,3)              num2cell(sbj_srg.srg_age) ...
                num2cell(sbj_dem.sbj_edu) num2cell(sbj_sze.sbj_age_ons) num2cell(sbj_sze.sbj_aed_num) num2cell(sbj_sze.sbj_sze_frq) sbj_srg.srg_sde ...
                sbj_srg.srg_typ           sbj_dem.sbj_sex               sbj_dem.sbj_hnd               sbj_sze.sbj_mts               sbj_srg.eng_out           dom_add ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical' '.csv'], [ cln_dta_nme ; cln_dta_out ]);






