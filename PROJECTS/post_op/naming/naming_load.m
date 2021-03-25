%% Get participants
% Load Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = cog_tst_inc
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
[ pst_cog_dta , ~ ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

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

%% Load MRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Load DTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Load fMRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Load rs-fMRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        if iM~=3
            [ dta_out , dta_lbl ] = ejk_create_laterality_index( rsf_dta(2:end,:) , rsf_dta(1,:));
            rsf_dta_lat = [ rsf_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
        end
        
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
        
        if iM~=3
            fcfg = [];
            fcfg.sbj_nme = rsf_dta_lat(2:end,1);
            fcfg.dta     = cell2mat(rsf_dta_lat(2:end,5:end));
            fcfg.dta_lbl = rsf_dta_lat(1,5:end);
            fcfg.out_dir = lat_dir;
            fcfg.out_pre_fix = 'rsfMRI_lat';
            ejk_qc_roi(fcfg)
        end
        
    end
end

%% Save Demographics & Cognitive Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save
cln_dta_nme = { 'SubjID'      'VisitID'           'SideOfSeizureFocus' 'FieldStrength' 'AgeAtSurgery' ...
                'Educ'        'AgeOfSeizureOnset' 'NumAEDs'            'SeizureFreq'   'SurgicalSide' ...
                'SurgeryType' 'Sex'               'Handedness'         'MTS'           'EngelOutcome' };
cln_dta_out = [ sbj_dem.sbj_nme           dti_dta(2:end,2)              sbj_sze.sbj_sde_ons           dti_dta(2:end,3)              num2cell(sbj_srg.srg_age) ...
                num2cell(sbj_dem.sbj_edu) num2cell(sbj_sze.sbj_age_ons) num2cell(sbj_sze.sbj_aed_num) num2cell(sbj_sze.sbj_sze_frq) sbj_srg.srg_sde ...
                sbj_srg.srg_typ           sbj_dem.sbj_sex               sbj_dem.sbj_hnd               sbj_sze.sbj_mts               sbj_srg.eng_out          ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical' '.csv'], [ cln_dta_nme ; cln_dta_out ]);

%% Surface load
rcn_hld = mmil_readtext(rcn_fle);

hms_hld = { 'lhs' 'rhs' };

% Put together load list %%%%%%%%%%%%%%%
for iS = 1:size(sbj_dem.sbj_nme, 1)
    rcn_row = strcmpi( rcn_hld(:,1), sbj_dem.sbj_nme{iS});
    
    sbj_srf_out{iS,1} = sbj_dem.sbj_nme{iS};
    sbj_srf_out{iS,2} = '/home/mmilmcdRSI/data/';
    sbj_srf_out{iS,3} = rcn_hld{ rcn_row, 3};
      
end

cell2csv( [ prj_dir '/' prj_nme '/' 'Data' '/' 'sbj_srf_out.csv' ], sbj_srf_out)
 
% Load data %%%%%%%%%%%%%%%
sbj_srf_out = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'sbj_srf_out.csv' ]);

% MRI
tic;
fcfg = [];

fcfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 313; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

for iH = 1:numel(hms_hld)
    srf_dta = nan(size(sbj_srf_out,1),163842);
    for iS = 1:size(size(sbj_srf_out,1),1)
        fcfg.hms         = hms_hld{iH}(1:2);
        fcfg.prc_dir     = [ sbj_srf_out{iS,2} '/' 'proc_dti' '/' ]; 
        fcfg.sbj_fsr_dir = sbj_srf_out{iS,3};
        srf_dta_hld = ejk_extract_vertices(fcfg);
        if ~isempty(srf_dta_hld); srf_dta(iS,:) = srf_dta_hld; end
    end
    srf_dta_sbj = sbj_srf_out(:,1);
    save( [ prj_dir '/' prj_nme '/' 'Data' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

% DTI
tic;
fcfg = [];

fcfg.mes_typ = 'wmparc_fa'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 313; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

for iH = 1:numel(hms_hld)
    srf_dta = nan(size(sbj_srf_out,1),163842);
    for iS = 1:size(size(sbj_srf_out,1),1)
        fcfg.hms         = hms_hld{iH}(1:2);
        fcfg.prc_dir     = [ sbj_srf_out{iS,2} '/' 'fsurf' '/' ]; 
        fcfg.sbj_fsr_dir = sbj_srf_out{iS,3};
        srf_dta_hld = ejk_extract_vertices(fcfg);
        if ~isempty(srf_dta_hld); srf_dta(iS,:) = srf_dta_hld; end
    end
    srf_dta_sbj = sbj_srf_out(:,1);
    save( [ prj_dir '/' prj_nme '/' 'Data' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc






