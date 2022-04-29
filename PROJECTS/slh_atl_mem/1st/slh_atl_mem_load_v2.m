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

fcfg = [];
fcfg.sbj_nme = cog_dta_out(:,1);
fcfg.dta     = cell2mat(cog_dta_out(:,2:end));
fcfg.dta_lbl = cog_tst_nme;
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' 'Cognitive' '/' ];
fcfg.out_pre_fix = 'Cognitive';
ejk_qc_roi(fcfg)

% Calculate Cognitive Categories
for iC = 2:3
    con_ind = string_find(cog_dta_out(:,1),'fc');
    
    cog_dta_out_cat(:,iC-1) = cog_dta_out(:,iC);
    
    con_men = nanmean(cell2mat(cog_dta_out_cat(con_ind,iC-1)));
    con_std = nanstd(cell2mat(cog_dta_out_cat(con_ind,iC-1)));
    
    zsc_hld = cellfun(@(x) (x-con_men) / con_std, cog_dta_out_cat(:,iC-1));
    
    cog_dta_out_cat( zsc_hld <= -1.5, iC-1) = {'Impaired'};
    cog_dta_out_cat( zsc_hld >  -1.5, iC-1) = {'NotImpaired'};
    cog_dta_out_cat( isnan(zsc_hld),  iC-1) = {NaN};
    
    cog_dta_zsc(:,iC-1) = num2cell(zsc_hld);
    
end

for iC = 4:5
    zsc_hld = cell2mat(cog_dta_out(:,iC));
    
    cog_dta_out_cat( zsc_hld <= -1.5, iC-1) = {'Impaired'};
    cog_dta_out_cat( zsc_hld >  -1.5, iC-1) = {'NotImpaired'};
    cog_dta_out_cat( isnan(zsc_hld),  iC-1) = {NaN};    
end   

out_dta_fin = [ ['SubjID' cog_tst_nme(1:2) strcat(cog_tst_nme(1:2),'_zsc') strcat(cog_tst_nme(1:2),'_cat') cog_tst_nme(3:4)   strcat(cog_tst_nme(3:4),'_cat')] ; ...
                [cog_dta_out(:,1:3)        cog_dta_zsc                     cog_dta_out_cat(:,1:2)          cog_dta_out(:,4:5) cog_dta_out_cat(:,3:4)] ];

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv'], out_dta_fin);

%% Load Neurobio
inp_fle     = { mri_238_fle dti_238_fle mri_dev_fle dti_dev_fle };
inp_pre     = { mri_238_pre dti_238_pre mri_dev_pre dti_dev_pre };
inp_suf     = { mri_238_suf dti_238_suf mri_dev_suf dti_dev_suf };
inp_mse     = { mri_238_mse dti_238_mse mri_dev_mse dti_dev_mse };
inp_roi     = { mri_238_roi dti_238_roi mri_dev_roi dti_dev_roi };
inp_lat     = { mri_238_lat dti_238_lat mri_dev_lat dti_dev_lat };
inp_icv     = { mri_238_icv dti_238_icv mri_dev_icv dti_dev_icv };

for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '/' ];
            else
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM}  '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '/' ];
            end
            
            if inp_roi{iF}{iM}(iR)==0
                roi_nme_hld = [];
                roi_hms_hld = [];
            else
                roi_nme_hld = { [ roi_loc '/' 'lh.' roi_nme{inp_roi{iF}{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{inp_roi{iF}{iM}(iR)}] };
                roi_hms_hld = { 'lh' 'rh' };
            end
            
            % Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fcfg = [];
            fcfg.sbj_nme = sbj_nme;
            fcfg.fle_nme = inp_fle{iF};
            fcfg.mes_nme = inp_mse{iF}{iM};
            fcfg.rcn_nme = rcn_fle;
            fcfg.roi_nme = roi_nme_hld;
            fcfg.roi_hms = roi_hms_hld;
            
            [ neu_bio_dta, neu_bio_wrn ] = ejk_extract_mmps_roi( fcfg );
            neu_bio_dta(1,:) = cellfun(@(x) strrep(x,[inp_mse{iF}{iM} '_'],''),neu_bio_dta(1,:),'uni',0);
            
            % Save
            cell2csv( roi_fle_nme, neu_bio_dta)
            if ~isempty(neu_bio_wrn); cell2csv( [roi_fle_nme(1:end-4) '_warnings.txt' ], neu_bio_wrn); end
            
            % QC
            fcfg = [];
            fcfg.sbj_nme = neu_bio_dta(2:end,1);
            fcfg.dta     = cell2mat(neu_bio_dta(2:end,5:end));
            fcfg.dta_lbl = neu_bio_dta(1,5:end);
            fcfg.out_dir     = qal_dir;
            fcfg.out_pre_fix = inp_pre{iF};
            ejk_qc_roi(fcfg)
            
            % ICV normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmpi(inp_icv{iF}{iM},'')
                
                fcfg = [];
                
                fcfg.dta     = neu_bio_dta;
                fcfg.cor_col = inp_icv{iF}{iM};
                
                neu_bio_dta_icv = ejk_cor_roi( fcfg );
                
                % Save
                cell2csv( [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '.csv'], neu_bio_dta_icv)
                
            end
            
            % Laterality Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if inp_lat{iF}(iM)==1
                
                [ dta_out , dta_lbl ] = ejk_create_laterality_index( neu_bio_dta(2:end,:) , neu_bio_dta(1,:));
                neu_bio_dta_lat = [ neu_bio_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
                
                % Save
                cell2csv( [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '.csv'], neu_bio_dta_lat)
                
                % QC
                fcfg = [];
                fcfg.sbj_nme = neu_bio_dta(2:end,1);
                fcfg.dta     = cell2mat(neu_bio_dta(2:end,5:end));
                fcfg.dta_lbl = neu_bio_dta(1,5:end);
                fcfg.out_dir     = [ qal_dir(1:end-1) '_' 'LateralityIndex' ];
                fcfg.out_pre_fix = inp_pre{iF};
                ejk_qc_roi(fcfg)
                
            end
            
            % Clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear neu_bio_dta neu_bio_dta_icv neu_bio_dta_lat neu_bio_wrn
            
        end
    end
end

%% Save Demographics & Cognitive Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iF = 2; iM = 1;
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = inp_fle{iF};
fcfg.mes_nme = inp_mse{iF}{iM};
fcfg.rcn_nme = rcn_fle;
fcfg.roi_nme = roi_nme_hld;
fcfg.roi_hms = roi_hms_hld;
[ neu_bio_dta, neu_bio_wrn ] = ejk_extract_mmps_roi( fcfg );

% Save
cln_dta_nme = { 'SubjID'      'VisitID'           'SideOfSeizureFocus' 'FieldStrength' 'AgeAtSurgery' ...
                'Educ'        'AgeOfSeizureOnset' 'NumAEDs'            'SeizureFreq'   'SurgicalSide' ...
                'SurgeryType' 'Sex'               'Handedness'         'MTS'           'EngelOutcome' };
cln_dta_out = [ sbj_dem.sbj_nme           neu_bio_dta(2:end,2)          sbj_sze.sbj_sde_ons           neu_bio_dta(2:end,3)          num2cell(sbj_srg.srg_age) ...
                num2cell(sbj_dem.sbj_edu) num2cell(sbj_sze.sbj_age_ons) num2cell(sbj_sze.sbj_aed_num) num2cell(sbj_sze.sbj_sze_frq) sbj_srg.srg_sde ...
                sbj_srg.srg_typ           sbj_dem.sbj_sex               sbj_dem.sbj_hnd               sbj_sze.sbj_mts               sbj_srg.eng_out          ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical' '.csv'], [ cln_dta_nme ; cln_dta_out ]);


