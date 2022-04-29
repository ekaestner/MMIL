cog_sbj_cln = {  };
cog_roi_cln = {  };

cln_hld_sbj.subcort_vol_238_sbj_cln = { }; % todo
cln_hld_col.subcort_vol_238_roi_cln = { }; % todo

cln_hld_sbj.cort_thick_ctx_aparc_annot_238_sbj_cln = { }; % todo
cln_hld_col.cort_thick_ctx_aparc_annot_238_roi_cln = { }; % todo

cln_hld_sbj.wmparc_FA_wm_aparc_annot_238_sbj_cln = { 'epd046'    'epd_ucsf007' 'epd_ucsf018' 'epd_ucsf026' 'fc041'};
cln_hld_col.wmparc_FA_wm_aparc_annot_238_roi_cln = { {'all_roi'} {'all_roi'}   {'all_roi'}   {'all_roi'}   {'all_roi'} };

cln_hld_sbj.fiber_FA_238_sbj_cln = { 'epd046'    'epd_ucsf007' 'epd_ucsf018' 'epd_ucsf026' 'fc041'};
cln_hld_col.fiber_FA_238_roi_cln = { {'all_roi'} {'all_roi'}   {'all_roi'}   {'all_roi'}   {'all_roi'} };

cln_hld_sbj.subcort_vol_dev_sbj_cln = { }; % todo
cln_hld_col.subcort_vol_dev_roi_cln = { }; % todo

cln_hld_sbj.cort_thick_ctx_aparc_annot_dev_sbj_cln = { }; % todo
cln_hld_col.cort_thick_ctx_aparc_annot_dev_roi_cln = { }; % todo

cln_hld_sbj.wmparc_FA_wm_aparc_annot_dev_sbj_cln = { 'epd046'    'epd_ucsf007' 'epd_ucsf018' 'epd_ucsf026' 'fc041'};
cln_hld_col.wmparc_FA_wm_aparc_annot_dev_roi_cln = { {'all_roi'} {'all_roi'}   {'all_roi'}   {'all_roi'}   {'all_roi'} };

cln_hld_sbj.fiber_FA_dev_sbj_cln = { 'epd046'    'epd_ucsf007' 'epd_ucsf018' 'epd_ucsf026' 'fc041'};
cln_hld_col.fiber_FA_dev_roi_cln = { {'all_roi'} {'all_roi'}   {'all_roi'}   {'all_roi'}   {'all_roi'} };

% DTI subjects to check: epd005, epd041, epd052, epd_ucsf035, fc023, fc028, fc043, fc052 

cln_sbj_nme = fieldnames(cln_hld_sbj);
cln_sbj_col = fieldnames(cln_hld_col);

%% Cognitive Data
cog_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '.csv']);

fcfg = [];

fcfg.dta     = cog_dta;
fcfg.dta_col = 2:size(cog_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = cog_sbj_cln;
fcfg.roi_nme = cog_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '_' 'QC' '.csv'], cln_dta)

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1);
fcfg.dta     = cell2mat(cln_dta(2:end,2:4));
fcfg.dta_lbl = cln_dta(1,2:4);
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'Cognitive' '/'];
fcfg.out_pre_fix = 'Cognitive';
ejk_qc_roi(fcfg)

%% Neurobio
inp_fle     = { mri_238_fle dti_238_fle mri_dev_fle dti_dev_fle };
inp_pre     = { mri_238_pre dti_238_pre mri_dev_pre dti_dev_pre };
inp_suf     = { mri_238_suf dti_238_suf mri_dev_suf dti_dev_suf };
inp_mse     = { mri_238_mse dti_238_mse mri_dev_mse dti_dev_mse };
inp_roi     = { mri_238_roi dti_238_roi mri_dev_roi dti_dev_roi };
inp_lat     = { mri_238_lat dti_238_lat mri_dev_lat dti_dev_lat };
inp_icv     = { mri_238_icv dti_238_icv mri_dev_icv dti_dev_icv };

lod_cnt = 1;
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            % Normal Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                nme_hld = [inp_mse{iF}{iM} '_' inp_suf{iF}{iM}];
            else
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                nme_hld = [inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM}];
            end
                        
            neu_dta = mmil_readtext(roi_fle_nme);
            
            fcfg = [];            
            fcfg.dta     = neu_dta;
            fcfg.dta_col = 5:size(neu_dta,2);
            fcfg.sbj_col = 1;            
            fcfg.sbj_nme = cln_hld_sbj.(cln_sbj_nme{string_find(cln_sbj_nme,nme_hld)});
            fcfg.roi_nme = cln_hld_col.(cln_sbj_col{string_find(cln_sbj_col,nme_hld)});         
            cln_dta = ejk_clean_roi(fcfg);
            
            cell2csv( [roi_fle_nme(1:end-4) '_' 'QC' '.csv'], cln_dta)
            
%             fcfg = [];
%             fcfg.sbj_nme = neu_dta_sbj;
%             fcfg.dta     = cell2mat(cln_dta(2:end,5:end));
%             fcfg.dta_lbl = cln_dta(1,5:end);
%             fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'fiber_FA' '/'];
%             fcfg.out_pre_fix = 'fiber_FA';
%             ejk_qc_roi(fcfg)
            
            lod_cnt = lod_cnt + 1;
            
            % ICV Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmpi(inp_icv{iF}{iM},'')
                
                fle_nme = [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '.csv'];
                
                neu_dta = mmil_readtext(fle_nme);
                
                fcfg = [];
                fcfg.dta     = neu_dta;
                fcfg.dta_col = 5:size(neu_dta,2);
                fcfg.sbj_col = 1;
                fcfg.sbj_nme = cln_hld_sbj.(cln_sbj_nme{string_find(cln_sbj_nme,nme_hld)});
                fcfg.roi_nme = cln_hld_col.(cln_sbj_col{string_find(cln_sbj_col,nme_hld)});
                cln_dta = ejk_clean_roi(fcfg);
            
                cell2csv( [fle_nme(1:end-4) '_' 'QC' '.csv'], cln_dta)
                
                lod_cnt = lod_cnt + 1;
            end
            
            % Laterality Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if inp_lat{iF}(iM)==1
                fle_nme = [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '.csv'];
                
                neu_dta = mmil_readtext(fle_nme);
                
                fcfg = [];
                fcfg.dta     = neu_dta;
                fcfg.dta_col = 5:size(neu_dta,2);
                fcfg.sbj_col = 1;
                fcfg.sbj_nme = cln_hld_sbj.(cln_sbj_nme{string_find(cln_sbj_nme,nme_hld)});
                fcfg.roi_nme = cln_hld_col.(cln_sbj_col{string_find(cln_sbj_col,nme_hld)});
                cln_dta = ejk_clean_roi(fcfg);
            
                cell2csv( [fle_nme(1:end-4) '_' 'QC' '.csv'], cln_dta)                
                
                lod_cnt = lod_cnt + 1;
            end
            
        end
    end
end

