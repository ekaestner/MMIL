
% Groups
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'])

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

% out structure
qal_str     = [ cln_dta_sbj num2cell(zeros(numel(cln_dta_sbj),3))];
qal_str_col = { 'sbj_nme' 'dti_vis_qal' 'dti_pst_prc_qal' 'low_pre_cog' };
cell2csv([ prj_dir '/' prj_nme '/' 'Data' '/' 'QC_overall.csv'], [qal_str_col ; qal_str]);

%% Get Visual QC
% UCSD/UCSF
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = red_cap_fle;
fcfg.sep     = '|';
[~ , ~, sbj_scn, ~, ~, ~] = mmil_load_redcap(fcfg);

dti_qal_hld     = cell(numel(cln_dta_sbj),3);
dti_qal_hld_col = { 'sbj_nme' 'vis_dti' 'vis_dti_nte' };
for iS=1:numel(cln_dta_sbj)
    sbj_ind = find(strcmpi(sbj_scn.sbj_nme,cln_dta_sbj{iS}));
    dti_qal_hld{iS,1} = cln_dta_sbj{iS};
    if ~isempty(sbj_ind)
        dti_qal_hld{iS,2} = sbj_scn.dti_pre_scn_qal(sbj_ind);
        dti_qal_hld{iS,3} = sbj_scn.dti_pre_scn_qal_nte{sbj_ind};
    end    
end

% Add in Emory
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Emory' '/' 'Emory_DTI_QC_Adam.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );

% Save
cell2csv( [ prj_dir '/' prj_nme '/' 'Data' '/' 'DTI_Visual_QC.csv'], [ dti_qal_hld_col ; dti_qal_hld ]);

%% Apply QC to grp
% Load QC exclusions
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'QC_overall_ejk.csv'];
fcfg.dta_col = 2;
[ qal_dta, qal_dta_sbj, qal_dta_col] = ejk_dta_frm( fcfg );

% Apply to groups
vis_qal_rmv = find(sum(cell2mat(qal_dta(:,1)),2)>0);

fcfg = [];
fcfg.rmv_ind = vis_qal_rmv;
fcfg.grp     = grp;
vis_qal_grp  = ejk_qal_grp(fcfg);

%% Neurobio Harmonize
inp_fle     = { mri_dev_fle dti_dev_fle };
inp_pre     = { mri_dev_pre dti_dev_pre };
inp_suf     = { mri_dev_suf dti_dev_suf };
inp_mse     = { mri_dev_mse dti_dev_mse };
inp_roi     = { mri_dev_roi dti_dev_roi };
inp_lat     = { mri_dev_lat dti_dev_lat };
inp_icv     = { mri_dev_icv dti_dev_icv };

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
            
            % Load Data
            if ~strcmpi(inp_icv{iF}{iM},'')
                roi_fle_nme = [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '.csv'];
                fcfg = [];
                fcfg.dta_loc = roi_fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col, msc_out] = ejk_dta_frm( fcfg );
            else
                fcfg = [];
                fcfg.dta_loc = roi_fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col, msc_out] = ejk_dta_frm( fcfg );
            end
                                    
            % If DTI, Remove vis-QC
            if ~isempty(strfind(inp_fle{iF},'DTI_all'))
                neu_dta(vis_qal_rmv,:) = {NaN};
            end
            
            % comBAT
            btc_dta = cell(size(neu_dta_sbj,1),1);
            btc_dta(grp.site.pre_cog.control) = {'UCSD'};
            btc_dta(grp.site.pre_cog.emory) = {'Emory'};
            btc_dta(grp.site.pre_cog.ucsf) = {'UCSF'};
            btc_dta(grp.site.pre_cog.ucsd) = {'UCSD'};
            
            cov_dta          = cell(size(neu_dta_sbj,1),1);
            cov_dta(grp.diagnosis.pre_cog_dti.ltle) = {'LTLE'};
            cov_dta(grp.diagnosis.pre_cog_dti.rtle) = {'RTLE'};
            cov_dta(grp.diagnosis.pre_cog_dti.control) = {'HC'};
%             cov_dta(cellfun(@isempty,cov_dta) & ~cellfun(@isnan,neu_dta(:,5))) = {'BTLE'};
            cov_dta(cellfun(@isempty,cov_dta)) = {'empty'};
            
            % ttt = cellfun(@isempty,cov_dta) & ~cellfun(@isnan,neu_dta(:,5));
            % ttt = find(cellfun(@isempty,cov_dta))
            %[ neu_dta_sbj(ttt) btc_dta(ttt) cln_dta(ttt,2) neu_dta(ttt,40) ]
            
            fcfg = [];            
            fcfg.sbj_nme = neu_dta_sbj;            
            fcfg.dta     = cell2mat(neu_dta);
            fcfg.dta_nme = neu_dta_col;            
            fcfg.btc     = btc_dta;
            fcfg.btc_nme = {'Site'};            
            fcfg.cov     = cov_dta;
            fcfg.cov_nme = {'Diagnosis'};            
            fcfg.plt = 1;           
            fcfg.out_dir = qal_dir;            
            com_bat_epd_dta = ejk_ComBat(fcfg);      
            
            % Save
            cell2csv( [roi_fle_nme(1:end-4) '_' 'ComBat' '.csv'], [ [ 'sbj_nme' ; neu_dta_sbj ] msc_out [ neu_dta_col ; num2cell(com_bat_epd_dta) ] ])
            
        end
    end
end

%% QC post-Harmonized Neurobio
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            % Filenames
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '/' ];
            else
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM}  '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '/' ];
            end
            
            % Load Data
            if ~strcmpi(inp_icv{iF}{iM},'')
                roi_lat_fle_nme = [roi_fle_nme(1:end-4) '_LateralityIndex.csv'];
                roi_fle_nme     = [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '_ComBat.csv'];                
                
                fcfg = [];
                fcfg.dta_loc = roi_fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col, msc_out] = ejk_dta_frm( fcfg );
                
                
                fcfg = [];
                fcfg.dta_loc = roi_lat_fle_nme;
                fcfg.dta_col = 5;
                [ neu_lat_dta, neu_lat_dta_sbj, neu_lat_dta_col, msc_out] = ejk_dta_frm( fcfg );
                
                
            else
                roi_lat_fle_nme = [roi_fle_nme(1:end-4) '_LateralityIndex.csv'];
                roi_fle_nme     = [roi_fle_nme(1:end-4) '_ComBat.csv']; 
                
                fcfg = [];
                fcfg.dta_loc = roi_fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col, msc_out] = ejk_dta_frm( fcfg );
                
                fcfg = [];
                fcfg.dta_loc = roi_lat_fle_nme;
                fcfg.dta_col = 5;
                [ neu_lat_dta, neu_lat_dta_sbj, neu_lat_dta_col, msc_out] = ejk_dta_frm( fcfg );
            end
                        
            % Bilateral QC
            fcfg = [];
            fcfg.sbj_nme = neu_dta_sbj;
            fcfg.dta     = cell2mat(neu_dta);
            fcfg.dta_lbl = neu_dta_col;
            fcfg.out_dir     = qal_dir(1:end-1);
            fcfg.out_pre_fix = [ inp_pre{iF} '_' 'ComBat'];
            ejk_qc_roi(fcfg)
            
            % Laterality QC
            fcfg = [];
            fcfg.sbj_nme = neu_lat_dta_sbj;
            fcfg.dta     = cell2mat(neu_lat_dta);
            fcfg.dta_lbl = neu_lat_dta_col;
            fcfg.out_dir     = [ qal_dir(1:end-1) '_' 'LateralityIndex' ];
            fcfg.out_pre_fix = inp_pre{iF};
            ejk_qc_roi(fcfg)
            
        end
    end
end

%% Apply Total Imaging QC to grp
% Load QC exclusions
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'QC_overall_ejk.csv'];
fcfg.dta_col = 2;
[ qal_dta, qal_dta_sbj, qal_dta_col] = ejk_dta_frm( fcfg );

% Apply to groups
img_qal_rmv = find(sum(cell2mat(qal_dta(:,2)),2)>0);

fcfg = [];
fcfg.rmv_ind = img_qal_rmv;
fcfg.grp     = grp;
img_qal_grp  = ejk_qal_grp(fcfg);

grp = img_qal_grp;
save([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'],'grp');

grp = vis_qal_grp;
save([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_vis_qal.mat'],'grp');

%% Apply Total Imaging & Cognitive QC to grp
% Load QC exclusions
tst_use     = find(strcmpi(cog_dta_col,'lm2_pre'));
cog_cut_off = 3;

low_cog_rmv = find(cell2mat(cog_dta(:,tst_use))<=cog_cut_off);

% Apply to groups
fcfg = [];
fcfg.rmv_ind = low_cog_rmv;
fcfg.grp     = img_qal_grp;
cog_grp      = ejk_qal_grp(fcfg);

grp = cog_grp;
save([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_cog_qal.mat'],'grp');










