
qal_dti_nme     = { 'fiber_FA' 'gwcsurf_FA_wm_ctx' };
qal_dti_roi_nme = { { 'fiber_FA_L_CgH' 'fiber_FA_R_CgH' 'fiber_FA_L_tSLF' 'fiber_FA_R_tSLF' 'fiber_FA_L_IFO' 'fiber_FA_R_IFO' 'fiber_FA_L_ILF' 'fiber_FA_R_ILF' 'fiber_FA_L_Unc' 'fiber_FA_R_Unc' 'fiber_FA_L_CST' 'fiber_FA_R_CST' } ...
                    { 'gwcsurf_FA_wm_ctx_lh_entorhinal' 'gwcsurf_FA_wm_ctx_rh_entorhinal' 'gwcsurf_FA_wm_ctx_lh_fusiform' 'gwcsurf_FA_wm_ctx_rh_fusiform' 'gwcsurf_FA_wm_ctx_lh_inferiortemporal' 'gwcsurf_FA_wm_ctx_rh_inferiortemporal' 'gwcsurf_FA_wm_ctx_lh_lateralorbitofrontal' 'gwcsurf_FA_wm_ctx_rh_lateralorbitofrontal' 'gwcsurf_FA_wm_ctx_lh_precentral' 'gwcsurf_FA_wm_ctx_rh_precentral' 'gwcsurf_FA_wm_ctx_lh_cuneus' 'gwcsurf_FA_wm_ctx_rh_cuneus' } };

%% Setup
[ ~, qal_dti_ind ] = intersect( dti_dev_mse, qal_dti_nme); 

neu_bio_qal_dta = cell( numel(all_sbj),numel(cat(2,qal_dti_roi_nme{:})));
neu_bio_qal_sbj = all_sbj;
neu_bio_qal_col = cell(1,numel(cat(2,qal_dti_roi_nme{:})));

%% Load DTI data & subselect gwcsurf & fibers
col_ind = 1;
for iQ = 1:numel(qal_dti_ind)
    
    if dti_dev_roi{qal_dti_ind(iQ)}(1)==0
        roi_nme_hld = [];
        roi_hms_hld = [];
    else
        roi_nme_hld = { [ roi_loc '/' 'lh.' roi_nme{dti_dev_roi{qal_dti_ind(iQ)}(1)}] [ roi_loc '/' 'rh.' roi_nme{dti_dev_roi{qal_dti_ind(iQ)}(1)}] };
        roi_hms_hld = { 'lh' 'rh' };
    end
    
    fcfg = [];
    fcfg.sbj_nme = all_sbj;
    fcfg.fle_nme = dti_dev_fle;
    fcfg.mes_nme = dti_dev_mse{qal_dti_ind(iQ)};
    fcfg.rcn_nme = rcn_fle;
    fcfg.roi_nme = roi_nme_hld;
    fcfg.roi_hms = roi_hms_hld;
    
    [ neu_bio_dta, neu_bio_wrn ] = ejk_extract_mmps_roi( fcfg );
    
    [ ~, int_ind ] = intersect(neu_bio_dta(1,:),qal_dti_roi_nme{iQ});
    neu_bio_qal_dta(:,col_ind:col_ind+numel(int_ind)-1) = neu_bio_dta(2:end,int_ind);
    neu_bio_qal_col(1,col_ind:col_ind+numel(int_ind)-1) = neu_bio_dta(1,int_ind);
    
    col_ind = col_ind + numel(int_ind);
end

%% ComBat
btc_col = strcmpi(all_col,'site');
sde_col = strcmpi(all_col,'sde_sze_ons');
typ_col = strcmpi(all_col,'sbj_nme');

% comBAT
btc_dta = cell(size(all_sbj,1),1);
btc_dta(strcmpi(all_dta(:,btc_col),'UCSD'))  = {'UCSD'};
btc_dta(strcmpi(all_dta(:,btc_col),'Emory')) = {'Emory'};
btc_dta(strcmpi(all_dta(:,btc_col),'UCSF'))  = {'UCSF'};

cov_dta                                    = cell(size(all_sbj,1),1);
cov_dta(strcmpi(all_dta(:,sde_col),'L'))    = {'LTLE'};
cov_dta(strcmpi(all_dta(:,sde_col),'R'))    = {'RTLE'};
cov_dta(string_find(all_dta(:,typ_col),'fc')) = {'HC'};
cov_dta(cellfun(@isempty,cov_dta))         = {'empty'};

fcfg = [];
fcfg.sbj_nme = all_sbj;
fcfg.dta     = cell2mat(neu_bio_qal_dta);
fcfg.dta_nme = neu_bio_qal_col;
fcfg.btc     = btc_dta;
fcfg.btc_nme = {'Site'};
fcfg.cov     = cov_dta;
fcfg.cov_nme = {'Diagnosis'};
fcfg.plt     = 1;
fcfg.out_dir = qal_dir;
com_bat_neu_bio_qal_dta = ejk_ComBat(fcfg);

%% Run QC
fcfg = [];
fcfg.sbj_nme     = neu_bio_qal_sbj;
fcfg.dta         = com_bat_neu_bio_qal_dta;
fcfg.dta_lbl     = neu_bio_qal_col;
fcfg.out_dir     = [ qal_dir '/' 'DTI' '/' ];
fcfg.out_pre_fix = 'DTI';
ejk_qc_roi(fcfg)






