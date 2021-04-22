
function ejk_qc_surfaces_roi(cfg)

%%
ejk_chk_dir( cfg.out_dir );

load( cfg.dta_hld ); % srf_dta ; srf_dta_sbj

[ prc_loc, prc_lbl, ~]=fs_read_annotation( [ cfg.prc_loc '/' cfg.hms cfg.prc_nme ] );

%% Make overall correlation
% Make average surf data of roi's
srf_roi_hld = nan( numel(cfg.roi_sbj), numel(cfg.roi_lbl));
for iS = 1:numel(cfg.roi_sbj)
    sbj_ind = find(strcmpi( srf_dta_sbj, cfg.roi_sbj{iS}));
    
    for iRO = 1:numel(cfg.roi_lbl)
        
        prc_ind = find(strcmpi( prc_lbl, cfg.roi_lbl{iRO}));
                
        srf_roi_hld(iS, iRO) = nanmean( srf_dta(sbj_ind, prc_loc==prc_ind) );
        
    end
end

% Correlate surf data and roi's
sbj_cor = nan(numel(cfg.roi_sbj),1);
for iS = 1:numel(cfg.roi_sbj)
    cor_hld       = corrcoef( srf_roi_hld(iS,:), cfg.roi_hld(iS,:) );
    sbj_cor(iS,1) = cor_hld(1,2);
end

% Save 
cell2csv( [ cfg.out_dir '/' cfg.out_pre_fix '_subject_correlation.csv' ] , [cfg.roi_sbj num2cell(sbj_cor)])

% Make Plot
fcfg = [];

fcfg.ydt = {[1] sbj_cor [1]};
fcfg.xdt = {0 1 2};

fcfg.edg_col = {rgb('white') rgb('black') rgb('white')};
fcfg.fce_col = {rgb('white') rgb('black') rgb('white')};

fcfg.ylm = [0 1];

fcfg.ylb = {'Surface/ROI correlation (r)'};

fcfg.out_dir = cfg.out_dir;
fcfg.out_nme = cfg.out_pre_fix;

ejk_scatter(fcfg)

%% Make ROIs for important ROIs
use_sbj = find(~isnan(sbj_cor));

roi_dff = nan( numel(sbj_cor), 1);
chg_rvl = nan( numel(sbj_cor), 1);
for iRO = 1:numel( cfg.prc_lbl )
    
    % Get Correlation 
    roi_ind = find(strcmpi( cfg.roi_lbl, cfg.prc_lbl{iRO}));   
    
    cor_hld = corrcoef( srf_roi_hld( use_sbj, roi_ind), cfg.roi_hld( use_sbj,roi_ind) );
        cor_hld = cor_hld(1,2);
    
    % Get Leave-out Correlation
    for iS = 1:numel(sbj_cor)
        if ~isnan(sbj_cor(iS))
            
            use_sbj_tmp = setxor(iS, use_sbj);
            
            cor_hld_tmp = corrcoef( srf_roi_hld( use_sbj_tmp, roi_ind), cfg.roi_hld( use_sbj_tmp,roi_ind) );
                cor_hld_tmp = cor_hld_tmp(1,2);
        
            chg_rvl(iS,1) = abs(cor_hld_tmp - cor_hld);
            
            roi_dff(iS,1) = abs( srf_roi_hld( iS, roi_ind) - cfg.roi_hld( iS, roi_ind) );
        end
    end
    
    % Make Plot
    fcfg = [];
    
    fcfg.ydt = { srf_roi_hld( use_sbj_tmp, roi_ind) };
    fcfg.xdt = { cfg.roi_hld( use_sbj_tmp,roi_ind) };
    
    fcfg.edg_col = { rgb('black') };
    fcfg.fce_col = { rgb('black') };
        
    fcfg.ylb = {'Surface'};
    fcfg.xlb = {'ROI'};
    
    fcfg.out_dir = cfg.out_dir;
    fcfg.out_nme = [cfg.out_pre_fix '_' cfg.prc_lbl{iRO}];
    
    ejk_scatter(fcfg)

    % Output 
    cell2csv( [ cfg.out_dir '/' cfg.out_pre_fix '_' cfg.prc_lbl{iRO} '_correlation.csv' ], [ cfg.roi_sbj num2cell([  srf_roi_hld( :, roi_ind) cfg.roi_hld( :,roi_ind) roi_dff chg_rvl ])])
    
end


end













