% 
% 
% 
% 
% 

function ejk_check_surfaces_roi(cfg)

ejk_chk_dir( cfg.out_dir )

%%
ejk_chk_dir( cfg.out_dir );

load( cfg.dta_hld ); % srf_dta ; srf_dta_sbj

[ prc_loc, prc_lbl, ~]=fs_read_annotation( [ cfg.prc_loc '/' cfg.hms cfg.prc_nme ] );

srf_hld = fs_read_surf( [ cfg.prc_loc '/' 'fsaverage' '/' 'surf' '/' cfg.hms '.' 'pial' ] );

rvl_srf = load( cfg.rvl_fle );

sbj_use = find(strcmpi( cfg.grp, cfg.grp_cmp));

%%
% Make average surf data of roi's
srf_roi_hld = nan( numel(sbj_use), numel(cfg.prc_lbl));
for iS = 1:numel(sbj_use)
    sbj_ind = find(strcmpi( srf_dta_sbj, cfg.sbj_nme{sbj_use(iS)}));
    
    sbj_ind_hld(iS,1) = sbj_ind;
    
    for iRO = 1:numel(cfg.prc_lbl)
        
        % ROI comparison
        prc_ind = find(strcmpi( prc_lbl, cfg.prc_lbl{iRO}));
                
        srf_roi_hld(iS, iRO) = nanmean( srf_dta(sbj_ind, prc_loc==prc_ind) );
        
        % Center Point
        pnt_ind = find(prc_loc==prc_ind);
        [ ~, ind_pnt] = min(mean(abs(srf_hld.vertices(prc_loc==prc_ind,:) - mean(srf_hld.vertices(prc_loc==prc_ind,:),1)),2));
        pnt_ind_hld(iRO) = pnt_ind(ind_pnt);
        
        %
        rvl_ovr_hld{iRO} = rvl_srf.rvalues(prc_loc==prc_ind);
        
    end
end

% Correlation Hld
cor_hld = cfg.cor(sbj_use);

%% Make scatter
for iRO = 1:numel(cfg.prc_lbl)
    
    roi_ind = find(strcmpi( cfg.roi_lbl, cfg.prc_lbl{iRO}));
    
    figure('Visible','off');
    
    % Original Correlation
    rvl_hld{1} = corrcoef(cor_hld, cfg.roi_hld( sbj_use, roi_ind), 'rows', 'pairwise' );
    subplot(3,2,1)
    scatter( cor_hld, cfg.roi_hld( sbj_use, roi_ind), 'r', 'MarkerFaceColor', 'r' );
    title(['ROI/Covariate r = ' num2str(rvl_hld{1}(1,2))])
    
    % Surface Correlation
    rvl_hld{2} = corrcoef(cor_hld, srf_roi_hld(sbj_use,iRO), 'rows', 'pairwise' );
    subplot(3,2,2)
    scatter( cor_hld, srf_roi_hld(sbj_use,iRO), 'm', 'MarkerFaceColor', 'm' );
    title(['Surface/Covariater = ' num2str(rvl_hld{2}(1,2))])
    
    % Original/Surface Correlation
    rvl_hld{3} = corrcoef(srf_roi_hld(sbj_use,iRO), cfg.roi_hld( sbj_use, roi_ind), 'rows', 'pairwise' );
    subplot(3,2,3)
    scatter( srf_roi_hld(sbj_use,iRO), cfg.roi_hld( sbj_use, roi_ind), 'k', 'MarkerFaceColor', 'k');
    title(['ROI/Surface r = ' num2str(rvl_hld{3}(1,2))])
   
    % Original/Surface Correlation
    rvl_hld{4} = corrcoef( cor_hld, srf_dta(sbj_ind_hld, pnt_ind_hld(iRO)), 'rows', 'pairwise' );
    subplot(3,2,4)
    scatter( cor_hld, srf_dta(sbj_ind_hld, pnt_ind_hld(iRO)), 'g', 'MarkerFaceColor', 'g');
    title(['ROI/Point r = ' num2str(rvl_hld{4}(1,2))])
    
    % Original/Surface Correlation
    subplot(3,1,3)
    hist( rvl_ovr_hld{iRO}, 1000 );
    xlim([-1 1]); ylm_hld = get(gca,'ylim');
    line( [ rvl_hld{1}(1,2) rvl_hld{1}(1,2)  ], ylm_hld, 'Color', 'r', 'LineWidth',2);
    line( [ rvl_hld{2}(1,2) rvl_hld{2}(1,2)  ], ylm_hld, 'Color', 'm', 'LineWidth',2);
    line( [ rvl_hld{4}(1,2) rvl_hld{4}(1,2)  ], ylm_hld, 'Color', 'g', 'LineWidth',2);
    
    tightfig();
    print( [ cfg.out_dir '/' cfg.hms '_' cfg.prc_lbl{iRO} '.png'], '-dpng')
    close all
    
end



end

% fcfg = []; %%% make plot
% 
% fcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/';
% fcfg.out_pre_fix = 'lhs_STG_midpoint';
% 
% fcfg.vtx_cor = { [ 160962          ]    [160962] }; % 
% fcfg.vtx_col = { { rgb('bright red') } {rgb('bright red')} };
% 
% fcfg.hme_wrk = 1;
% 
% ejk_highlight_vertex(fcfg)













