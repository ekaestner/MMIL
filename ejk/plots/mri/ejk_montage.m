
function ejk_montage(cfg)

dim_hld = ceil(sqrt(numel(cfg.slc)));

figure('Visible','off')
for iSL = 1:numel(cfg.slc)
    
    sub_plt_hld = subplot(dim_hld,dim_hld,iSL,'Visible','off');
    
    %
    avg_brn_hld = rot90(squeeze(cfg.brn_dta(:,cfg.slc(iSL),:)));
    
    dta_hld = rot90(squeeze(cfg.plt_dta(:,cfg.slc(iSL),:)));
    pvl_hld = rot90(squeeze(cfg.pvl_dta(:,cfg.slc(iSL),:)));
      
    %
    dta_hld(pvl_hld>cfg.pvl_msk) = nan;
    dta_hld(avg_brn_hld<cfg.brn_cut_off) = nan;

    avg_brn_hld(avg_brn_hld<cfg.brn_cut_off) = nan;

    %
    imAlpha_msk = ones(size(avg_brn_hld));
    imAlpha_msk(isnan(avg_brn_hld))=0;
        
    imAlpha = ones(size(dta_hld));
    imAlpha(isnan(dta_hld))=0;
    
    %
    ax1 = axes('Position',sub_plt_hld.Position);
    img_one = imagesc( ax1, avg_brn_hld, 'AlphaData', imAlpha_msk, [ 0.15 0.75 ] ); hold on;
    axis off
    hold all;
    
    %
    ax2 = axes('Position',sub_plt_hld.Position);
    img_two = imagesc( ax2, dta_hld, 'AlphaData', imAlpha, [ cfg.ylm_low cfg.ylm_hgh ] );
    %     img_two.AlphaData = 0.5;
    axis off;
    %     set(gcf, 'InvertHardcopy', 'off')
    
    %
    linkaxes([ax1,ax2])
    ax1.Visible = 'off';
    ax1.XTick = [];
    ax1.YTick = [];
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    %
    colormap(ax1,'gray')
    colormap(ax2,'cool')
    
end

cb1 = colorbar(ax2,'Position',[ .94 .1 .0175 .815]);

print(gcf,[cfg.plt_dir '/' cfg.plt_nme '.png'],'-dpng');
close all;

end
