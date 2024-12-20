

function ejk_bar(cfg)

if ~isfield(cfg,'fce_wdt'); cfg.fce_wdt = repmat(0.4,1,numel(cfg.xdt)); end

if ~isfield(cfg,'sbp')
    axe_plot = figure('Visible','off'); hold on;
end

for iR = 1:numel(cfg.xdt)
    if ~isempty(cfg.ydt{iR})
        
                rectangle('Position',[ cfg.xdt{iR}-(cfg.fce_wdt(iR)/2) 0 cfg.fce_wdt(iR) cfg.ydt{iR} ],'FaceColor',cfg.fce_col{iR},'EdgeColor',[0   0   0],'LineWidth',3);
        rectangle('Position',[ cfg.xdt{iR}-(cfg.fce_wdt(iR)/2) 0 cfg.fce_wdt(iR) cfg.ydt{iR} ],'FaceColor',cfg.fce_col{iR},'EdgeColor',[.70 .70 .70],'LineWidth',1.5);       
        
        if isfield(cfg,'ydt_err')
            err_bar(iR) = errorbar(cfg.xdt{iR},cfg.ydt{iR},nan,cfg.ydt_err{iR},'LineWidth',4);
            err_bar(iR).Color = cfg.fce_col{iR};
            err_bar(iR).CapSize = 20;
        end
    end
    
end

if ~isfield(cfg,'ylm')
    ylim([min(cell2mat(cfg.ydt))-(min(cell2mat(cfg.ydt)) * 0.3)  max(cell2mat(cfg.ydt))+(max(cell2mat(cfg.ydt)) * 0.3) ])
else
    ylim(cfg.ylm)
end
edg_lng = (cfg.fce_wdt(iR)/2) + (cfg.fce_wdt(iR)/4); 
xlim([ min(cell2mat(cfg.xdt))-edg_lng max(cell2mat(cfg.xdt))+edg_lng ]);
set(gca,'XTick',[]);
set(gca,'XTick',cell2mat(cfg.xdt),'XTickLabel',cfg.xlb);
xtickangle(45)

xlm_hld = xlim;
ylm_hld = ylim;
line( [xlm_hld], [ylm_hld(1) ylm_hld(1)], 'LineWidth', 3, 'Color', 'k');

if isfield(cfg,'ttl'); title(mmil_spec_char(cfg.ttl,{'_'})); end
ylabel(mmil_spec_char(cfg.ylb,{'_'}));

axe_hld = gca;
axe_hld.LineWidth = 3;

if isfield(cfg,'out_dir')
    print([cfg.out_dir '/' cfg.out_nme '.png'],'-dpng')
    close all
end

end