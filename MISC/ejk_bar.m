

function ejk_bar(cfg)

if ~isfield(cfg,'sbp')
    axe_plot = figure('Visible','off');
end

for iR = 1:numel(cfg.xdt)
    if ~isempty(cfg.ydt{iR})
        
        rectangle('Position',[ cfg.xdt{iR}-0.2 0 0.4 cfg.ydt{iR} ],'FaceColor',cfg.fce_col{iR});
        
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
xlim([ min(cell2mat(cfg.xdt))-0.5 max(cell2mat(cfg.xdt))+0.5 ]);
set(gca,'XTick',[]);
set(gca,'XTick',cell2mat(cfg.xdt),'XTickLabel',cfg.xlb);
xtickangle(45)

if isfield(cfg,'ttl'); title(mmil_spec_char(cfg.ttl,{'_'})); end
ylabel(mmil_spec_char(cfg.ylb,{'_'}));

if isfield(cfg,'out_dir')
    print([cfg.out_dir '/' cfg.out_nme '.png'],'-dpng')
    close all
end

end