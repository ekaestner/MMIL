

function ejk_bar_vert(cfg)

if ~isfield(cfg,'sbp')
    axe_plot = figure('Visible','off');
end

hold on;

for iR = 1:numel(cfg.xdt)
    if ~isempty(cfg.ydt{iR})
        
        rectangle('Position',[ 0 cfg.ydt{iR}-0.2 cfg.xdt{iR} 0.4  ],'FaceColor',cfg.fce_col{iR},'EdgeColor',rgb('black'));
        
        if isfield(cfg,'xdt_err')
            err_bar(iR) = errorbar(cfg.xdt{iR},cfg.ydt{iR},nan,cfg.xdt_err{iR},'horizontal','LineWidth',4);
            err_bar(iR).Color = cfg.fce_col{iR};
            err_bar(iR).CapSize = 20;
        end
    end
    
end

if ~isfield(cfg,'xlm')
    xlim([min(cell2mat(cfg.xdt))-(min(cell2mat(cfg.xdt)) * 0.3)  max(cell2mat(cfg.ydt))+(max(cell2mat(cfg.xdt)) * 0.3) ])
else
    xlim(cfg.xlm)
end
ylim([ min(cell2mat(cfg.ydt))-0.5 max(cell2mat(cfg.ydt))+0.5 ]);
set(gca,'YTick',[]);
set(gca,'YTick',cell2mat(cfg.ydt),'YTickLabel', cellfun(@(x) mmil_spec_char(x,{'_'}),cfg.ylb,'uni',0) );
ytickangle(45)

if isfield(cfg,'ttl'); title(mmil_spec_char(cfg.ttl,{'_'})); end
xlabel(mmil_spec_char(cfg.xlb,{'_'}));

if isfield(cfg,'out_dir')
    print([cfg.out_dir '/' cfg.out_nme '.png'],'-dpng')
    close all
end

end