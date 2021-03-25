

function ejk_scatter(cfg)

if ~isfield(cfg,'sbp')
    if ~isfield(cfg,'out_dir')
        figure();
    else
        figure('Visible','off');
    end
end

if all( cellfun( @numel , cfg.xdt ) == 1 )
    org_xdt = cfg.xdt;
    for iX = 1:numel(cfg.xdt)
        min_hld = cfg.xdt{iX} - 0.15;
        max_hld = cfg.xdt{iX} + 0.15;
        cfg.xdt{iX} = (max_hld-min_hld) .* rand(numel(cfg.ydt{iX}),1) + min_hld;
    end
end

for iR = 1:numel(cfg.xdt)
    if ~isempty(cfg.ydt{iR})
        scatter(cfg.xdt{iR},cfg.ydt{iR},48,'o','MarkerEdgeColor',cfg.edg_col{iR},'MarkerFaceColor',cfg.fce_col{iR},'Linewidth',0.4); hold on;
    end    
end

if isfield(cfg,'trd_lne')
    for iL = 1:numel(cfg.trd_lne)
        if cfg.trd_lne == 1
            lne_obj = lsline();
            lne_obj.LineWidth = 2;
            lne_obj.Color = [0 0 0];
        end
    end
end

if ~all(cat(1,cfg.xdt{:})==0); xlim([ roundsd(min(cat(1,cfg.xdt{:}))-abs(min(cat(1,cfg.xdt{:}))*0.1),2) roundsd(max(cat(1,cfg.xdt{:}))+(max(cat(1,cfg.xdt{:}))*0.1),2)]); end

ylm_hld = [ roundsd(min(cat(1,cfg.ydt{:}))-abs(min(cat(1,cfg.ydt{:}))*0.1),2) roundsd(max(cat(1,cfg.ydt{:}))+(max(cat(1,cfg.ydt{:}))*0.1),2)];
if ~isfield(cfg,'ylm') && ~(all(ylm_hld==0)) && ~(all(isnan(ylm_hld)))
    if ~all(cat(1,cfg.ydt{:})==0); ylim(ylm_hld); end
elseif isfield(cfg,'ylm')
    ylim(cfg.ylm)
end

if isfield(cfg,'xlb') && numel(cfg.xlb)==1
    xlabel(cfg.xlb)
elseif isfield(cfg,'xlb') && numel(cfg.xlb)>1
    set(gca,'XTick',unique([org_xdt{:}]));
    set(gca,'XTickLabel',cfg.xlb);
    xtickangle(45)
end

ylabel(mmil_spec_char(cfg.ylb{1},{'_'},{' '}))

if isfield(cfg,'ttl')
    title(mmil_spec_char(cfg.ttl,{'_'},' '))
end

if isfield(cfg,'out_dir')
    ejk_chk_dir(cfg.out_dir)
    tightfig();
    print([cfg.out_dir '/' cfg.out_nme '.png'],'-dpng')
    close all
end

end