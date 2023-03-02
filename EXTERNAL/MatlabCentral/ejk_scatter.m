

function ejk_scatter(cfg)

if ~isfield(cfg,'mkr_typ'); cfg.mkr_typ = repmat({'o'},1,numel(cfg.xdt)); end 
if ~isfield(cfg,'mkr_sze'); cfg.mkr_sze = repmat(48,1,numel(cfg.xdt)); end 
if ~isfield(cfg,'jtr'); cfg.jtr = 1; end 


if ~isfield(cfg,'sbp')
    if ~isfield(cfg,'out_dir')
        figure();
    else
        figure('Visible','off');
    end
end
hold on;

if isfield(cfg,'box_plt') && all( cellfun( @numel , cfg.xdt(cfg.box_plt) ) == 1 )
    for iL = 1:numel(cfg.box_plt)
        if cfg.box_plt(iL) == 1
            boxchart( repmat(cfg.xdt{iL},1,numel(cfg.ydt{iL})), cfg.ydt{iL}', 'BoxFaceColor', cfg.box_plt_col{iL}, 'BoxFaceAlpha', 0.1, ...
                      'WhiskerLineColor', cfg.box_plt_col{iL}, 'LineWidth', 2, 'MarkerStyle', 'none' ) 
        end
    end
end

if all( cellfun( @numel , cfg.xdt ) == 1 ) && cfg.jtr
    org_xdt = cfg.xdt;
    for iX = 1:numel(cfg.xdt)
        min_hld = cfg.xdt{iX} - 0.15;
        max_hld = cfg.xdt{iX} + 0.15;
        cfg.xdt{iX} = (max_hld-min_hld) .* rand(numel(cfg.ydt{iX}),1) + min_hld;
    end
end

for iR = 1:numel(cfg.xdt)
    if ~isempty(cfg.ydt{iR})
        
        dta_plt = scatter(cfg.xdt{iR},cfg.ydt{iR},cfg.mkr_sze(iR),cfg.mkr_typ{iR},'MarkerEdgeColor',cfg.edg_col{iR},'MarkerFaceColor',cfg.fce_col{iR},'Linewidth',0.4); hold on;
        if isfield(cfg,'aph_val'); dta_plt.MarkerFaceAlpha = cfg.aph_val; dta_plt.MarkerEdgeAlpha = cfg.aph_val; end
        if isfield(cfg,'trd_lne')
            if cfg.trd_lne(iR) == 1
                use_idx = ~isnan(cfg.xdt{iR}) & ~isnan(cfg.ydt{iR});
                trd_lne_fit = polyfit(cfg.xdt{iR}(use_idx),cfg.ydt{iR}(use_idx),1);
                trd_lne_fit = polyval(trd_lne_fit,cfg.xdt{iR}(use_idx));
                plot( cfg.xdt{iR}(use_idx), trd_lne_fit, 'Color', cfg.fce_col{iR},'LineWidth',2)
            end
        end
    end
end

xdt_tot = cat(1,cfg.xdt{:}); ydt_tot = cat(1,cfg.ydt{:});
rmv_ind = isnan(xdt_tot) | isnan(ydt_tot);
xdt_tot(rmv_ind) = []; ydt_tot(rmv_ind) = []; 

if isfield(cfg,'xlm')
   xlim(cfg.xlm);   
elseif ~all(xdt_tot==0); 
    xlim([ roundsd(min(xdt_tot)-abs(min(xdt_tot)*0.1),2) roundsd(max(xdt_tot)+(max(xdt_tot)*0.1),2)]);
end

ylm_hld = [ roundsd(min(ydt_tot)-abs(min(ydt_tot)*0.1),2) roundsd(max(ydt_tot)+(max(ydt_tot)*0.1),2)];
if ~isfield(cfg,'ylm') && ~(all(ylm_hld==0)) && ~(all(isnan(ylm_hld)))
    if ~all(ydt_tot==0); ylim(ylm_hld); end
elseif isfield(cfg,'ylm')
    ylim(cfg.ylm)
end

if isfield(cfg,'xlb') && numel(cfg.xlb)==1 
    xlabel(mmil_spec_char(cfg.xlb{1},{'_'},{' '}))
elseif isfield(cfg,'xlb') && numel(cfg.xlb)>1 && all( cellfun( @numel , org_xdt ) == 1 )
    set(gca,'XTick',unique([org_xdt{:}]));
    set(gca,'XTickLabel',cellfun(@(x) mmil_spec_char(x,{'_'},{' '}),cfg.xlb,'UniformOutput',false));
    xtickangle(45)
elseif all(cellfun(@numel,cfg.xdt)==numel(cfg.xlb))
    set(gca,'XTick',unique(round([cfg.xdt{:}])));
    set(gca,'XTickLabel',cfg.xlb);
    xtickangle(45)
end

if isfield(cfg,'ylb'); ylabel(mmil_spec_char(cfg.ylb{1},{'_'},{' '})); end

if isfield(cfg,'ttl'); title(mmil_spec_char(cfg.ttl,{'_'},{' '})); end

if isfield(cfg,'hln')
    line( get(gca,'xlim') ,[cfg.hln cfg.hln], 'Color', cfg.hln_col);
    lne_ord = get(gca, 'Children');
    set(gca,'Children',lne_ord([2:end 1]))
end

if isfield(cfg,'vln')
    line( [cfg.vln cfg.vln], get(gca,'ylim'), 'Color', cfg.vln_col);
    lne_ord = get(gca, 'Children');
    set(gca,'Children',lne_ord([2:end 1]))
end

if isfield(cfg,'out_dir')
    ejk_chk_dir(cfg.out_dir)
    tightfig();
    print([cfg.out_dir '/' cfg.out_nme '.png'],'-dpng')
    close all
end

end