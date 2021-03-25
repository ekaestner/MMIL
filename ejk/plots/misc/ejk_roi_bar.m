
function ejk_roi_bar(cfg)

for iDT = 1:numel(cfg.dta_lbl)
    
    if ~exist([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/']);                             mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/']); end
    if ~exist([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'BAR' '/' cfg.dta_lbl{iDT} ]); mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'BAR' '/' cfg.dta_lbl{iDT} ]); end
    
    fst_nme = fieldnames(cfg.ydt{iDT}); fst_nme(strcmpi(fst_nme,'sbj_nme')) = [];
    
    int_num_hld = [];
    int_nme_hld = cell(0);
    
    for iFN = 1:numel(fst_nme)
        
        %% Gather Data
        for iST = 1:numel(cfg.grp_col)
            
            if ~isnan(cfg.grp_col(iST))
                for iXY = 1:numel(cfg.grp_nme{iST})
                    
                    ydt{iST}{iXY} = cfg.ydt{iDT}.(fst_nme{iFN})(strcmpi(cfg.grp(2:end,cfg.grp_col(iST)),cfg.grp_nme{iST}{iXY}),1);
                    rmv_ind = unique( find(isnan(ydt{iST}{iXY})) );
                    ydt{iST}{iXY}(rmv_ind) = [];
                    
                    num_sbj{iST}(iXY) = numel(ydt{iST}{iXY});
                    
                end
            end
            
        end
        
        
        %% Make Stats
        for iST = 1:numel(cfg.grp_col)
             if isfield(cfg,'plt_anv') && ~isempty(cfg.plt_anv{iST}) && ~all(isnan(cat(1,ydt{iST}{:})))

                 scfg = [];
                 
                 scfg.dta     = ydt{iST}(cfg.plt_anv{iST});
                 scfg.grp_nme = cfg.grp_nme{iST}(cfg.plt_anv{iST});
                 
                 [ ~ , stt{iST}{iFN,1} , fll{iST}(iFN,:) ] = mmil_roi_anv(scfg); %
                 
             end
        end
        
        %% Make Plot
        plt_dim = ceil(sqrt(numel(cfg.grp_col)));
        
        figure('Visible','off','Position',[0 0 1080 1080])
        for iPL = 1:numel(cfg.grp_col)
            
            sbp(iPL) = subplot(plt_dim,plt_dim,iPL); hold on;
            
            if ~isnan(cfg.grp_col(iPL))
                
                if ~all(isnan(cat(1,ydt{iPL}{:})))
                    
                    fcfg = [];
                    
                    fcfg.sbp = sbp(iPL);
                    
                    fcfg.ydt     = num2cell(cellfun(@mean,ydt{iPL}));
                    fcfg.ydt_err = num2cell(cellfun(@std,ydt{iPL})); % ./ sqrt(num_sbj{iPL}
                    fcfg.xdt     = num2cell(cfg.xdt{iPL});
                    
                    fcfg.fce_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                    fcfg.edg_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                    
                    fcfg.ttl = cfg.nme_col{iPL};
                    fcfg.xlb = cfg.grp_nme{iPL};
                    fcfg.ylb = fst_nme{iFN};
                    
                    if ~all(cell2mat(fcfg.ydt_err)==0)
                        ejk_bar(fcfg)
                    end
                    
                end
            end
            
        end
        
        %% Plot Stats
        for iPL = 1:numel(cfg.grp_col)
            if  ~isnan(cfg.grp_col(iPL)) && ~all(isnan(cat(1,ydt{iPL}{:}))) && ~all(cell2mat(num2cell(cellfun(@std,ydt{iPL})))==0)
                            
                scfg = [];
                
                scfg.sub_plt = sbp(iPL);
                scfg.scr_top = max(cellfun( @mean , ydt{iPL} )) + max(cellfun( @std , ydt{iPL} ));
                
                scfg.dta     = ydt{iPL};
                scfg.xdt     = cfg.xdt{iPL};
                
                scfg.plt_stt_cmp = cfg.plt_cmp{iPL};
                scfg.tst_typ = 'ttest';
                scfg.sig_val = .05;
                scfg.sig_oly = 1;
                
                mmil_bar_stt(scfg);
            
            end
        end

        %% Save Figure
        set(gcf,'Position',[0 0 1080 1080]);
        tightfig();
        print([ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'BAR' '/' cfg.dta_lbl{iDT} '/' fst_nme{iFN} '_' cfg.dta_lbl{iDT} '.png'],'-dpng')
        close all
        
        clear ydt num_sbj sbp
        
    end
    
    %% Save Stats
    for iST = 1:numel(cfg.grp_col)
        if isfield(cfg,'plt_anv') && ~isempty(cfg.plt_anv{iST})
            
            col_hed = nchoosek(1:numel(cfg.grp_nme{iST}),2);
            for iCH = 1:size(col_hed,1)
               col_hed_nme{iCH} = [ cfg.grp_nme{iST}{col_hed(iCH,1)} '_VS_' cfg.grp_nme{iST}{col_hed(iCH,2)} ]; 
            end
            
            cell2csv([ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'BAR' '/' cfg.dta_lbl{iDT} '/' 'z' cfg.dta_lbl{iDT} '_' cfg.nme_col{iST} '_' 'PValues' '.csv'],[ { '' 'ANOVA' col_hed_nme{:} } ; fst_nme num2cell(fll{iST}) ]);
            cell2csv([ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'BAR' '/' cfg.dta_lbl{iDT} '/' 'z' cfg.dta_lbl{iDT} '_' cfg.nme_col{iST} '_' 'ANOVA' '.csv'],[ { '' 'ANOVA' } ; fst_nme stt{iST} ]);
            
            clear col_hed_nme
            
        end
    end
    
    clear pnv stt fll
    
    
end

end