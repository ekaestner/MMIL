%

function ejk_roi_scatter(cfg)

if ~isfield(cfg,'pvl_cut'); cfg.pvl_cut = 0.05; end
if strcmpi(cfg.grp{1,1},'SubjID'); cfg.grp = cfg.grp(2:end,:); end

for iDT = 1:numel(cfg.dta_lbl)
    
    fst_nme = fieldnames(cfg.ydt{iDT}); fst_nme(strcmpi(fst_nme,'sbj_nme')) = [];
    if isstruct(cfg.xdt{iDT}); scd_nme = fieldnames(cfg.xdt{iDT}); scd_nme(strcmpi(scd_nme,'sbj_nme')) = []; end
    
    int_num_hld = [];
    int_nme_hld = cell(0);
    
    for iFN = 1:numel(fst_nme)
        
        if isstruct(cfg.xdt{iDT})
            if isfield(cfg,'out_dir')
                out_dir_hld = [ cfg.out_dir '/' cfg.dta_lbl{iDT} '/' fst_nme{iFN} ];
            else
                out_dir_hld = [ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'SCATTER' '/' cfg.dta_lbl{iDT} '/' fst_nme{iFN} ];
            end
        else
            if isfield(cfg,'out_dir')
                out_dir_hld = [ cfg.out_dir '/' cfg.dta_lbl{iDT} '/' ];
            else
                out_dir_hld = [ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'SCATTER' '/' cfg.dta_lbl{iDT} '/' ];
            end
        end
        
        ejk_chk_dir(out_dir_hld)
        
        %% Stats
        if isstruct(cfg.xdt{iDT})
            
            for iSN = 1:numel(scd_nme)
                
                out_csv_ind = 1;
                for iST = 1:numel(cfg.grp_col)
                    
                    for iXY = 1:numel(cfg.grp_nme{iST})
                        
                        xdt{iXY} = cfg.xdt{iDT}.(scd_nme{iSN})(strcmpi(cfg.grp(:,cfg.grp_col(iST)),cfg.grp_nme{iST}{iXY}),1);
                        ydt{iXY} = cfg.ydt{iDT}.(fst_nme{iFN})(strcmpi(cfg.grp(:,cfg.grp_col(iST)),cfg.grp_nme{iST}{iXY}),1);
                        
                        rmv_ind = unique( [ find(isnan(xdt{iXY})) ; find(isnan(ydt{iXY})) ] );
                        
                        xdt{iXY}(rmv_ind) = [];
                        ydt{iXY}(rmv_ind) = [];
                        
                        num_sbj = numel(xdt{iXY});
                        [ cor_hld , pvl_hld ] = corrcoef( xdt{iXY} , ydt{iXY} );
                        
                        out_cor_csv(iSN,out_csv_ind) = roundsd(cor_hld(1,2),2);
                        out_pvl_csv(iSN,out_csv_ind) = roundsd(pvl_hld(1,2),2);
                        
                        out_csv_ind = out_csv_ind + 1;
                        
                    end
                    
                    clear xdt ydt
                    
                end
                
                if any(out_pvl_csv(iSN,:)<cfg.pvl_cut)
                    int_num_hld(end+1,:) = out_pvl_csv(iSN,:);
                    int_nme_hld(end+1,:) = [ fst_nme(iFN)  scd_nme(iSN) ];
                end
                
            end
            
            cell2csv( [ out_dir_hld '/' 'corr' '_' fst_nme{iFN} '.csv']    , [ [ {''} cat(2,cfg.grp_nme{:})] ; [scd_nme num2cell(out_cor_csv)] ] )
            cell2csv( [ out_dir_hld '/' 'pval' '_' fst_nme{iFN} '.csv']    , [ [ {''} cat(2,cfg.grp_nme{:})] ; [scd_nme num2cell(out_pvl_csv)] ] )
            
            clear out_cor_csv out_pvl_csv
            
        end
        
        %% Scatter Plot
        if isstruct(cfg.xdt{iDT})
            
            for iSN = 1:numel(scd_nme)
                
                plt_dim = ceil(sqrt(numel(cfg.grp_col)));
                
                figure('Visible','off','Position',[0 0 1080 1080])
                for iPL = 1:numel(cfg.grp_col)
                    
                    sbp = subplot(plt_dim,plt_dim,iPL);
                    
                    for iXY = 1:numel(cfg.grp_nme{iPL})
                        xdt{iXY} = cfg.xdt{iDT}.(scd_nme{iSN})(strcmpi(cfg.grp(:,cfg.grp_col(iPL)),cfg.grp_nme{iPL}{iXY}),1);
                        ydt{iXY} = cfg.ydt{iDT}.(fst_nme{iFN})(strcmpi(cfg.grp(:,cfg.grp_col(iPL)),cfg.grp_nme{iPL}{iXY}),1);
                    end
                    
                    if ~all(isnan(cat(1,xdt{:}))) && ~all(isnan(cat(1,ydt{:})))
                        fcfg = [];
                        fcfg.sbp = sbp;
                        fcfg.xdt = xdt;
                        fcfg.ydt = ydt;
                        fcfg.fce_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                        fcfg.edg_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                        fcfg.xlb = scd_nme{iSN};
                        fcfg.ylb = fst_nme{iFN};
                        ejk_scatter(fcfg)
                        
                        clear xdt ydt
                    end
                    
                end
                
                print([ out_dir_hld '/' scd_nme{iSN} '__' 'BY' '__' fst_nme{iFN} '.png'],'-dpng')
                close all
                
            end
            
        else
            
            plt_dim = ceil(sqrt(numel(cfg.grp_col)));
            
            figure('Visible','off','Position',[0 0 1080 1080])
            for iPL = 1:numel(cfg.grp_col)
                
                sbp = subplot(plt_dim,plt_dim,iPL);
                if ~isnan(cfg.grp_col(iPL))
                    for iXY = 1:numel(cfg.grp_nme{iPL})
                        xdt{iXY} = cfg.xdt{iPL}(iXY);
                        ydt{iXY} = cfg.ydt{iDT}.(fst_nme{iFN})(strcmpi(cfg.grp(:,cfg.grp_col(iPL)),cfg.grp_nme{iPL}{iXY}),1);
                    end
                    
                    if ~all(isnan(cat(1,xdt{:}))) && ~all(isnan(cat(1,ydt{:})))
                        fcfg = [];
                        fcfg.sbp = sbp;
                        fcfg.xdt = xdt;
                        fcfg.ydt = ydt;
                        fcfg.fce_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                        fcfg.edg_col = cellfun(@rgb,cfg.grp_clr{iPL},'uni',0);
                        fcfg.xlb = '';
                        fcfg.ylb = fst_nme{iFN};
                        ejk_scatter(fcfg)
                        
                        clear xdt ydt
                    end
                    
                end
                
            end
            
            print([ out_dir_hld '/' fst_nme{iFN} '.png'],'-dpng')
            close all
            
        end
        
    end
end

if isfield(cfg,'out_dir')
    out_dir_hld = [ cfg.out_dir '/' cfg.dta_lbl{iDT} '/' ];
else
    out_dir_hld = [ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'ROI' '/' 'SCATTER' '/' cfg.dta_lbl{iDT} '/' fst_nme{iFN} ];
end
cell2csv([ out_dir_hld '/' 'overall_pvl.csv' ] , [ [ {''} {''}  cat(2,cfg.grp_nme{:})] ; int_nme_hld num2cell(int_num_hld) ]);

end




