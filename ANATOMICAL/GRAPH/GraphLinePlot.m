% cfg = [];
%
% cfg.grp_nme = phe_grp;
% cfg.grp_col = phe_col;
% cfg.leg_loc = 'southwest';
% cfg.bin_cut = bin_cut_off;
%
% cfg.grp_dta = out_dta;
%
% cfg.con_spr_plt = 1;
% cfg.con_nme = pvl_con;
%
% cfg.con_cmp_stt = 1;
% cfg.cmp_nme = pvl_cmp;
% cfg.con_typ = 'con_shf'; 'dff_shf';
%
% cfg.pvl_lvl = pvl_lvl;
%
% cfg.reg_plt = 1;
% cfg.thk_nme = phe_wht_nme(:,15:end);
%
% cfg.ttl = 'Modularity Score';
% cfg.sve_loc = [out_dir '/' 'GRAPH' '/' prj_nme '/' mes_off_int{iB} '/' ];
% cfg.sve_nme = 'p011_Modularity_Score.png';
%
% GraphLinePlot2(cfg)

%%
function GraphLinePlot(cfg)

if ~isfield(cfg,'plt_fig'); cfg.plt_fig = 0; end

%% GROUP PLOT
if cfg.grp_plt
    
    for iC = 1:numel(cfg.grp_ovr)
        
        for iG = 1:numel(cfg.grp_nme{iC})
            cfg.grp_dta.dta{iC}{iG}(cfg.grp_dta.dta{iC}{iG}==Inf) = nan;
        end
        
        figure('Visible','off','Position',[0 0 1080 1080]); hold on;
        
        % Plot Limits
        if ~isfield(cfg,'ylm')
            if size(cfg.grp_dta.dta{iC}{1},1)==1
                ylm_max = max([cfg.grp_dta.dta{iC}{:}])+(max([cfg.grp_dta.dta{iC}{:}]) * .125);
            elseif size([cfg.grp_dta.dta{1}],1)>1
                ylm_max_hld = cat(1,cfg.grp_dta.dta{:});
                ylm_max_hld(ylm_max_hld==inf) = [];
                ylm_max = max(ylm_max_hld(:))+(max(ylm_max_hld(:)) * .125);
            end
            set(gca,'ylim',[0 ylm_max]);
        else
            ylm_max = cfg.ylm(2) + ( cfg.ylm(2) * .125 );
            set(gca,'ylim',[cfg.ylm(1) ylm_max]);
        end
        
        if size(cfg.grp_dta.dta{iC}{1},1)==1
            ylm_cel = max([cfg.grp_dta.dta{iC}{:}])+(max([cfg.grp_dta.dta{iC}{:}]) * [.025:.025:.125]);
        elseif size([cfg.grp_dta.dta{1}],1)>1
            ylm_max_hld = cat(1,cfg.grp_dta.dta{:});
            ylm_max_hld(ylm_max_hld==inf) = [];
            ylm_max = max(ylm_max_hld(:))+(max(ylm_max_hld(:)) *  [.025:.025:.125]);
        end
        
        set(gca,'xlim',[cfg.bin_cut(1)-3 cfg.bin_cut(end)+3]);
        
        % Control Variance Plot
        if isfield(cfg,'con_spr_plt') && cfg.con_spr_plt
            for iBC = 1:numel(cfg.bin_cut)
                pct_plt(1,iBC) = prctile(cfg.grp_dta.con_spr{iC}(:,iBC),100-((cfg.pvl_lvl/2)*100));
                pct_plt(2,iBC) = prctile(cfg.grp_dta.con_spr{iC}(:,iBC),(cfg.pvl_lvl/2)*100);
            end
            patch([cfg.bin_cut fliplr(cfg.bin_cut)],[pct_plt(1,:) fliplr(pct_plt(2,:))],cfg.grp_col{iC}{strcmpi(cfg.grp_nme{iC},cfg.con_nme{iC})},'FaceAlpha',0.1)
        end
        
        % Significance Plot
        if isfield(cfg,'con_typ') && ~isempty(cfg.con_typ)
            
            if strcmpi(cfg.con_typ,'dff_shf')
                
                for iSG = 1:numel(cfg.cmp_nme{iC})
                    
                    cmp_nme = mmil_spec_char(cfg.cmp_nme{iC}{iSG},{' ' '&'});
                    con_loc = strcmpi(cfg.grp_nme,cfg.con_nme);
                    cmp_loc = strcmpi(cfg.grp_nme,cfg.cmp_nme{iC}{iSG});
                    
                    for iBC = 1:numel(cfg.bin_cut)
                        
                        dff_hst = cfg.grp_dta.cmp_stt.(cmp_nme).con_dff_dta{iBC} - cfg.grp_dta.cmp_stt.(cmp_nme).cmp_dff_dta{iBC};
                        dff_act = cfg.grp_dta.dta{con_loc}(iBC) - cfg.grp_dta.dta{cmp_loc}(iBC);
                        
                        pvl_min = prctile(dff_hst,(cfg.pvl_lvl/2)*100);
                        pvl_max = prctile(dff_hst,100-(cfg.pvl_lvl/2)*100);
                        
                        if dff_act<pvl_min || pvl_max<dff_act
                            scatter(cfg.bin_cut(iBC),ylm_cel(iSG),150,cfg.grp_col(cmp_loc,:),'o','filled','MarkerEdgeColor',rgb('light grey'),'LineWidth',2);
                        end
                    end
                end
                
                cfg.sve_nme{iC} = [cfg.sve_nme{iC} '_' 'diff_shuffle'];
                
            elseif  strcmpi(cfg.con_typ,'con_shf')
                
                cfg.grp_dta.con_spr{iC}(cfg.grp_dta.con_spr{iC}==Inf) = nan;
                
                for iSG = 1:numel(cfg.cmp_nme{iC})
                    
                    cmp_loc = strcmpi(cfg.grp_nme{iC},cfg.cmp_nme{iC}{iSG});
                    
                    for iBC = 1:numel(cfg.bin_cut)
                        if cfg.grp_dta.dta{iC}{cmp_loc}(iBC) < pct_plt(2,iBC) || ...
                                pct_plt(1,iBC) < cfg.grp_dta.dta{iC}{cmp_loc}(iBC)
                            scatter(cfg.bin_cut(iBC),ylm_cel(iSG),250,rgb('dark grey'),'o','filled'); hold on;
                            if isfield(cfg,'grp_col_hgh_lgh')
                                scatter(cfg.bin_cut(iBC),ylm_cel(iSG),150,cfg.grp_col{iC}{cmp_loc},'o','filled','MarkerEdgeColor',cfg.grp_col_hgh_lgh(cmp_loc,:),'LineWidth',1.5);
                            else
                                scatter(cfg.bin_cut(iBC),ylm_cel(iSG),150,cfg.grp_col{iC}{cmp_loc},'o','filled');
                            end
                        end
                    end
                end
                
                cfg.sve_nme{iC} = [cfg.sve_nme{iC} '_' 'control_shuffle'];
                
            elseif  strcmpi(cfg.con_typ,'sbj_cmp')
                
                %%% %%% %%% %%% %%% %%%
                
            end
        end
        
        % Plot Data
        for iP = 1:numel(cfg.grp_nme{iC}) % Group Plot
            
            if size(cfg.grp_dta.dta{1},1)==1
                plot(cfg.bin_cut,cfg.grp_dta.dta{iC}{iP},'Color',rgb('dark grey'),'LineWidth',7.75);
                if isfield(cfg,'grp_col_hgh_lgh')
                    plot(cfg.bin_cut,cfg.grp_dta.dta{iC}{iP},'Color',cfg.grp_col_hgh_lgh(plt_loc(iP),:),'LineWidth',6.5);
                end
                leg(iP) = plot(cfg.bin_cut,cfg.grp_dta.dta{iC}{iP},'Color',cfg.grp_col{iC}{iP},'LineWidth',4.25);
            elseif size(cfg.grp_dta.dta{1},1)>1
                dta_hld = mean(cfg.grp_dta.dta{plt_loc(iP)},1);
                
                plot(cfg.bin_cut,dta_hld,'Color',rgb('dark grey'),'LineWidth',7.75);
                if isfield(cfg,'grp_col_hgh_lgh')
                    plot(cfg.bin_cut,dta_hld,'Color',cfg.grp_col_hgh_lgh(plt_loc(iP),:),'LineWidth',6.5);
                end
                leg(iP) = plot(cfg.bin_cut,dta_hld,'Color',cfg.grp_col(plt_loc(iP),:),'LineWidth',4.25);
            end
            
        end
        
        % Scatter if subjects
        if size(cfg.grp_dta.dta{iC}{1},1)>1
            for iP = 1:numel(tot_sbj) % Group Plot
                
                plt_loc(iP) = find(strcmpi(cfg.grp_nme,tot_sbj{iP}));
                
                for iBC = 1:numel(cfg.bin_cut)
                    
                    xvl_hld = (((cfg.bin_cut(iBC)+0.4) - (cfg.bin_cut(iBC)-0.4)) .* rand(1,size(cfg.grp_dta.dta{iP},1))) + (cfg.bin_cut(iBC)-0.4);
                    scatter(xvl_hld,cfg.grp_dta.dta{iP}(:,iBC),75,rgb('black'),'o','filled')
                    scatter(xvl_hld,cfg.grp_dta.dta{iP}(:,iBC),65,cfg.grp_col_hgh_lgh(plt_loc(iP),:),'o','filled')
                    scatter(xvl_hld,cfg.grp_dta.dta{iP}(:,iBC),50,cfg.grp_col(plt_loc(iP),:),'o','filled')
                    
                end
            end
        end
        
        % Legend
        if ~cfg.plt_fig
            lob = legend(leg,cfg.grp_nme{iC},'Location',cfg.leg_loc,'FontSize',20);
            lne_lob = findobj(lob,'type','line');
            set(lne_lob,'linewidth',100);
            legend('boxoff')
        end
        
        % Title
        ylabel(mmil_spec_char(cfg.ttl,{'_'}))
        
        % Save
        set(gcf,'color','white')
        
        if ~cfg.plt_fig
            print([cfg.sve_loc{iC} '/' cfg.sve_nme{iC} '.png'],'-dpng','-r250')
        elseif cfg.plt_fig
            tightfig(); set(gcf,'Position',[0 0 1080 1080]);
            cfg = [];
            cfg.fle_nme = [cfg.sve_loc '/' cfg.sve_nme];
            cfg.prn_typ = 'eps';
            cfg.prn_typ_ext = 'png';
            mmil_print_plot2(cfg)
        end
        
        close all
        
    end
end

%% Region Plots
if isfield(cfg,'reg_lne_plt') && cfg.reg_lne_plt
    for iC = 1:numel(cfg.grp_ovr)
        
        % Individual Line Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot Limits
        if ~isfield(cfg,'ylm')
            dta_hld = [cfg.grp_dta.dta{iC}{:}];
            ylm_max = max(dta_hld(:))+(max(dta_hld(:)) * .125);
        else
            ylm_max = cfg.ylm(2) + ( cfg.ylm(2) * .125 );
            set(gca,'ylim',[cfg.ylm(1) ylm_max]);
        end
        
        ylm_cel = max(dta_hld(:))+(max(dta_hld(:)) * [.025:.025:.125]);
        
        % Split Data
        hms_nme = {'lhs' 'rhs'};
        
        dta_reg_hld{1} = cellfun(@(x) x(string_find(cfg.thk_nme,'lhs'),:),cfg.grp_dta.dta_reg{iC},'uni',0);
        dta_reg_hld{2} = cellfun(@(x) x(string_find(cfg.thk_nme,'rhs'),:),cfg.grp_dta.dta_reg{iC},'uni',0);
        
        if isfield(cfg,'con_spr_plt') && cfg.con_spr_plt
            dta_reg_con_spr_hld{1} = cellfun(@(x) x( : , string_find(cfg.thk_nme,'lhs') , : ),cfg.grp_dta.con_spr_reg,'uni',0);
            dta_reg_con_spr_hld{1} = dta_reg_con_spr_hld{1}{iC};
            dta_reg_con_spr_hld{2} = cellfun(@(x) x( : , string_find(cfg.thk_nme,'rhs') , : ),cfg.grp_dta.con_spr_reg,'uni',0);
            dta_reg_con_spr_hld{2} = dta_reg_con_spr_hld{2}{iC};
        end
        
        thk_nme_hld{1} = cfg.thk_nme(string_find(cfg.thk_nme,'lhs'));
        thk_nme_hld{2} = cfg.thk_nme(string_find(cfg.thk_nme,'rhs'));
        
        if ~isfield(cfg,'reg_inc_plt')
            reg_inc_plt{1} = thk_nme_hld{1};
            reg_inc_plt{2} = thk_nme_hld{2};
        else
            reg_inc_plt{1} = cfg.reg_inc_plt(string_find(cfg.reg_inc_plt,'lhs'));
            reg_inc_plt{2} = cfg.reg_inc_plt(string_find(cfg.reg_inc_plt,'rhs'));
        end
        
        % Plot
        for iH = 1:2
            
            if ~cfg.plt_fig; figure('Visible','off','Position',[0 0 1080 1080]); hold on; end
            
            plt_dim = ceil(sqrt(numel(reg_inc_plt{iH})));
            
            for iR = 1:numel(reg_inc_plt{iH})
                
                if cfg.plt_fig; figure('Visible','off','Position',[0 0 1080 1080]); hold on; end
                
                reg_plt_loc = find(strcmpi(thk_nme_hld{iH},reg_inc_plt{iH}{iR}));
                
                if ~cfg.plt_fig; subplot(plt_dim,plt_dim,iR); hold on; end
                
                title(mmil_spec_char(reg_inc_plt{iH}{iR},{'_'}),'FontSize',9)
                set(gca,'ylim',[0 ylm_max]);
                set(gca,'xlim',[cfg.bin_cut(1)-3 cfg.bin_cut(end)+3]);
                
                if ~cfg.plt_fig
                    set(gca,'xticklabel',[])
                    set(gca,'yticklabel',[])
                end
                
                % Control Variance Plot
                if isfield(cfg,'con_spr_plt') && cfg.con_spr_plt
                    for iBC = 1:numel(cfg.bin_cut)
                        pct_plt(1,iBC) = prctile(dta_reg_con_spr_hld{iH}(:,reg_plt_loc,iBC),100-((cfg.pvl_lvl/2)*100));
                        pct_plt(2,iBC) = prctile(dta_reg_con_spr_hld{iH}(:,reg_plt_loc,iBC),(cfg.pvl_lvl/2)*100);
                    end
                    patch([cfg.bin_cut fliplr(cfg.bin_cut)],[pct_plt(1,:) fliplr(pct_plt(2,:))],cfg.grp_col{iC}{strcmpi(cfg.grp_nme{iC},cfg.con_nme{iC})},'FaceAlpha',0.2)
                end
                
                % Significance Plot
                %             if isfield(cfg,'con_typ') && ~isempty(cfg.con_typ)
                %
                %                 if  strcmpi(cfg.con_typ,'con_shf')
                %
                %                     for iSG = 1:numel(cfg.cmp_nme)
                %
                %                         cmp_loc = strcmpi(cfg.grp_nme,cfg.cmp_nme{iSG});
                %
                %                         for iBC = 1:numel(cfg.bin_cut)
                %                             if dta_reg_hld{iH}{cmp_loc}(reg_plt_loc,iBC) < pct_plt(2,iBC) || ...
                %                                     pct_plt(1,iBC) < dta_reg_hld{iH}{cmp_loc}(reg_plt_loc,iBC)
                %                                 scatter(cfg.bin_cut(iBC),ylm_cel(iSG),250,rgb('dark grey'),'o','filled'); hold on;
                %                                 if isfield(cfg,'grp_col_hgh_lgh')
                %                                     scatter(cfg.bin_cut(iBC),ylm_cel(iSG),150,cfg.grp_col(strcmpi(cfg.grp_nme,cfg.cmp_nme{iSG}),:),'o','filled','MarkerEdgeColor',cfg.grp_col_hgh_lgh(cmp_loc,:),'LineWidth',1.5);
                %                                 else
                %                                     scatter(cfg.bin_cut(iBC),ylm_cel(iSG),150,cfg.grp_col(strcmpi(cfg.grp_nme,cfg.cmp_nme{iSG}),:),'o','filled');
                %                                 end
                %                             end
                %                         end
                %                     end
                %                 end
                %
                %             end
                
                % Line Plot
                for iP = 1:numel(cfg.grp_nme{iC}) % Line Plot
                    if isfield(cfg,'plt_fig') && cfg.plt_fig == 1
                        plot(cfg.bin_cut,dta_reg_hld{iH}{plt_loc(iP)}(reg_plt_loc,:),'Color',rgb('dark grey'),'LineWidth',7.75);
                        if isfield(cfg,'grp_col_hgh_lgh')
                            plot(      cfg.bin_cut , dta_reg_hld{iH}{iP}(reg_plt_loc,:),'Color',cfg.grp_col_hgh_lgh{iC}{iP},'LineWidth',6.5);
                        end
                        leg(iP) = plot(cfg.bin_cut , dta_reg_hld{iH}{iP}(reg_plt_loc,:),'Color',cfg.grp_col{iC}{iP},'LineWidth',4.25);
                    elseif plt_dim<=3
                        leg(iP) = plot(cfg.bin_cut , dta_reg_hld{iH}{iP}(reg_plt_loc,:),'Color',cfg.grp_col{iC}{iP},'LineWidth',3);
                    else
                        leg(iP) = plot(cfg.bin_cut , dta_reg_hld{iH}{iP}(reg_plt_loc,:),'Color',cfg.grp_col{iC}{iP},'LineWidth',1);
                    end
                end
                
                if cfg.plt_fig
                    tightfig(); set(gcf,'Position',[0 0 1080 1080]);
                    cfg = [];
                    cfg.fle_nme = [cfg.sve_loc '/' cfg.sve_nme{iC} '_' 'region_line_plot' '_' reg_inc_plt{iH}{iR} '_' hms_nme{iH} ];
                    cfg.prn_typ = 'eps';
                    cfg.prn_typ_ext = 'png';
                    mmil_print_plot2(cfg)
                    close all
                end
                
            end
            
            if ~cfg.plt_fig
                tightfig(); set(gcf,'Position',[0 0 1080 1080]);
                print(gcf,[cfg.sve_loc{iC} '/' cfg.sve_nme{iC} '_' 'region_line_plot' '_' hms_nme{iH} '.png'],'-dpng','-r200')
                close all
            end
            
        end
        
    end
end

%% Surface Plot
if isfield(cfg,'reg_srf_plt') && cfg.reg_srf_plt
    
    plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};    
    
    for iGO = 1:numel(cfg.grp_ovr)
        
        % Regional Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Views
        sph     = {'lh'  'rh'};
        sph_vew = {'lat' 'med' };
        
        % Load Surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~strcmpi( cfg.fsr_nme , 'fsaverage' )
            sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.fsr_dir));
            sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.fsr_nme '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
        else
            sbj_dir_lst = dir(sprintf('%s/*',cfg.fsr_dir));
            sbj_dir_lst = regexp({sbj_dir_lst.name},['' cfg.fsr_nme],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
        end
        
        srf_brn{1} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'lh.pial']);
        srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
        srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
        srf_brn{2} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'rh.pial']);
        srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
        srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;
        
        % Data Setup
        dta_reg_hld{1} = cellfun(@(x) x(string_find(cfg.thk_nme,'lhs'),:),cfg.grp_dta.dta_reg{iGO},'uni',0);
        dta_reg_hld{2} = cellfun(@(x) x(string_find(cfg.thk_nme,'rhs'),:),cfg.grp_dta.dta_reg{iGO},'uni',0);
        
        thk_nme_hld{1} = cfg.thk_nme(string_find(cfg.thk_nme,'lhs'));
        thk_nme_hld{2} = cfg.thk_nme(string_find(cfg.thk_nme,'rhs'));
        
        % Colorbar
        top_pct = 1;
        col{1} = rgb('orange');
        col{2} = rgb('neon red');
        col{3} = rgb('red');
        col{4} = rgb('dark red');
        col{5} = rgb('medium grey')-0.15;
        col{6} = rgb('dark blue');
        col{7} = rgb('blue');
        col{8} = rgb('bright blue');
        col{9} = rgb('cyan');
        col_map = [];
        for iC = 1:numel(col)-1
            col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
        end
        col_map = flipud(col_map);
           
        for iCM = 1:numel( cfg.cmp_nme{iGO} )
            
            % Plot
            fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
            
            cmp_loc = find( strcmpi( cfg.grp_nme{iGO} , cfg.cmp_nme{iGO}{iCM} ) );
            
            for iH = 1:numel(sph)
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,albl,~]=fs_read_annotation([ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' sph{iH} '.aparc'  cfg.prc_nme '.annot']);
                albl_hld = albl;
                albl = cellfun( @(x) mmil_spec_char( x , {'-'}) , albl , 'uni' , 0 );
                pct_hld = [albl num2cell(nan(size(albl,1),1))];
                                                
                for iF = 1:size(pct_hld,1)
                    
                    roi_loc = find(strcmpi( cfg.thk_nme , [ sph{iH} 's' '_'  pct_hld{iF,1} ] ));
                    
                    con_hld = nanmean( squeeze(cfg.grp_dta.con_spr_reg{iGO}( : , roi_loc , cfg.reg_col )),2);
                    dta_hld = nanmean( cfg.grp_dta.dta{iGO}{cmp_loc}( roi_loc , cfg.reg_col ));
                    
                    if ~isempty(roi_loc)
                        
                        pct_plt(1) = prctile( con_hld , 100-((cfg.pvl_lvl/2)*100));
                        pct_plt(2) = prctile( con_hld , (cfg.pvl_lvl/2)*100);
                        
                        if pct_plt(1) < dta_hld || ...
                           pct_plt(2) > dta_hld
                            
                            men_roi = nanmean( con_hld );
                            std_roi = nanstd( con_hld );
                            pct_hld{iF,2} = ( dta_hld - men_roi ) / std_roi;
                            
                        else
                            
                            pct_hld{iF,2} = 0;
                            
                        end
                        
                    else
                        pct_hld{iF,2} = nan;
                    end
                end
                
                pct_hld( cell2mat(pct_hld(:,2))>5 , 2 )  = {10};
                pct_hld( cell2mat(pct_hld(:,2))<-5 , 2 ) = {-10};
                
                pct_hld(:,2) = num2cell( (cell2mat(pct_hld(:,2)) / 23) + .5);
                
                pct_hld(:,1) = albl_hld;
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isfield(cfg,'reg_inc_plt') && ~isempty(cfg.reg_inc_plt)
                    iIN = ismember(pct_hld(:,1),unique(cellfun(@(x) x(4:end),cfg.reg_inc_plt,'uni',0)) );
                    pct_hld = pct_hld(iIN,:);
                end
                
                for iSP = 1:numel(sph_vew)
                    
                    if cfg.plt_fig; fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]); end
                    
                    pcfg = [];
                    
                    pcfg.surf_brain  = srf_brn{iH};
                    pcfg.aparc       = [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' sph{iH} '.aparc'  cfg.prc_nme '.annot'];
                    
                    pcfg.sph         = sph{iH};
                    pcfg.sph_vew     = sph_vew{iSP};
                    
                    pcfg.label       = 0;
                    pcfg.radius      = [];
                    pcfg.alpha       = 1;
                    
                    pcfg.non_ele     = [];
                    pcfg.sve_img     = 0; % ###
                    
                    pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
                    
                    pcfg.fig_hdl = fig_hld(1);
                    
                    pcfg.col_map = col_map;
                    
                    pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
                    pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
                    
                    pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
                    
                    pcfg.top_pct = 1;
                    
                    nyu_plot2(pcfg);
                    
                end
            end
                       
            % Add Colorbar
            axes('OuterPosition',[.05 .05 .9 .05],'visible','off')
            colormap(col_map)
            clb = colorbar('south','Position',[.05 .05 .9 .05]);
            clb.TickLength = 0;
            clb.TickLabels = cellfun(@num2str,num2cell(roundsd(linspace(-5,5,11),2)),'uni',0);
            ylabel(clb,cfg.ttl,'FontSize',20)
            
            % Save
            print(gcf,[ cfg.sve_loc{iGO} '/' cfg.sve_nme{iGO} '_' 'region_surface_plot' '_' cfg.cmp_nme{iGO}{iCM} '.png'],'-dpng','-r200')
            close all
            
        end
    end
end

end





