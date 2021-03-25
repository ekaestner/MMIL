% fcfg = [];
%
% fcfg.grp_nme = phe_grp;
% fcfg.grp_mbs = grp_mbm;
% fcfg.thk_nme = phe_wht_nme(:,15:end);
% fcfg.thk_dta = cell2mat(phe_wht_dta(:,15:end));
%
% fcfg.clc_con_smp_mtx = 1;
% fcfg.con_smp_nme = pvl_con;
% fcfg.con_smp_rep = num_rep;
%
% fcfg.clc_dff_mtx = 1;
% fcfg.con_dff_nme = pvl_con;
% fcfg.cmp_dff_nme = pvl_cmp;
% fcfg.dff_rep     = num_rep;

function [ cor_mtx , con_smp , dff_smp ] = GraphThickMatrix( cfg )

%%
fcfg = [];
fcfg.grp_fle = cfg.sbj_grp;
grp_fle      = ejk_collate_groups( fcfg );

%% Calculate Group Matrices
for iD = 1:numel(cfg.sbj_grp_col)
    
    cor_mtx{iD} = repmat( { zeros( size(cfg.roi_nme,2) , size(cfg.roi_nme,2)) } , 1 , numel(cfg.sbj_grp_nme{iD}) );
    
    for iG = 1:numel(cfg.sbj_grp_nme{iD})
        
        grp_ind = grp_fle.table.(cfg.sbj_grp_col{iD}){ strcmpi( grp_fle.table.(cfg.sbj_grp_col{iD})(:,1) , cfg.sbj_grp_nme{iD}{iG} ) , 2 };
        sbj_row = grp_fle.(cfg.sbj_grp_col{iD}) == grp_ind;
        
        for iR = 1:size(cfg.roi_nme,2)
            for iC = 1:size(cfg.roi_nme,2)
                
                if iR ~= iC
                    
                    hld = corrcoef( cfg.dta( sbj_row , iR ) , cfg.dta( sbj_row , iC ) );
                    cor_mtx{iD}{iG}(iR,iC) = hld(1,2);
                    
                else
                    
                    cor_mtx{iD}{iG}(iR,iC) = nan;
                    
                end
                
            end
        end
    end
    
end

%% Make Plot
if isfield(cfg,'plt') && cfg.plt
    
    ejk_chk_dir(cfg.out_dir);
    
    for iD = 1:numel(cfg.sbj_grp_col)
        
        ejk_chk_dir([ cfg.out_dir '/' cfg.sbj_grp_col{iD} ]);
        ejk_chk_dir([ cfg.out_dir '/' cfg.sbj_grp_col{iD} '/' 'covariance_plots' ]);
        
        tot_nme_num = ceil(sqrt(numel(cfg.roi_nme)));
        
        for iG = 1:numel(cfg.sbj_grp_nme{iD})
            
            ejk_chk_dir([ cfg.out_dir '/' cfg.sbj_grp_col{iD} '/' 'covariance_plots' '/' cfg.sbj_grp_nme{iD}{iG} ]);
            
            grp_ind = grp_fle.table.(cfg.sbj_grp_col{iD}){ strcmpi( grp_fle.table.(cfg.sbj_grp_col{iD})(:,1) , cfg.sbj_grp_nme{iD}{iG} ) , 2 };
            sbj_row = grp_fle.(cfg.sbj_grp_col{iD}) == grp_ind;
            
            ttt = cfg.dta( sbj_row , : );
            dta_max = roundsd(max(ttt(:)),2);
            dta_min = roundsd(min(ttt(:)),2);
            
            for iR = 1:size(cfg.roi_nme,2)
                
                figure('Visible','off','Position',[0 0 1920 1080])
                
                for iC = 1:size(cfg.roi_nme,2)
                    if iR ~= iC
                        
                        subplot(tot_nme_num,tot_nme_num,iC)
                        [hld,hld_pvl] = corrcoef(cfg.dta(sbj_row,iR),cfg.dta(sbj_row,iC));
                        
                        scatter(cfg.dta(sbj_row,iR),cfg.dta(sbj_row,iC),'filled','MarkerFaceColor',rgb('black'),'MarkerEdgeColor',rgb('grey'))
                        
                        ply_fit = polyfit(cfg.dta(sbj_row,iR),cfg.dta(sbj_row,iC),1);
                        ply_fit = ply_fit(1)*cfg.dta(sbj_row,iR)+ply_fit(2);
                        hold on; plot(cfg.dta(sbj_row,iR),ply_fit,'r');
                        
                        set(gca,'xlim',[dta_min dta_max]);
                        set(gca,'ylim',[dta_min dta_max]);
                        title(cfg.roi_nme{iC})
                        set(gca,'xticklabel',[])
                        set(gca,'yticklabel',[])
                        
                        if hld_pvl(1,2)>.05
                            text(dta_min,dta_max-dta_max*0.10,num2str(roundsd(hld(1,2),2)),'Color',rgb('black'))
                        elseif hld_pvl(1,2)<=.05
                            text(dta_min,dta_max-dta_max*0.10,num2str(roundsd(hld(1,2),2)),'Color',rgb('red'))
                        end
                    else
                        subplot(tot_nme_num,tot_nme_num,iC)
                        axis off
                        text(0.1,0.5,cfg.roi_nme{iR},'FontSize',16)
                    end
                    
                end
                
                tightfig();
                set(gcf,'Position',[0 0 1920 1080]);
                
                print(gcf, [ cfg.out_dir '/' cfg.sbj_grp_col{iD} '/' 'covariance_plots' '/' cfg.sbj_grp_nme{iD}{iG} '/' cfg.roi_nme{iR} '.png'] , '-dpng' , '-r200' )
                close all
                
            end
        end
    end
end

%% Calculate Control Resample Matrix
if cfg.clc_con_smp_mtx
    
    for iD = 1:numel(cfg.sbj_grp_col)
        
        grp_ind     = grp_fle.table.(cfg.sbj_grp_col{iD}){ strcmpi( grp_fle.table.(cfg.sbj_grp_col{iD})(:,1) , cfg.con_smp_nme{iD} ) , 2 };
        sbj_row     = find( grp_fle.(cfg.sbj_grp_col{iD}) == grp_ind );
        tot_num_con = numel( sbj_row );
        
        num_rsm     = tot_num_con-cfg.lve_num_out;
        
        rep_mtx     = zeros( cfg.con_smp_rep , num_rsm );
        sbj_row_smp = zeros( cfg.con_smp_rep , num_rsm );
        con_smp_mtx = zeros( cfg.con_smp_rep , numel(cfg.roi_nme) , numel(cfg.roi_nme) );
        
        tic
        for iRP = 1:cfg.con_smp_rep
            
            rep_mtx(iRP,:)     = randsample(tot_num_con,num_rsm);
            sbj_row_smp(iRP,:) = sbj_row( rep_mtx(iRP,:) );
            
            for iR = 1:size(cfg.roi_nme,2)
                for iC = 1:size(cfg.roi_nme,2)
                    if iR ~= iC
                        hld = corrcoef(cfg.dta( sbj_row_smp(iRP,:) , iR ),cfg.dta( sbj_row_smp(iRP,:) , iC ));
                        con_smp_mtx(iRP,iR,iC) = hld(1,2);
                    else
                        con_smp_mtx(iRP,iR,iC) = nan;
                    end
                end
            end
            
            
            switch iRP
                case 250
                    fprintf('Finished %i : %f\n',250,toc);
                case 500
                    fprintf('Finished %i : %f\n',500,toc);
                case 750
                    fprintf('Finished %i : %f\n',750,toc);
                case 1000
                    fprintf('Finished %i : %f\n',1000,toc);
                case 1250
                    fprintf('Finished %i : %f\n',1250,toc);
                case 1500
                    fprintf('Finished %i : %f\n',1500,toc);
                case 1750
                    fprintf('Finished %i : %f\n',1750,toc);
                case 2000
                    fprintf('Finished %i : %f\n',2000,toc);
                case 2250
                    fprintf('Finished %i : %f\n',2250,toc);
                case 2500
                    fprintf('Finished %i : %f\n',2500,toc);
                case 2750
                    fprintf('Finished %i : %f\n',2750,toc);
                case 3000
                    fprintf('Finished %i : %f\n',3000,toc);
                case 3250
                    fprintf('Finished %i : %f\n',3250,toc);
                case 3500
                    fprintf('Finished %i : %f\n',3500,toc);
                case 3750
                    fprintf('Finished %i : %f\n',3750,toc);
                case 4000
                    fprintf('Finished %i : %f\n',4000,toc);
                case 4250
                    fprintf('Finished %i : %f\n',4250,toc);
                case 4500
                    fprintf('Finished %i : %f\n',4500,toc);
                case 4750
                    fprintf('Finished %i : %f\n',4750,toc);
                case 5000
                    fprintf('Finished %i : %f\n',5000,toc);
            end
        end
        
        con_smp.tme = toc;
        
        con_smp.con_smp_mtx{iD} = con_smp_mtx;
        con_smp.sbj_row_smp{iD} = sbj_row_smp;
        
        clear rep_mtx sbj_row_smp con_smp_mtx toc
        
    end
end

%% Calculate Stat Resamples
if cfg.clc_dff_mtx
    
    error('fix me! ;P')
    
    for iD = 1:numel(cfg.cmp_dff_nme)
        
        men_num_pat = [];%ceil(mean(cellfun(@numel,cfg.grp_mbs(~iCO(1:end-1)))));
        
        iCM = strcmpi(cfg.grp_nme,cfg.cmp_dff_nme{iD});
        tot_num_cmp = numel(cfg.grp_mbs{iCM});
        
        tot_num_mbm = men_num_pat + tot_num_cmp;
        
        fke_con = zeros(cfg.dff_rep,men_num_pat);
        fke_cmp = zeros(cfg.dff_rep,tot_num_cmp);
        
        con_dff_mtx = repmat({zeros(size(cfg.thk_nme,2),size(cfg.thk_nme,2))},1,cfg.con_smp_rep);
        cmp_dff_mtx = repmat({zeros(size(cfg.thk_nme,2),size(cfg.thk_nme,2))},1,cfg.con_smp_rep);
        
        tic
        for iRP = 1:cfg.dff_rep
            
            dff_hld = [sbj_row_smp(iRP,:) cfg.grp_mbs{iCM}'];
            
            rep_mtx        = randperm(tot_num_mbm);
            
            fke_con(iRP,:) = dff_hld(rep_mtx(1:men_num_pat));
            fke_cmp(iRP,:) = dff_hld(rep_mtx(men_num_pat+1:end));
            
            for iR = 1:size(cfg.thk_nme,2)
                for iC = 1:size(cfg.thk_nme,2)
                    if iR ~= iC
                        hld = corrcoef(cfg.thk_dta(fke_con(iRP,:),iR),cfg.thk_dta(fke_con(iRP,:),iC));
                        con_dff_mtx{iRP}(iR,iC) = hld(1,2);
                        hld = corrcoef(cfg.thk_dta(fke_cmp(iRP,:),iR),cfg.thk_dta(fke_cmp(iRP,:),iC));
                        cmp_dff_mtx{iRP}(iR,iC) = hld(1,2);
                    else
                        con_dff_mtx{iRP}(iR,iC) = nan;
                        cmp_dff_mtx{iRP}(iR,iC) = nan;
                    end
                end
            end
            
            switch iRP
                case 250
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},250,toc);
                case 500
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},500,toc);
                case 750
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},750,toc);
                case 1000
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},1000,toc);
                case 1250
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},1250,toc);
                case 1500
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},1500,toc);
                case 1750
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},1750,toc);
                case 2000
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},2000,toc);
                case 2250
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},2250,toc);
                case 2500
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},2500,toc);
                case 2750
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},2750,toc);
                case 3000
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},3000,toc);
                case 3250
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},3250,toc);
                case 3500
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},3500,toc);
                case 3750
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},3750,toc);
                case 4000
                    fprintf('Finished %s %i : %f\n',cfg.cmp_dff_nme{iD},4000,toc);
            end
            
        end
        
        dff_smp.(mmil_spec_char(cfg.cmp_dff_nme{iD},{' ' '&'})).tme = toc;
        
        dff_smp.(mmil_spec_char(cfg.cmp_dff_nme{iD},{' ' '&'})).con_dff_mtx = con_dff_mtx;
        dff_smp.(mmil_spec_char(cfg.cmp_dff_nme{iD},{' ' '&'})).fke_con     = fke_con;
        dff_smp.(mmil_spec_char(cfg.cmp_dff_nme{iD},{' ' '&'})).cmp_dff_mtx = cmp_dff_mtx;
        dff_smp.(mmil_spec_char(cfg.cmp_dff_nme{iD},{' ' '&'})).fke_cmp     = fke_cmp;
        
        clear fke_cmp fke_con con_dff_mtx cmp_dff_mtx toc
        
    end
else
    dff_smp = [];
end

end

















