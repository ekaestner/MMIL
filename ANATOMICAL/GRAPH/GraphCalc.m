% cfg = [];
%
% cfg.mes = 'PathEfficiency';
%
% cfg.dta = cor_mtx_bin;
%
% cfg.con_spr_stt = 1;
% cfg.con_spr_dta = con_smp_bin;
%
% cfg.dff_cmp_stt = 1;
% cfg.dff_cmp_dta = dff_smp_bin;

function out_dta = GraphCalc(cfg)

%% Data - ????
for iC = 1:numel(cfg.dta)
    
    dst_mtx{iC} = repmat({nan( size(cfg.dta{iC}{1}))},1,numel(cfg.dta{1}));
    dta_hld{iC} = repmat({nan( [ size(cfg.dta{iC}{1},2) size(cfg.dta{iC}{1},1) ] )},1,numel(cfg.dta{1}));
    
    if strcmpi(cfg.mes,'PathEfficiency') || strcmpi(cfg.mes,'CharacteristicPathLength') || strcmpi(cfg.mes,'Transitivity') || strcmpi(cfg.mes,'Modularity')
        if size(cfg.dta{iC}{1},4)>1
            dta(iC) = repmat({nan(size(cfg.dta{iC}{1},3),numel(cfg.dta{iC}))},1,1);
        elseif size(cfg.dta{iC}{1},4)==1
            dta{iC} = repmat({nan( 1 , size(cfg.dta{iC}{1},1) )},1,numel(cfg.dta{iC}));
        end
    elseif strcmpi(cfg.mes,'LocalPathEfficiency') || strcmpi(cfg.mes,'ClusteringCoefficient') || strcmpi(cfg.mes,'Degree') || strcmpi(cfg.mes,'BetweennessCentrality')
        if size(cfg.dta{iC}{1},4)>1
            dta(iC) = repmat({nan(size(cfg.dta{iC}{1},1),numel(cfg.dta{iC}),size(cfg.dta{1}{iC},3))},1,1);
        elseif size(cfg.dta{iC}{1},4)==1
            dta{iC} = repmat( {nan( [ size(cfg.dta{iC}{1},2) size(cfg.dta{iC}{1},1) ])} , 1 , numel(cfg.dta{iC}) );
        end
    end    
end

%% Calculate
for iC = 1:numel(cfg.dta)
    
    for iG = 1:numel(cfg.dta{iC})
        for iBC = 1:size(cfg.dta{iC}{iG},1)
            switch cfg.mes
                case 'PathEfficiency' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                        for iS = 1:size(cfg.dta{iC}{iBC},3)
                            dst_mtx{iC}{iBC}{iS} = distance_bin(squeeze(cfg.dta{iC}{iBC}(:,:,iS)));
                            [~,dta{iC}(iS,iBC)] = charpath(dst_mtx{iC}{iBC}{iS});
                        end
                    elseif size(cfg.dta{iC}{iG},4)==1
                        dst_mtx{iC}{iG}(iBC,:,:) =  distance_bin( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                        [~,dta{iC}{iG}(iBC)] = charpath( squeeze(dst_mtx{iC}{iG}(iBC,:,:)) );
                    end
                case 'CharacteristicPathLength' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                        for iS = 1:size(cfg.dta{iC}{iBC},3)
                            dst_mtx{iC}{iS,iBC} = distance_bin(squeeze(cfg.dta{iC}{iBC}(:,:,iS)));
                            [dta{iC}(iS,iBC),~] = charpath(dst_mtx{iC}{iS,iBC});
                            for iRW = 1:size(dst_mtx{iC}{iS,iBC},1)
                                dst_mtx{iC}{iS,iBC}(iRW,dst_mtx{iC}{iS,iBC}(iRW,:)==inf) = nan;
                                dst_chr_mtx{iC}(iRW,iBC) = nanmean(dst_mtx{iC}{iS,iBC}(iRW,dst_mtx{iC}{iS,iBC}(iRW,:)>0));
                            end
                        end
                    elseif size(cfg.dta{iC}{iG},4)==1
                        dst_mtx{iC}{iG}(iBC,:,:) = distance_bin( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                        [dta{iC}{iG}(iBC),~] = charpath( squeeze(dst_mtx{iC}{iG}(iBC,:,:)) );
                        for iRW = 1:size(dst_mtx{iC}{iG}(iBC,:,:),2)
                            dst_chr_mtx{iC}{iG}(iRW,iBC) = mean(dst_mtx{iC}{iG}( iBC , iRW , dst_mtx{iC}{iG}( iBC , iRW , : )>0));
                        end
                    end
                case 'Transitivity' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                        for iS = 1:size(cfg.dta{iC}{iBC},3)
                            dta{iC}(iS,iBC) = transitivity_bu(squeeze(cfg.dta{iC}{iBC}(:,:,iS)));
                        end
                    elseif size(cfg.dta{iC}{iG},4)==1
                        dta{iC}{iG}(iBC) = transitivity_bu( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                    end
                case 'LocalPathEfficiency' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me');
                    elseif size(cfg.dta{iC}{iG},4)==1
                        [ dta{iC}{iG}(:,iBC) ] = efficiency_bin( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) , 1 );
                    end    
                case 'ClusteringCoefficient'  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                    elseif size(cfg.dta{iC}{iG},4)==1
                        dta{iC}{iG}(:,iBC) = clustering_coef_bu( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                    end
                case 'Modularity'  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                    elseif size(cfg.dta{iC}{iG},4)==1
                        [ dta_hld{iC}{iG}(:,iBC) , dta{iC}{iG}(iBC) ] = modularity_und( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                    end
                case 'Degree'  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if size(cfg.dta{iC}{iG},4)>1
                        error('fix me')
                        for iS = 1:size(cfg.dta{iC}{iBC},3)
                            dta{iC}(:,iBC,iS) = degrees_und(squeeze(cfg.dta{iC}{iBC}(:,:,iS)));
                        end
                    elseif size(cfg.dta{iC}{iG},4)==1
                        dta{iC}{iG}(:,iBC) = degrees_und( squeeze(cfg.dta{iC}{iG}(iBC,:,:)) );
                    end
                case 'BetweennessCentrality'
                    if size(cfg.dta{iC}{iBC},3)>1
                        error('')
                    elseif size(cfg.dta{iC}{iBC},3)==1
                        dta{iC}(:,iBC) = betweenness_bin(cfg.dta{iC}{iBC});
                    end
            end
        end
    end
end

%% Create Output Structure   
switch cfg.mes
    case 'PathEfficiency'
        out_dta.dta = dta;
    case 'LocalPathEfficiency' %
        if size(cfg.dta{iC}{iG},4)>1
            out_dta.dta     = cellfun(@median,dta,'uni',0);
        elseif size(cfg.dta{iC}{iG},4)==1
            out_dta.dta = dta;
        end
        out_dta.dta_reg = dta;
    case 'CharacteristicPathLength'
        out_dta.dta     = dta;
        out_dta.dta_reg = dst_chr_mtx;
    case 'Transitivity'
        out_dta.dta = dta;
    case 'ClusteringCoefficient'
        if size(cfg.dta{iC}{iG},4)>1
            out_dta.dta = cellfun(@median,dta,'uni',0);
        elseif size(cfg.dta{iC}{iG},4)==1
            out_dta.dta = dta;
        end
        out_dta.dta_reg = dta;
    case 'Modularity'
        out_dta.dta     = dta;
        out_dta.dta_cmm = dta_hld;
    case 'Degree'
        if size(cfg.dta{iC}{iG},4)>1
            out_dta.dta = cellfun(@median,dta,'uni',0);
        elseif size(cfg.dta{iC}{iG},4)==1
            out_dta.dta = cellfun( @(x) cellfun(@median,x,'uni',0),dta,'uni',0);
        end
        out_dta.dta_reg = dta;
    case 'BetweennessCentrality'
        if size(cfg.dta{iC}{iG},4)>1
            out_dta.dta = cellfun(@median,dta,'uni',0);
        elseif size(cfg.dta{iC}{iG},4)==1
            
        end
        out_dta.dta_reg = dta;
end

%% CONTROL SPREAD
if cfg.con_spr_stt
    
    %% Data Setup
    for iC = 1:numel(cfg.dta)
        if strcmpi(cfg.mes,'PathEfficiency') || strcmpi(cfg.mes,'CharacteristicPathLength') || strcmpi(cfg.mes,'Transitivity') || strcmpi(cfg.mes,'Modularity')  || strcmpi(cfg.mes,'Modularity')
            con_spr{iC}     = nan(size(cfg.con_spr_dta{iC},2),size(cfg.con_spr_dta{iC},1));
            con_spr_hld{iC} = nan(size(cfg.con_spr_dta{iC},2),size(cfg.con_spr_dta{iC},3),size(cfg.con_spr_dta{iC},1));
        elseif strcmpi(cfg.mes,'LocalPathEfficiency') || strcmpi(cfg.mes,'ClusteringCoefficient') || strcmpi(cfg.mes,'Degree') || strcmpi(cfg.mes,'BetweennessCentrality')
            con_spr{iC}     = nan( size(cfg.con_spr_dta{iC},2) , size(cfg.con_spr_dta{iC},3) , size(cfg.con_spr_dta{iC},1) ); %repmat({zeros(size(cfg.con_spr_dta{1}{1},1),numel(cfg.con_spr_dta{1}))},1,numel(cfg.con_spr_dta));
        end
    end
    
    %% Calculate
    for iC = 1:numel(cfg.dta)
        
        for iBC = 1:size(cfg.con_spr_dta{iC},1)
            for iRP = 1:size(cfg.con_spr_dta{iC},2)
                
                switch cfg.mes
                    case 'PathEfficiency' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        dst_con_mtx = distance_bin( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                        [~,con_spr{iC}(iRP,iBC)] = charpath(dst_con_mtx);
                    case 'CharacteristicPathLength' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        dst_mtx          = distance_bin( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                        con_spr{iC}(iRP,iBC) = charpath(dst_mtx);
                        for iRW = 1:size(dst_mtx,1)
                            con_spr_mtx{iC}(iRP,iRW,iBC) = mean(dst_mtx( iRW , dst_mtx( iRW , : )>0));
                        end
                    case 'Transitivity' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        con_spr{iC}(iRP,iBC) = transitivity_bu( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                    case 'LocalPathEfficiency' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        con_spr{iC}(iRP,:,iBC) = efficiency_bin( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) , 1 );
                    case 'ClusteringCoefficient' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        con_spr{iC}(iRP,:,iBC) = clustering_coef_bu( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                    case 'Modularity' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [ con_spr_hld{iBC}(iRP,:,iBC) , con_spr{iC}(iRP,iBC) ] = modularity_und( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                    case 'Degree' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Surface Covariance Updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        con_spr{iC}(iRP,:,iBC) = degrees_und( squeeze(cfg.con_spr_dta{iC}(iBC,iRP,:,:)) );
                    case 'BetweennessCentrality'
                        con_spr{iBC}(:,iRP) = betweenness_bin(cfg.con_spr_dta{iBC}{iRP});
                end
                
            end
        end
    end
    
    %% Create measure
    switch cfg.mes
        case 'PathEfficiency'
            out_dta.con_spr = con_spr;
        case 'LocalPathEfficiency'
            out_dta.con_spr = cellfun(@median,con_spr,'uni',0);
            out_dta.con_spr_reg = con_spr;
        case 'CharacteristicPathLength'
            out_dta.con_spr{iC}     = con_spr{iC};
            out_dta.con_spr_reg{iC} = con_spr_mtx{iC};
        case 'Transitivity'
            out_dta.con_spr = con_spr;
        case 'ClusteringCoefficient'
            out_dta.con_spr = cellfun(@median,con_spr,'uni',0);
            out_dta.con_spr_reg = con_spr;
        case 'Modularity'
            out_dta.con_spr = con_spr;
            out_dta.con_spr_cmm = con_spr_hld;
        case 'Degree'
            out_dta.con_spr     = cellfun(@squeeze,cellfun(@(x) median(x,2),con_spr,'uni',0),'uni',0);
            out_dta.con_spr_reg = con_spr;
        case 'BetweennessCentrality'
            out_dta.con_spr = cellfun(@median,con_spr,'uni',0);
            out_dta.con_spr_reg = con_spr;
    end
end

%% DIFFERENCE STATS
if cfg.dff_cmp_stt
    
    stt_nme = fieldnames(cfg.dff_cmp_dta);
    
    for iG = 1:numel(stt_nme)
        
        if strcmpi(cfg.mes,'PathEfficiency') || strcmpi(cfg.mes,'CharacteristicPathLength') || strcmpi(cfg.mes,'Transitivity') || strcmpi(cfg.mes,'Modularity')  || strcmpi(cfg.mes,'Modularity')
            out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta = repmat({zeros(1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
            out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta = repmat({zeros(1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
        elseif strcmpi(cfg.mes,'LocalPathEfficiency') || strcmpi(cfg.mes,'ClusteringCoefficient') || strcmpi(cfg.mes,'Degree') || strcmpi(cfg.mes,'BetweennessCentrality')
            out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta = repmat({zeros(1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
            out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta = repmat({zeros(1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
            
            out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta_reg = repmat({zeros(size(cfg.dta{1}{1},1),numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
            out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta_reg = repmat({zeros(size(cfg.dta{1}{1},1),numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{1}))},1,numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin));
            
        end
        
        for iBC = 1:numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin)
            for iRP = 1:numel(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{iBC})
                switch cfg.mes
                    case 'PathEfficiency'
                        
                        dff_con_dst_mtx = distance_bin(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{iBC}{iRP});
                        [~,pth_con_eff] = charpath(dff_con_dst_mtx);
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP) = pth_con_eff;
                        
                        dff_cmp_dst_mtx = distance_bin(cfg.dff_cmp_dta.(stt_nme{iG}).cmp_dff_bin{iBC}{iRP});
                        [~,pth_cmp_eff] = charpath(dff_cmp_dst_mtx);
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP) = pth_cmp_eff;
                        
                    case 'LocalPathEfficiency'
                        
                        error(':p')
                        
                    case 'CharacteristicPathLength'
                        
                    case 'Transitivity'
                        
                        trn_con = transitivity_bu(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP) = trn_con;
                        
                        trn_cmp = transitivity_bu(cfg.dff_cmp_dta.(stt_nme{iG}).cmp_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP) = trn_cmp;
                        
                    case 'ClusteringCoefficient'
                        
                        con_cls = clustering_coef_bu(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP) = median(con_cls);
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta_reg{iBC}(:,iRP) = con_cls;
                        
                        cmp_cls = clustering_coef_bu(cfg.dff_cmp_dta.(stt_nme{iG}).cmp_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP) = median(cmp_cls);
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta_reg{iBC}(:,iRP) = cmp_cls;
                        
                    case 'Modularity'
                        
                        [~, mod_con] = modularity_und(cfg.dff_cmp_dta.(stt_nme{iG}).con_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP) = mod_con;
                        
                        [~, mod_cmp] = modularity_und(cfg.dff_cmp_dta.(stt_nme{iG}).cmp_dff_bin{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP) = mod_cmp;
                        
                    case 'Degree'
                        
                        deg                                                       = degrees_und(cfg.con_spr_dta{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP)       = median(deg);
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta_reg{iBC}(:,iRP) = deg;
                        
                        deg                                                       = degrees_und(cfg.con_spr_dta{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP)       = median(deg);
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta_reg{iBC}(:,iRP) = deg;
                        
                    case 'BetweennessCentrality'
                        
                        btw_cen                                                   = betweenness_bin(cfg.con_spr_dta{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta{iBC}(iRP)       = median(btw_cen);
                        out_dta.cmp_stt.(stt_nme{iG}).con_dff_dta_reg{iBC}(:,iRP) = btw_cen;
                        
                        btw_cen                                                   = betweenness_bin(cfg.con_spr_dta{iBC}{iRP});
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta{iBC}(iRP)       = median(btw_cen);
                        out_dta.cmp_stt.(stt_nme{iG}).cmp_dff_dta_reg{iBC}(:,iRP) = btw_cen;
                end
            end
            
        end
        
    end
    
end

end