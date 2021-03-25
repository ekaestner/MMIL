% cfg = [];
% 
% cfg.grp_nme = phe_grp; % 
% cfg.bin_dir = 'pos'; % 
% cfg.bin_cut_off = bin_cut_off; % 
% 
% cfg.cor_mtx = cor_mtx; % 
% 
% cfg.clc_con_smp_mtx = 1; % 
% cfg.con_smp = con_smp; % 
% 
% cfg.clc_dff_mtx = 1; % 
% cfg.dff_smp = dff_smp; % 

function [ cor_mtx_bin , con_smp_bin , dff_smp_bin ] = GraphBinzarize( cfg )

%% Binarize Group Networks
for iD = 1:numel(cfg.sbj_grp_col)
    for iG = 1:numel(cfg.sbj_grp_nme{iD})
        
        msk = tril(true(size(squeeze(cfg.cor_mtx{iD}{1}))),-1);

        cor_mtx_bin{iD}{iG} = nan( [ numel(cfg.bin_cut_off) size(squeeze(cfg.cor_mtx{iD}{1}))] );
        
        if strcmpi(cfg.bin_dir,'pos')
            for iBC = 1:numel(cfg.bin_cut_off)
                
                cor_mtx_bin{iD}{iG}(iBC,:,:) = cfg.cor_mtx{iD}{iG};
                
                cor_mtx_hld = cor_mtx_bin{iD}{iG}(iBC,:,:);
                lim = prctile( cor_mtx_hld(msk) , 100-cfg.bin_cut_off(iBC) );
                
                cor_mtx_hld( isnan(cor_mtx_hld)  ) = 0;
                cor_mtx_hld( cor_mtx_hld<lim(1)  ) = 0;
                cor_mtx_hld( cor_mtx_hld>=lim(1) ) = 1;
                
                cor_mtx_bin{iD}{iG}(iBC,:,:) = cor_mtx_hld;
                
            end
        end
    end
    
end

%% Binarize Control Spread
if cfg.clc_con_smp_mtx
    for iD = 1:numel(cfg.sbj_grp_col)
        
        con_smp_bin{iD} = nan([ numel(cfg.bin_cut_off) size(cfg.con_smp.con_smp_mtx{iD})] );
        
        for iRP = 1:size(cfg.con_smp.con_smp_mtx{iD},1)
            
            msk = tril(true(size(squeeze(cfg.con_smp.con_smp_mtx{iD}(iRP,:,:)))),-1);
            
            if strcmpi(cfg.bin_dir,'pos')
                for iBC = 1:numel(cfg.bin_cut_off)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    cor_mtx_hld = squeeze( cfg.con_smp.con_smp_mtx{iD}( iRP , : , : ) );
                    lim = prctile( cor_mtx_hld(msk) , 100-cfg.bin_cut_off(iBC) );
                    
                    cor_mtx_hld( isnan(cor_mtx_hld)  ) = 0;
                    cor_mtx_hld( cor_mtx_hld<lim(1)  ) = 0;
                    cor_mtx_hld( cor_mtx_hld>=lim(1) ) = 1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    con_smp_bin{iD}( iBC , iRP , : , : ) = cor_mtx_hld;
                    
                end
            end
            
        end
    end
else
    con_smp_bin = [];
end

%% Binzarize Difference Matrices
if cfg.clc_dff_mtx
    
    grp_nme = fieldnames(cfg.dff_smp);
    
    for iG = 1:numel(grp_nme)
        
        dff_smp_bin.(grp_nme{iG}).con_dff_bin = repmat({cell(1,numel(cfg.dff_smp.(grp_nme{iG}).con_dff_mtx))},1,numel(cfg.bin_cut_off));
        dff_smp_bin.(grp_nme{iG}).cmp_dff_bin = repmat({cell(1,numel(cfg.dff_smp.(grp_nme{iG}).con_dff_mtx))},1,numel(cfg.bin_cut_off));
        
        for iRP = 1:numel(cfg.dff_smp.(grp_nme{iG}).con_dff_mtx)
            
            msk_con = tril(true(size(cfg.dff_smp.(grp_nme{iG}).con_dff_mtx{iRP})),-1);
            msk_cmp = tril(true(size(cfg.dff_smp.(grp_nme{iG}).cmp_dff_mtx{iRP})),-1);
            
            if strcmpi(cfg.bin_dir,'pos')
                for iBC = 1:numel(cfg.bin_cut_off)
                    
                    %
                    dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP} = cfg.dff_smp.(grp_nme{iG}).con_dff_mtx{iRP};
                    lim_con = prctile(dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}(msk),100-cfg.bin_cut_off(iBC));
                                        
                    dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}(dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}>lim_con(1))=1;
                    dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}(dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}<lim_con(1))=0;
                    dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}(isnan(dff_smp_bin.(grp_nme{iG}).con_dff_bin{iBC}{iRP}))=0;
                    
                    %
                    dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP} = cfg.dff_smp.(grp_nme{iG}).cmp_dff_mtx{iRP};
                    lim_cmp = prctile(dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}(msk),100-cfg.bin_cut_off(iBC));
                    
                    dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}(dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}>lim_cmp(1))=1;
                    dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}(dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}<lim_cmp(1))=0;
                    dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}(isnan(dff_smp_bin.(grp_nme{iG}).cmp_dff_bin{iBC}{iRP}))=0;
                    
                end
            end
            
        end
        
    end
else
    dff_smp_bin = [];
end

end