

function ejk_grp_stt( cfg , dta )

% Make Stats
for iST = 1:size( dta , 2 )
    
    for iGR = 1:numel(cfg.ind)
        ydt{iGR} = dta( cfg.ind{iGR} , iST );
    end
    
    if ~all( isnan(cat(1,ydt{:})) )
        
        scfg = [];
        
        scfg.dta     = ydt;
        scfg.grp_nme = cfg.ind_nme;
        
        [ ~ , stt{iST} , fll{iST} ] = mmil_roi_anv(scfg); %
    else
        
        stt{iST} = '';
        fll{iST} = nan(1,numel(fll{1}));
        
    end
    
end

%% Save Stats
col_hed = nchoosek(1:numel(cfg.ind),2);
for iCH = 1:size(col_hed,1)
    col_hed_nme{iCH} = [ cfg.ind_nme{col_hed(iCH,1)} '_VS_' cfg.ind_nme{col_hed(iCH,2)} ];
end

ejk_chk_dir(cfg.out_dir)

cell2csv([ cfg.out_dir '/' cfg.out_nme '_' 'PValues' '.csv'] , [ { '' 'ANOVA' col_hed_nme{:} } ; cfg.dta_lbl num2cell(cat(1,fll{:})) ] );
cell2csv([ cfg.out_dir '/' cfg.out_nme '_' 'ANOVA' '.csv']   , [ { '' 'ANOVA' } ; [ cfg.dta_lbl stt(:) ] ]);

clear pnv stt fll