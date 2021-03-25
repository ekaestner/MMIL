function dta_out = ejk_kmeans( cfg )

tmp_out_dir = [ cfg.out_dir '/' 'Kmeans' '/' ];

%% Save Data
dep_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    dep_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD);
end

ejk_chk_dir(tmp_out_dir)
save( [ tmp_out_dir '/' 'dep_var.mat' ], 'dep_var' )

pos_col = distinguishable_colors(15);

save( [ tmp_out_dir '/' 'pos_col.mat' ], 'pos_col' )

%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_kmeans.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_kmeans( ''%s/dep_var.mat'', ''%s/pos_col.mat'', ''%s'' )', sve_cmd, tmp_out_dir, tmp_out_dir, tmp_out_dir);
cell2csv( [tmp_out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' tmp_out_dir '/example_R_script.r' ] );

%% Load 
dta_out = mmil_readtext( [ tmp_out_dir '/' 'model_cluster.csv' ] );
dta_out = cell2mat( dta_out(2:end, 2:end) );

%% Setup 
max_num = max(dta_out(:,end));
mdl_num = max_num-2:max_num;
col = distinguishable_colors( max_num );
for iC = 1:size(col,1)
   fce_col{iC} = col(iC,:);
   edg_col{iC} = rgb('black');
end



for iTP = 1:size(dta_out,2)
    
    dta_dta{iTP}     = cell( numel(cfg.dta_nme), mdl_num(iTP) );
    
    avg_dta{iTP}     = nan( numel(cfg.dta_nme), mdl_num(iTP) );
    std_err_dta{iTP} = nan( numel(cfg.dta_nme), mdl_num(iTP) );
    
    for iCT = 1:max(dta_out(:,iTP))      
        for iDN = 1:numel( cfg.dta_nme )
            
            dta_dta{iTP}{iDN, iCT } = cfg.dta( dta_out(:,iTP)==iCT, iDN);
            
            avg_dta{iTP}(iDN, iCT)     = mean(dta_dta{iTP}{iDN, iCT});
            std_err_dta{iTP}(iDN, iCT) = std(dta_dta{iTP}{iDN, iCT}) / sqrt(numel(dta_dta{iTP}{iDN, iCT}));
            
        end
    end
end

%% Omniubs Graph
for iTP = 1:size(dta_out,2)

    figure('Visible','off','Position',[0 0 1080 1080]); hold on;
    for iCT = 1:mdl_num(iTP)
        %     plot( 1:numel( cfg.dta_nme ), avg_dta(:,iCT)', 'Color', fce_col{iCT} )
        errorbar( avg_dta{iTP}(:,iCT)', std_err_dta{iTP}(:,iCT)', 'Color', fce_col{iCT} )
    end
    legend('Location','north','Orientation','horizontal')
    set(gcf,'PaperSize',[80 80],'PaperUnits','inches');
    xlim([0 numel( cfg.dta_nme )+1]);
    xticklabels( [{''} cfg.dta_nme {''}])
    tightfig();
    print( gcf, [ tmp_out_dir '/' 'Omnibus_all_' num2str(mdl_num(iTP)) '.png'], '-dpng' );
    close all

end

%% Region Graph
for iTP = 1:size(dta_out,2)
    
    ejk_chk_dir([tmp_out_dir '/' 'Scatters' '_' num2str(mdl_num(iTP)) '/'])
    
    for iDN = 1:numel( cfg.dta_nme )
        
        %
        ylm_hgh = max(cat(1,dta_dta{iTP}{iDN, :})) + ( max(cat(1,dta_dta{iTP}{iDN, :})) * .05);
        ylm_low = min(cat(1,dta_dta{iTP}{iDN, :})) - abs( min(cat(1,dta_dta{iTP}{iDN, :})) * .05);
        
        %
        fcfg = [];
        
        fcfg.ydt     = dta_dta{iTP}( iDN, :)';
        fcfg.xdt     = num2cell(1:mdl_num(iTP));
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = edg_col;
        
        fcfg.xlb = { strcat( 'Group', cellfun( @num2str, fcfg.xdt,'uni',0) )};
        fcfg.ylb = 'Confidence';
        
        fcfg.ylm = [ylm_low ylm_hgh];
        
        fcfg.ttl = cfg.dta_nme;
        
        fcfg.out_dir = [tmp_out_dir '/' 'Scatters' '_' num2str(mdl_num(iTP)) '/'];
        fcfg.out_nme = cfg.dta_nme{iDN};
        
        ejk_scatter(fcfg)
        
    end
    
end

end