function dta_out = ejk_lpa( cfg )

tmp_out_dir = [ cfg.out_dir '/' ];

%% Save Data
dep_var.sbj_nme = cfg.sbj_nme;
for iD = 1:numel(cfg.dta_nme)
    dep_var.(cfg.dta_nme{iD}) = cfg.dta(:,iD);
end

ejk_chk_dir(tmp_out_dir)
save( [ tmp_out_dir '/' 'dep_var.mat' ], 'dep_var' )

for iM = cfg.mdl_num-1:cfg.mdl_num+1
    ejk_chk_dir([ tmp_out_dir '/' 'model' '_' num2str(iM) ])
end
    
%% Command
sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_lpa.r'')\n',sve_cmd);
sve_cmd = sprintf('%sejk_lpa( ''%s/dep_var.mat'', %i, ''%s'' )', sve_cmd, tmp_out_dir, cfg.mdl_num, tmp_out_dir);
cell2csv( [tmp_out_dir '/example_R_script.r'], {sve_cmd} );
unix( [ 'Rscript ' tmp_out_dir '/example_R_script.r' ] );

%% Load
cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    dta_out{cnt} = mmil_readtext( [ tmp_out_dir '/' 'model' '_' num2str(iM) '/' 'model_' num2str(iM) '.csv' ] );
    dta_out{cnt} = cell2mat( dta_out{cnt}(2:end, 2:end) );
    cnt = cnt+1;
end

%% Setup 
col = distinguishable_colors( cfg.mdl_num+1 );
for iC = 1:size(col,1)
    fce_col{iC} = col(iC,:);
    edg_col{iC} = rgb('black');
end

cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    
    cnf_dta{cnt} = cell( iM, 1 );
    dta_dta{cnt} = cell( numel(cfg.dta_nme), iM);
    
    avg_dta{cnt}     = nan( numel(cfg.dta_nme), iM );
    std_err_dta{cnt} = nan( numel(cfg.dta_nme), iM );
    
    hgh_cnf_avg_dta{cnt}     = nan( numel(cfg.dta_nme), iM );
    hgh_cnf_std_err_dta{cnt} = nan( numel(cfg.dta_nme), iM );
    
    for iCT = 1:iM
        
        cnf_dta{cnt}{iCT,1} = dta_out{cnt}( dta_out{cnt}(:,1)==iCT ,2 );
        hgh_cnf_ind{cnt}{iCT,1} = cnf_dta{cnt}{iCT,1} > .90;
        
        for iDN = 1:numel( cfg.dta_nme )
            
            dta_dta{cnt}{iDN, iCT } = cfg.dta( dta_out{cnt}(:,1)==iCT, iDN);
            
            avg_dta{cnt}(iDN, iCT)     = mean(dta_dta{cnt}{iDN, iCT});
            std_err_dta{cnt}(iDN, iCT) = std(dta_dta{cnt}{iDN, iCT}) / sqrt(numel(dta_dta{cnt}{iDN, iCT}));
            
            hgh_cnf_avg_dta{cnt}(iDN, iCT)     = mean(dta_dta{cnt}{iDN, iCT}(hgh_cnf_ind{cnt}{iCT,1}));
            hgh_cnf_std_err_dta{cnt}(iDN, iCT) = std(dta_dta{cnt}{iDN, iCT}(hgh_cnf_ind{cnt}{iCT,1})) / sqrt(numel(dta_dta{cnt}{iDN, iCT}(hgh_cnf_ind{cnt}{iCT,1})));
            
        end
    end
    
    cnt = cnt+1;
end

%% Confidence Scatter
%
cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    
    fcfg = [];
    
    fcfg.ydt     = cnf_dta{cnt}';
    fcfg.xdt     = num2cell(1:iM);
    
    fcfg.fce_col = fce_col(1:iM);
    fcfg.edg_col = edg_col(1:iM);
    
    fcfg.xlb = { strcat( 'Group', cellfun( @num2str, fcfg.xdt,'uni',0) )};
    fcfg.ylb = 'Confidence';
    
    fcfg.ylm = [-0.1 1.1];
    
    fcfg.out_dir = [ tmp_out_dir '/'  'model' '_' num2str(iM) '/'];
    fcfg.out_nme = 'ConfidenceResults';
    
    ejk_scatter(fcfg)
    
    cnt = cnt+1;
end

%% Omniubs Graph
cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    
    figure('Visible','off','Position',[0 0 1080 1080]); hold on;
    for iCT = 1:iM
        %     plot( 1:numel( cfg.dta_nme ), avg_dta(:,iCT)', 'Color', fce_col{iCT} )
        errorbar( avg_dta{cnt}(:,iCT)', std_err_dta{cnt}(:,iCT)', 'Color', fce_col{iCT} )
    end
    legend('Location','north','Orientation','horizontal')
    set(gcf,'PaperSize',[80 80],'PaperUnits','inches');
    xlim([0 numel( cfg.dta_nme )+1]);
    xticklabels( [{''} cfg.dta_nme {''}])
    tightfig();
    print( gcf, [ tmp_out_dir '/' 'model' '_' num2str(iM) '/' 'Omnibus_all.png'], '-dpng' );
    close all
    
    cnt = cnt+1;
end

%% Omniubs Graph - High confidence
cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    
    figure('Visible','off','Position',[0 0 1080 1080]); hold on;
    for iCT = 1:iM
        %     plot( 1:numel( cfg.dta_nme ), avg_dta(:,iCT)', 'Color', fce_col{iCT} )
        errorbar( hgh_cnf_avg_dta{cnt}(:,iCT)', hgh_cnf_std_err_dta{cnt}(:,iCT)', 'Color', fce_col{iCT} )
    end
    legend('Location','north','Orientation','horizontal')
    set(gcf,'PaperSize',[80 80],'PaperUnits','inches');
    xlim([0 numel( cfg.dta_nme )+1]);
    tightfig();
    print( gcf, [ tmp_out_dir '/' 'model' '_' num2str(iM) '/' 'Omnibus_high_confidence.png'], '-dpng' );
    close all
    
    cnt = cnt+1;
end

%% Region Graph
cnt = 1;
for iM = cfg.mdl_num-1:cfg.mdl_num+1
    
    ejk_chk_dir([tmp_out_dir '/' 'model' '_' num2str(iM) '/' 'Scatters' '/'])
    
    for iDN = 1:numel( cfg.dta_nme )
        
        %
        ylm_hgh = max(cat(1,dta_dta{cnt}{iDN, :})) + ( max(cat(1,dta_dta{cnt}{iDN, :})) * .05);
        ylm_low = min(cat(1,dta_dta{cnt}{iDN, :})) - abs( min(cat(1,dta_dta{cnt}{iDN, :})) * .05);
        
        %
        fcfg = [];
        
        fcfg.ydt     = dta_dta{cnt}( iDN, :)';
        fcfg.xdt     = num2cell(1:iM);
        
        fcfg.fce_col = fce_col;
        fcfg.edg_col = edg_col;
        
        fcfg.xlb = { strcat( 'Group', cellfun( @num2str, fcfg.xdt,'uni',0) )};
        fcfg.ylb = 'Confidence';
        
        fcfg.ylm = [ylm_low ylm_hgh];
        
        fcfg.ttl = cfg.dta_nme{iDN};
        
        fcfg.out_dir = [tmp_out_dir '/' 'model' '_' num2str(iM) '/' 'Scatters' '/'];
        fcfg.out_nme = cfg.dta_nme{iDN};
        
        ejk_scatter(fcfg)
        
    end
   
    cnt = cnt+1;
end

end