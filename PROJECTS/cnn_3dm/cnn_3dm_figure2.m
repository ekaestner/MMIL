
plt_out_dir = [ out_dir '/' 'Figures' '/' 'Figure2' '/']; ejk_chk_dir(plt_out_dir);

load([ dta_dir '/' 'performance.mat' ])

%% Box plot
ylm_hld = { [ 40 100 ] [ 40 100 ] [ 30 100 ] [ 30 100 ] [ 40 100 ] [ 40 100 ] [ 40 100 ] [ 30 100 ] [ 40 100 ] };

for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { cell2mat(dta_hld.fcn_3dm_org.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_3dm_shf.(mes_int_nme{iM})(:,3)) };
    fcfg.xdt = { 1 2 };
    
    fcfg.fce_col     = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_org')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_shf')) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_org')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_shf')) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { 'fcn_3dm_org' 'fcn_3dm_shf' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+3.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = ylm_hld{iM};
    
    fcfg.jtr_wdt = 0.20;
    fcfg.box_wdt = 0.30;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure2_boxplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end

%% Difference Plot
for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { cell2mat(dta_hld.fcn_3dm_org.(mes_int_nme{iM})(1:5,3))-cell2mat(dta_hld.fcn_3dm_shf.(mes_int_nme{iM})(:,3)) };
    fcfg.xdt = { 1 };
    
    fcfg.fce_col     = [ repmat({[0.2 0.2 0.2]},1,numel(fcfg.ydt)) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ repmat({[0.7 0.7 0.7]},1,numel(fcfg.ydt)) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { 'Diff' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+4.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = [ -round(max(cellfun(@max,fcfg.ydt)))-2 round(max(cellfun(@max,fcfg.ydt)))+2];
    
    fcfg.hln     = 0;    
    fcfg.jtr_wdt = 0.20;
    fcfg.box_wdt = 0.30;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure2_diffplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end
