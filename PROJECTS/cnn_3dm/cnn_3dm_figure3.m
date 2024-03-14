
plt_out_dir = [ out_dir '/' 'Figures' '/' 'Figure3' '/']; ejk_chk_dir(plt_out_dir);

load([ dta_dir '/' 'performance_flip_reihaneh.mat' ])

%% Box plot
ylm_hld = { [ 75 95 ]  [ 75 95 ] [ 60 100 ] [ 60 100 ] [ 70 100 ] [ 70 100 ] [ 70 100 ] [ 65 100 ] [ 65 100 ] };

for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { cell2mat(dta_hld.fcn_3dm_org.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_3dm_cbt.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_2dm_org.(mes_int_nme{iM})(:,3)) };
    fcfg.xdt = { 1 2 4};
    
    fcfg.fce_col     = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_org')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_cbt')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_org')) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_org')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_cbt')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_org')) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { 'fcn_3dm_org' 'fcn_3dm_shf' '' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+2.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = ylm_hld{iM};
    
    fcfg.jtr_wdt = 0.20;
    fcfg.box_wdt = 0.30;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure3_boxplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end

%% Difference Plot
for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { cell2mat(dta_hld.fcn_3dm_org.(mes_int_nme{iM})(:,3))-cell2mat(dta_hld.fcn_3dm_cbt.(mes_int_nme{iM})(:,3)) ...
                 cell2mat(dta_hld.fcn_3dm_org.(mes_int_nme{iM})(:,3))-cell2mat(dta_hld.fcn_2dm_org.(mes_int_nme{iM})(:,3)) };
    fcfg.xdt = { 1 2 };
    
    fcfg.fce_col     = [ repmat({[0.2 0.2 0.2]},1,numel(fcfg.ydt)) repmat({[0.2 0.2 0.2]},1,numel(fcfg.ydt)) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ repmat({[0.7 0.7 0.7]},1,numel(fcfg.ydt)) repmat({[0.7 0.7 0.7]},1,numel(fcfg.ydt)) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { 'org - cbt' '3dm - 2dm' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+3.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = [ -round(max(cellfun(@max,fcfg.ydt)))-2 round(max(cellfun(@max,fcfg.ydt)))+2];
    
    fcfg.hln     = 0;    
    fcfg.jtr_wdt = 0.20;
    fcfg.box_wdt = 0.30;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure3_diffplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end
