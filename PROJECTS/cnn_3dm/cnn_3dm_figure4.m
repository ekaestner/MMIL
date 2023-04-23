
plt_out_dir = [ out_dir '/' 'Figures' '/' 'Figure4' '/']; ejk_chk_dir(plt_out_dir);

load([ dta_dir '/' 'performance.mat' ])

%% Box plot
ylm_hld = { [ 60 95 ] [ 60 95 ] [ 50 100 ] [50 95 ] [50 100 ] [ 60 95 ] [ 50 100 ] [ 50 100 ] [ 50 100 ] };

for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { cell2mat(dta_hld.fcn_2dm_300.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_3dm_300.(mes_int_nme{iM})(:,3)) ...
        cell2mat(dta_hld.fcn_2dm_750.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_3dm_750.(mes_int_nme{iM})(:,3)) ...
        cell2mat(dta_hld.fcn_2dm_fll.(mes_int_nme{iM})(:,3)) cell2mat(dta_hld.fcn_3dm_fll.(mes_int_nme{iM})(:,3)) };
    fcfg.xdt = { 0.75 1.25 2.75 3.25 4.75 5.25 };
    
    fcfg.fce_col     = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_300')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_300')) ...
        mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_750')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_750')) ...
        mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_fll')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_fll')) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_300')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_300')) ...
        mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_750')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_750')) ...
        mdl_int_col(strcmpi(mdl_int_nme,'fcn_2dm_fll')) mdl_int_col(strcmpi(mdl_int_nme,'fcn_3dm_fll')) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { 'fcn_2dm_300' 'fcn_3dm_300' 'fcn_2dm_750' 'fcn_3dm_750' 'fcn_2dm_fll' 'fcn_3dm_fll' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = ylm_hld{iM};
    
    fcfg.jtr_wdt = 0.10;
    fcfg.box_wdt = 0.20;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure4_boxplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end

%% Difference Plot
for iM = 1:numel(mes_int_nme)

    fcfg = [];
    
    fcfg.ydt = { -1*(cell2mat(dta_hld.fcn_2dm_300.(mes_int_nme{iM})(:,3))-cell2mat(dta_hld.fcn_3dm_300.(mes_int_nme{iM})(:,3))) ...
                 -1*(cell2mat(dta_hld.fcn_2dm_750.(mes_int_nme{iM})(:,3))-cell2mat(dta_hld.fcn_3dm_750.(mes_int_nme{iM})(:,3))) ...
                 -1*(cell2mat(dta_hld.fcn_2dm_fll.(mes_int_nme{iM})(:,3))-cell2mat(dta_hld.fcn_3dm_fll.(mes_int_nme{iM})(:,3))) };
    fcfg.xdt = { 1 3 5 };
    
    fcfg.fce_col     = [ repmat({[0.2 0.2 0.2]},1,numel(fcfg.ydt)) ];
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
    fcfg.box_plt_col = [ repmat({[0.7 0.7 0.7]},1,numel(fcfg.ydt)) ];
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { '300' '750' 'fll' };
    fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
    fcfg.ylb = mes_int_nme(iM);
    fcfg.ylm = [ -round(max(cellfun(@max,fcfg.ydt)))-2 round(max(cellfun(@max,fcfg.ydt)))+2];
    
    fcfg.hln     = 0;
    fcfg.jtr_wdt = 0.10;
    fcfg.box_wdt = 0.20;
    
    fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.out_dir = plt_out_dir;
    fcfg.out_nme = [ 'figure4_diffplot' '_' num2str(iM) '_' mes_int_nme{iM} ];
    
    ejk_scatter(fcfg)

end
