clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm/';

out_dir = [ prj_dir '/' 'Figures/2023_02_17/' ]; % 2022_12_12

dta_dir = [ prj_dir '/' 'Data/2023_02_17/' ]; % 2022_12_12

dta_grp_fle = 'results_Feb17.csv';

%% Load Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_grp_fle ];
[ dta_grp_dta, dta_grp_dta_sbj, dta_grp_dta_col ] = ejk_dta_frm(fcfg);

%% Split Data
mdl_int     = { '3D FCNet (non-harmonized data)' '3D FCNet (harmonized data)' '2D FCNet (non-harmonized data)' '3D FCNet-shuffled labels (non-harmonized data)' };
mdl_int_nme = { 'ThreeD_5Lay_org'                'ThreeD_5Lay_ComBat'         'TwoD_5Lay_org'                  'ThreeD_5lay_shuffled' };
mdl_int_col = { rgb('bright teal')                rgb('light green')          rgb('maroon')                    rgb('grey') };     
 
mes_int     = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };
mes_int_nme = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };

% Add model number
num_mdl = zeros(1,numel(mdl_int)); 
add_mdl_num = [];
for iMDL = 1:numel(mdl_int) 
    num_mdl(iMDL) = sum(strcmpi(dta_grp_dta(:,1),mes_int{1}) & strcmpi(dta_grp_dta(:,2),mdl_int{iMDL})); 
    ttt = repmat(1:num_mdl(iMDL),numel(mes_int_nme),1);
    add_mdl_num = [ add_mdl_num ; ttt(:)];
end

dta_grp_dta     = [ dta_grp_dta num2cell(add_mdl_num)];
dta_grp_dta_col = [ dta_grp_dta_col 'ModelRun' ];

% Take Data
for iMDL = 1:numel(mdl_int)
    for iMES = 1:numel(mes_int)
        
        dta_ind = strcmpi(dta_grp_dta(:,strcmpi(dta_grp_dta_col,'Model')),mdl_int{iMDL}) & ...
                  strcmpi(dta_grp_dta(:,strcmpi(dta_grp_dta_col,'Metric')),mes_int{iMES}); 
        
        dta_hld.(mdl_int_nme{iMDL}).(mes_int_nme{iMES}) = dta_grp_dta( dta_ind, [1 4 3] );
        
    end
    
    % Add Precision
    dta_hld.(mdl_int_nme{iMDL}).Precision = dta_hld.(mdl_int_nme{iMDL}).PPV;
     
    % Add Recall
    dta_hld.(mdl_int_nme{iMDL}).Recall = dta_hld.(mdl_int_nme{iMDL}).Sensitivity;
    
    % Add F1
    prc     = cell2mat(dta_hld.(mdl_int_nme{iMDL}).Precision(:,3));
    rcl     = cell2mat(dta_hld.(mdl_int_nme{iMDL}).Recall(:,3));
    row_sze = size(dta_hld.(mdl_int_nme{iMDL}).Recall(:,3),1);
    dta_hld.(mdl_int_nme{iMDL}).F1 = [ repmat({'F1'},row_sze,1) num2cell([1:row_sze]') num2cell(2 * ((prc.*rcl)./(prc+rcl))) ];
    
end

mes_int     = [ mes_int 'Precision' 'Recall' 'F1' ];
mes_int_nme = mes_int ;

%% Preliminary Figures
plt_nme = { 'Dimensionality' 'Harmonization' 'Figure2' 'Figure3'            };
box_ind = { [ 1 3 ]          [ 1 2 ]         [ 1 4]     [1 2 3]            };
xdt_ind = { [ 1 2 ]          [ 1 2 ]         [ 1 2]     [1 2 4] };
dff_ind = { {[1 3]}          {[1 2] }        []         {[1 2] [1 3] }    };
ylm_use = { [ 70 100 ]       [ 70 100 ]      [ 0 100 ] [ 70 100 ] };
cmp_tbl_ind = [ 1 2 ];

for iPL = 1:numel(plt_nme)
    ejk_chk_dir([ out_dir '/' plt_nme{iPL} '/' ])
    
    for iMES = 1:numel(mes_int)
      
    % Box plots
    clear fcfg
    
    fcfg = [];
    
    for iD = 1:numel(box_ind{iPL})
        fcfg.ydt{iD} = cell2mat(dta_hld.(mdl_int_nme{box_ind{iPL}(iD)}).(mes_int_nme{iMES})(:,3));
        fcfg.xdt{iD} = xdt_ind{iPL}(iD);
    end
    
    fcfg.fce_col     = { mdl_int_col{box_ind{iPL}} };
    fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(box_ind{iPL})) ];
    fcfg.box_plt_col = { mdl_int_col{box_ind{iPL}}  };
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = mdl_int_nme(box_ind{iPL});
    fcfg.xlm = [ 0.5 numel(box_ind{iPL})+5.5 ];
    fcfg.ylb = mes_int_nme(iMES);
    fcfg.ylm = ylm_use{iPL};
    
    fcfg.mkr_sze = repmat(20,1,numel(box_ind{iPL}));
    fcfg.aph_val = 0.45;
           
    fcfg.out_dir = [ out_dir '/' plt_nme{iPL} '/' ];
    fcfg.out_nme = [ 'box_plot' '_' mes_int_nme{iMES}];
    
    ejk_scatter(fcfg)
        
    % Diff plots
    if ~isempty(dff_ind{iPL})
        clear fcfg
        
        fcfg = [];
        
        fcfg.xdt = num2cell(1:numel(dff_ind{iPL}));
        for iD = 1:numel(dff_ind{iPL})
            fcfg.ydt{iD} = cell2mat(dta_hld.(mdl_int_nme{dff_ind{iPL}{iD}(1)}).(mes_int_nme{iMES})(:,3)) - ...
                cell2mat(dta_hld.(mdl_int_nme{dff_ind{iPL}{iD}(2)}).(mes_int_nme{iMES})(:,3));
            
            fcfg.fce_col{iD}     = mdl_int_col{dff_ind{iPL}{iD}(2)};
            fcfg.box_plt_col{iD} = mdl_int_col{dff_ind{iPL}{iD}(2)};
            fcfg.xlb(iD)         = mdl_int_nme(dff_ind{iPL}{iD}(2));
        end
        
        fcfg.edg_col = [ repmat({[0 0 0]},1,numel(dff_ind{iPL})) ];
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        
        fcfg.xlm = [ 0.5 numel(dff_ind{iPL})+5.5 ];
        fcfg.ylb = mes_int_nme(iMES);
        fcfg.mkr_sze = repmat(20,1,numel(dff_ind{iPL}));
        fcfg.aph_val = 0.45;
        
        fcfg.out_dir = [ out_dir '/' plt_nme{iPL} '/' ];
        fcfg.out_nme = [ 'diff_plot' '_' mes_int_nme{iMES}];
        
        fcfg.hln = 0;
        fcfg.hln_col = rgb('black');
        
        ejk_scatter(fcfg)
    end
    
  % Line plots
    figure('visible','off'); hold on;
    xdt_hld = [];
    for iD = 1:numel(box_ind{iPL})
        plot( cell2mat(dta_hld.(mdl_int_nme{box_ind{iPL}(iD)}).(mes_int_nme{iMES})(:,2)), ...
              cell2mat(dta_hld.(mdl_int_nme{box_ind{iPL}(iD)}).(mes_int_nme{iMES})(:,3)), ...
              'Color', mdl_int_col{box_ind{iPL}(iD)}, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', mdl_int_col{box_ind{iPL}(iD)});
        xdt_hld = [ xdt_hld ; cell2mat(dta_hld.(mdl_int_nme{box_ind{iPL}(iD)}).(mes_int_nme{iMES})(:,2))];
    end
    ylim(ylm_use{iPL}); xlim([0.5 max(xdt_hld)+0.5])
    print([ out_dir '/' plt_nme{iPL} '/'  'line_plot' '_' mes_int_nme{iMES} '.png'],'-dpng')
    close all
    
    end
end

%% Preliminary Tables
% Performance Table
avg_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
med_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
for iC = 1:numel(mes_int_nme)
    for iR = 1:numel(mdl_int_nme)
        med_val = nanmedian(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_val = nanmean(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        std_val = nanstd(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        iqr_val = iqr(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        max_val = nanmax(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        min_val = nanmin(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',avg_val,std_val, max_val);
        med_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',med_val,iqr_val, max_val);
    end
end
cell2csv([ out_dir '/' 'PerformanceTable_mean.csv'],   [ {''}  mes_int_nme; mdl_int_nme' avg_out_tbl ])
cell2csv([ out_dir '/' 'PerformanceTable_median.csv'], [ {''}  mes_int_nme; mdl_int_nme' med_out_tbl ])

% Diff table
cmp_grp = cat(1,dff_ind{cmp_tbl_ind});
avg_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));
med_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));
fdc_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));
cmp_nme     = cell(numel(cmp_grp),1);
for iC = 1:numel(mes_int_nme)
    for iDF = 1:numel(cmp_grp)
        cmp_nme{iDF,1} = [ mdl_int_nme{cmp_grp{iDF}(1)} ' - ' mdl_int_nme{cmp_grp{iDF}(2)}];
        
        med_val = nanmedian(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3))) - ...
                  nanmedian(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        avg_val = nanmean(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3))) - ...
                  nanmean(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        std_val = nanstd(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3)) - ...
                  cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        iqr_val = iqr(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3)) - ...
                  cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        max_val = nanmax(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3)) - ...
                  cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        min_val = nanmin(cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3)) - ...
                  cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3)));
        fdc_val = (sum( cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3)) > ...
                       cell2mat(dta_hld.(mdl_int_nme{cmp_grp{iDF}(2)}).(mes_int_nme{iC})(:,3))) / ...
                       numel(dta_hld.(mdl_int_nme{cmp_grp{iDF}(1)}).(mes_int_nme{iC})(:,3))) * 100;
              
        avg_out_tbl{iDF,iC} = sprintf('%.1f (%.1f); %.1f',avg_val,std_val, max_val);
        med_out_tbl{iDF,iC} = sprintf('%.1f (%.1f); %.1f',med_val,iqr_val, max_val);
        fdc_out_tbl{iDF,iC} = sprintf('%.1f',fdc_val);
    end
end
cell2csv([ out_dir '/' 'DiffTable_mean.csv'],   [ {''}  mes_int_nme; cmp_nme avg_out_tbl ])
cell2csv([ out_dir '/' 'DiffTable_median.csv'], [ {''}  mes_int_nme; cmp_nme med_out_tbl ])
cell2csv([ out_dir '/' 'DiffTable_FDIC.csv'],   [ {''}  mes_int_nme; cmp_nme fdc_out_tbl ])
