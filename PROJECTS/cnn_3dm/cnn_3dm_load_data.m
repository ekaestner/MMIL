
%% Load Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_sub_grp_fle ];
[ dta_sub_grp_dta, dta_sub_grp_sbj, dta_sub_grp_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_fll_grp_fle ];
[ dta_fll_grp_dta, dta_fll_grp_sbj, dta_fll_grp_col ] = ejk_dta_frm(fcfg);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organize Subset Data
mdl_int_sub         = { '3D FCNet (full data)' '3D FCNet (750 samples)' '3D FCNet (300 samples)' '2D FCNet (full data)' '2D FCNet (750 samples)' '2D FCNet (300 samples)'   };
mdl_int_sub_nme     = { 'fcn_3dm_fll'          'fcn_3dm_750'            'fcn_3dm_300'            'fcn_2dm_fll'          'fcn_2dm_750'            'fcn_2dm_300'            };
mdl_int_sub_col = { rgb('bright teal')     rgb('teal')              rgb('dark teal')         rgb('light maroon')   rgb('maroon')            rgb('dark maroon')      };     
 
mes_int_sub     = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };
mes_int_sub_nme = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };

% Add model number
num_mdl = zeros(1,numel(mdl_int_sub)); 
add_mdl_num = [];
for iMDL = 1:numel(mdl_int_sub) 
    num_mdl(iMDL) = sum(strcmpi(dta_sub_grp_dta(:,1),mes_int_sub{1}) & strcmpi(dta_sub_grp_dta(:,2),mdl_int_sub{iMDL})); 
    ttt = repmat(1:num_mdl(iMDL),numel(mes_int_sub_nme),1);
    add_mdl_num = [ add_mdl_num ; ttt(:)];
end

dta_sub_grp_dta     = [ dta_sub_grp_dta num2cell(add_mdl_num)];
dta_sub_grp_col = [ dta_sub_grp_col 'ModelRun' ];

% Take Data
for iMDL = 1:numel(mdl_int_sub)
    for iMES = 1:numel(mes_int_sub)
        
        dta_ind = strcmpi(dta_sub_grp_dta(:,strcmpi(dta_sub_grp_col,'Model')),mdl_int_sub{iMDL}) & ...
                  strcmpi(dta_sub_grp_dta(:,strcmpi(dta_sub_grp_col,'Metric')),mes_int_sub{iMES}); 
        
        dta_hld.(mdl_int_sub_nme{iMDL}).(mes_int_sub_nme{iMES}) = dta_sub_grp_dta( dta_ind, [1 4 3] );
        
    end
    
    % Add Precision
    dta_hld.(mdl_int_sub_nme{iMDL}).Precision = dta_hld.(mdl_int_sub_nme{iMDL}).PPV;
     
    % Add Recall
    dta_hld.(mdl_int_sub_nme{iMDL}).Recall = dta_hld.(mdl_int_sub_nme{iMDL}).Sensitivity;
    
    % Add F1
    prc     = cell2mat(dta_hld.(mdl_int_sub_nme{iMDL}).Precision(:,3));
    rcl     = cell2mat(dta_hld.(mdl_int_sub_nme{iMDL}).Recall(:,3));
    row_sze = size(dta_hld.(mdl_int_sub_nme{iMDL}).Recall(:,3),1);
    dta_hld.(mdl_int_sub_nme{iMDL}).F1 = [ repmat({'F1'},row_sze,1) num2cell([1:row_sze]') num2cell(2 * ((prc.*rcl)./(prc+rcl))) ];
    
end

mes_int_sub     = [ mes_int_sub 'Precision' 'Recall' 'F1' ];
mes_int_sub_nme = mes_int_sub ;

%% Organize Full data
mdl_int     = { '3D FCNet (non-harmonized data)' '3D FCNet (harmonized data)' '2D FCNet (non-harmonized data)' '3D FCNet-shuffled labels (non-harmonized data)' };
mdl_int_nme = { 'fcn_3dm_org'                    'fcn_3dm_cbt'            'fcn_2dm_org'                    'fcn_3dm_shf' };
mdl_int_col = { rgb('bright teal')               rgb('light green')           rgb('maroon')                    rgb('grey') };

mes_int     = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };
mes_int_nme = { 'Accuracy' 'AUC' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };

% Add model number
num_mdl = zeros(1,numel(mdl_int)); 
add_mdl_num = [];
for iMDL = 1:numel(mdl_int) 
    num_mdl(iMDL) = sum(strcmpi(dta_fll_grp_dta(:,1),mes_int{1}) & strcmpi(dta_fll_grp_dta(:,2),mdl_int{iMDL})); 
    ttt = repmat(1:num_mdl(iMDL),numel(mes_int_nme),1);
    add_mdl_num = [ add_mdl_num ; ttt(:)];
end

dta_fll_grp_dta     = [ dta_fll_grp_dta num2cell(add_mdl_num)];
dta_fll_grp_col = [ dta_fll_grp_col 'ModelRun' ];

% Take Data
for iMDL = 1:numel(mdl_int)
    for iMES = 1:numel(mes_int)
        
        dta_ind = strcmpi(dta_fll_grp_dta(:,strcmpi(dta_fll_grp_col,'Model')),mdl_int{iMDL}) & ...
                  strcmpi(dta_fll_grp_dta(:,strcmpi(dta_fll_grp_col,'Metric')),mes_int{iMES}); 
        
        dta_hld.(mdl_int_nme{iMDL}).(mes_int_nme{iMES}) = dta_fll_grp_dta( dta_ind, [1 4 3] );
        
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



%% Save
% Combine
mes_int_nme = mes_int_sub_nme;
mdl_int_nme = fieldnames(dta_hld)';
mdl_int_col = [ mdl_int_sub_col mdl_int_col ];

% Save
save([ dta_dir '/' 'performance.mat' ],'dta_hld','mes_int_nme', 'mdl_int_nme', 'mdl_int_col');
