
%% Load Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_sub_grp_fle ];
[ dta_sub_grp_dta, dta_sub_grp_sbj, dta_sub_grp_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_fll_grp_rev_fle ];
[ dta_fll_grp_dta, dta_fll_grp_sbj, dta_fll_grp_col ] = ejk_dta_frm(fcfg);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organize Full data
mdl_int     = { '3D FCNet (non-harmonized data)' '3D FCNet (harmonized data)' '2D FCNet (non-harmonized data)' '2D FCNet (harmonized data)' 'SVM (hippocampus left and right)' 'SVM (all features)' };
mdl_int_nme = { 'fcn_3dm_org'                    'fcn_3dm_cbt'                'fcn_2dm_org'                    'fcn_2dm_cbt'                'svm_hip'                          'svm_all'};
mdl_int_col = { rgb('bright teal')               rgb('light green')           rgb('maroon')                    rgb('light maroon')          rgb('dark grey')                   rgb('light grey')};

mes_int     = { 'Accuracy' 'Sensitivity' 'Specificity' 'PPV' 'NPV' };
mes_int_nme = { 'Accuracy' 'Specificity' 'Sensitivity' 'PPV' 'NPV' };

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
mes_int_nme = [ mes_int 'Precision' 'Recall' 'F1' ];
mdl_int_nme = fieldnames(dta_hld)';
mdl_int_col = mdl_int_col;

% Save
save([ dta_dir '/' 'performance_flip_reihaneh_revision.mat' ],'dta_hld','mes_int_nme', 'mdl_int_nme', 'mdl_int_col');
