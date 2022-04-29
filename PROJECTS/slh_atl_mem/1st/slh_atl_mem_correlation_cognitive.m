out_put = [ prj_dir '/' prj_nme '/' 'InitialAnalysis'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% Group Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Group Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_grp = {  { 'controls_pre_3T_allSurg_all'   'tle_pre_3T_allSurg_left'   'tle_pre_3T_allSurg_right'  } ...
            { 'controls_pre_allT_allSurg_all' 'tle_pre_allT_allSurg_left' 'tle_pre_allT_allSurg_right' }...
            { 'tle_pre_3T_allSurg_left'       'tle_pre_3T_allSurg_right' } ...
            { 'tle_pre_allT_allSurg_left'     'tle_pre_allT_allSurg_right' }...
            { 'tle_post_3T_ATLonly_left'      'tle_post_3T_ATLonly_right' } ...
            { 'tle_post_allT_ATLonly_left'    'tle_post_allT_ATLonly_right' } };
cmp_nme = { { 'HC_pre_3T'                     'LTLE_pre_3T'               'RTLE_pre_3T' } ...
            { 'HC_pre_allT'                   'LTLE_pre_allT'             'RTLE_pre_allT'} ...
            { 'LTLE_pre_3T'                   'RTLE_pre_3T' } ...
            { 'LTLE_pre_allT'                 'RTLE_pre_allT'  } ...
            { 'LTLE_post_3T'                  'RTLE_post_3T' } ...
            { 'LTLE_post_allT'                'RTLE_post_allT' } };
cmp_out = { 'TLE_Controls_pre_3T' ...
            'TLE_Controls_pre_allT' ...
            'TLE_Only_pre_3T' ...
            'TLE_Only_pre_allT' ...
            'TLE_Only_post_3T' ...
            'TLE_Only_post_allT' };   
        
cog_tst_lm2_pre = find(strcmpi(cog_dta_col,'log_mem_nor_scr_two'));
cog_tst_vpa_pre = find(strcmpi(cog_dta_col,'vp2_nor_scr'));
cog_tst_lm2_pst = find(strcmpi(cog_dta_col,'log_mem_nor_scr_two_pst'));
cog_tst_vpa_pst = find(strcmpi(cog_dta_col,'vp2_nor_scr_pst'));
        
cog_col_use = { [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre cog_tst_lm2_pst cog_tst_vpa_pst] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre cog_tst_lm2_pst cog_tst_vpa_pst]};
        
fsh_use = [ 1 2 3 4 5 6 ];
ttb_use = [ 3 4 5 6 ];
ttw_use = [ 1 2 5 6 ];
anv_use = [ 1 2 ];    

% Correlation Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_out_cor = { 'TLE_Controls_pre_cog_3T'  ...
                ...'LTLE_Controls_pre_cog_3T' ...
                ...'RTLE_Controls_pre_cog_3T' ...
                'LTLE_post_cog_3T'         ...  
                'RTLE_post_cog_3T'         ...
                ...'TLE_post_cog_3T' ...
                'TLE_Controls_pre_cog_allT'  ...
                ...'LTLE_Controls_pre_cog_allT' ...
                ...'RTLE_Controls_pre_cog_allT' ...
                'LTLE_post_cog_allT'         ...  
                'RTLE_post_cog_allT'         ...
                };... 'TLE_post_cog_allT' };

cmp_nme_cor = { 'tle_controls_pre_3T_allSurg_all' ... 
                ...'tle_controls_pre_3T_allSurg_left'  ...
                ...'tle_controls_pre_3T_allSurg_right' ...
                'tle_post_3T_ATLonly_left'  ...  
                'tle_post_3T_ATLonly_right' ...
                ...'tle_post_3T_ATLonly_all' ...
                'tle_controls_pre_allT_allSurg_all' ... 
                ...'tle_controls_pre_allT_allSurg_left'  ...
                ...'tle_controls_pre_allT_allSurg_right' ...
                'tle_post_allT_ATLonly_left'  ...  
                'tle_post_allT_ATLonly_right' ...
                };...'tle_post_allT_ATLonly_all' };
        
cog_col_one = { [ 1 2 ] ...
                [ 1 2 3 4 ] ...
                [ 1 2 3 4 ] ...
                [ 1 2 ] ...
                [ 1 2 3 4 ] ...
                [ 1 2 3 4 ] };

cog_col_two = cog_col_one;

%% Fishers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_typ_use = { { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' } ...
                { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' } ...
                { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' } ...
                { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' } ...
                { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' 'log_mem_nor_scr_two_pst_cat' 'vp2_nor_scr_pst_cat' } ...
                { 'log_mem_nor_scr_two_cat' 'vp2_nor_scr_cat' 'log_mem_nor_scr_two_pst_cat' 'vp2_nor_scr_pst_cat' } };

for iC = 1:numel(fsh_use)
    
    [~, use_dta_col ] = intersect( cog_dta_col, dta_typ_use{fsh_use(iC)} );
    
    fcfg = [];
    fcfg.grp     = grp;
    fcfg.grp_inc = cmp_grp(fsh_use(iC));
    fcfg.grp_nme = cmp_nme(fsh_use(iC));
    fcfg.dta = cog_dta(:,use_dta_col);
    fcfg.sbj = cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};
    
    % Caluclate Fishers
    fcfg = [];
    
    fcfg.sbj = grp_sbj{1};
    
    fcfg.dta_one = grp_dta{1};
    fcfg.lbl_one = cog_dta_col(use_dta_col);
    
    fcfg.dta_two = repmat(grp_typ{1},1,numel(use_dta_col));
    fcfg.lbl_two = strcat( 'group_', cog_dta_col(use_dta_col));
    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Fisher' '/' cmp_out{fsh_use(iC)}];
    
    ejk_fisher_test( fcfg );
    
end

%% t-tests independent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_typ_use = { { } ...
                { } ...
                { 'log_mem_nor_scr_two' 'vp2_nor_scr' } ...
                { 'log_mem_nor_scr_two' 'vp2_nor_scr' } ...
                { 'log_mem_nor_scr_two' 'vp2_nor_scr' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' } ...
                { 'log_mem_nor_scr_two' 'vp2_nor_scr' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' } };
 
for iC = 1:numel(ttb_use)
    
    [~, use_dta_col ] = intersect( cog_dta_col, dta_typ_use{ttb_use(iC)} );
    
    fcfg = [];
    fcfg.grp     = grp;
    fcfg.grp_inc = cmp_grp(ttb_use(iC));
    fcfg.grp_nme = cmp_nme(ttb_use(iC));
    fcfg.dta = cog_dta(:,use_dta_col);
    fcfg.sbj = cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    % Test
    fcfg = [];    
    fcfg.sbj_nme = grp_sbj{1};    
    fcfg.dta     = grp_dta{1};
    fcfg.dta_nme = cog_dta_col(:,use_dta_col);    
    fcfg.grp     = grp_typ{1};
    fcfg.grp_nme = cmp_out(1);    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ttest2' '/' cmp_out{ttb_use(iC)} ];    
    ejk_ttest2_independent( fcfg );
    
end

%% t-tests single %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_typ_use = { { 'log_mem_nor_scr_two_zsc' 'vp2_nor_scr_zsc' } ...
                { 'log_mem_nor_scr_two_zsc' 'vp2_nor_scr_zsc' } ...
                { } ...
                { } ...
                { 'log_mem_nor_scr_two_zsc' 'vp2_nor_scr_zsc' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' } ...
                { 'log_mem_nor_scr_two_zsc' 'vp2_nor_scr_zsc' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst' } };

for iC = 1:numel(ttw_use)
    for iG = 1:numel(cmp_nme{ttw_use(iC)})
        
    [~, use_dta_col ] = intersect( cog_dta_col, dta_typ_use{ttw_use(iC)} );
    
    fcfg = [];
    fcfg.grp     = grp;
    fcfg.grp_inc = {cmp_grp{ttw_use(iC)}(iG)};
    fcfg.grp_nme = {cmp_nme{ttw_use(iC)}(iG)};
    fcfg.dta = cog_dta(:,use_dta_col);
    fcfg.sbj = cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    
    fcfg = [];    
    fcfg.sbj_nme = grp_sbj{1};    
    fcfg.dta     = grp_dta{1};
    fcfg.dta_nme = cog_dta_col(:,use_dta_col);
    fcfg.men     = 0;    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ttest1' '/' cmp_nme{ttw_use(iC)}{iG} '/'];    
    ejk_ttest1( fcfg );
    
    end    
end

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_typ_use = { { 'log_mem_nor_scr_two' 'vp2_nor_scr' } ...
                { 'log_mem_nor_scr_two' 'vp2_nor_scr' } ...
                {} ...
                {} ...
                {} ...
                {}};

for iC = 1:numel(anv_use)
    
    [~, use_dta_col ] = intersect( cog_dta_col, dta_typ_use{anv_use(iC)});
    
    fcfg = [];
    fcfg.grp     = grp;
    fcfg.grp_inc = cmp_grp(anv_use(iC));
    fcfg.grp_nme = cmp_nme(anv_use(iC));
    fcfg.dta = cog_dta(:,use_dta_col);
    fcfg.sbj = cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    fcfg = [];
    fcfg.sbj_nme = grp_sbj{1};
    fcfg.dta     = grp_dta{1};
    fcfg.dta_nme = cog_dta_col(use_dta_col);
    fcfg.grp     = grp_typ{1};
    fcfg.grp_nme = cmp_out(1);
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ANOVA' '/'  cmp_out{anv_use(iC)}];
    ejk_1way_anova( fcfg )
    
end

%% Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out_cor)

    fcfg = [];
    
    fcfg.sbj_nme = cog_dta_sbj( grp.(cmp_nme_cor{iG}), 1);
    
    fcfg.dta_two = cell2mat(cog_dta( grp.(cmp_nme_cor{iG}), cog_col_one{iG}));
    fcfg.lbl_two = cog_dta_col(cog_col_one{iG});
    
    fcfg.cor_typ = 'pearson';
    
    fcfg.dta_one = cell2mat(cog_dta( grp.(cmp_nme_cor{iG}), cog_col_two{iG}));
    fcfg.lbl_one = cog_dta_col(cog_col_two{iG});
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.20;
    
    fcfg.force_plot = 1;
    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Correlation' '/' cmp_out_cor{iG} '/' ];
    
    ejk_cross_cor( fcfg );
    
end


