out_put = [ prj_dir '/' prj_nme '/' 'SpecificCor'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% Group Comparisons
% Pre-operative ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
cmp_out = { 'TLE_Controls_pre_anova' };
cmp_grp = { { 'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right' } };
cmp_nme = { { 'HC'                          'LTLE'                    'RTLE' } };

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cog_dta;
fcfg.sbj = cog_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

% Test
for iG = 1:numel(cmp_out)
    cfg.sbj_nme = grp_sbj{iG};
    
    cfg.dta     = grp_dta{iG}(:,1:2);
    cfg.dta_nme = cog_dta_col(1:2);
    
    cfg.grp     = grp_typ{iG};
    cfg.grp_nme = cmp_out(iG);
    
    cfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ANOVA' '/' 'Pre_Omnibus_anova' '/'];
    
    ejk_1way_anova( cfg )
end

% Pre-operative deviation from 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
cmp_out = { 'TLE_post_ttest_left' 'TLE_post_ttest_right' };
cmp_grp = { { 'tle_post_3T_ATLonly_left' } { 'tle_post_3T_ATLonly_right' } };
cmp_nme = { { 'LTLE' } { 'RTLE' } };

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cog_dta;
fcfg.sbj = cog_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

for iG = 1:numel(grp_dta)
    fcfg = [];
    
    fcfg.sbj_nme = grp_sbj{iG};
    
    fcfg.dta     = grp_dta{iG};
    fcfg.dta_nme = cog_dta_col;
    
    fcfg.men     = 0;
    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ttest' '/' cmp_out{iG} '/'];
    
    ejk_ttest1( fcfg );
end

% Post-operative RTLE / LTLE ttest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
cmp_out = { 'TLE_post_ttest' };
cmp_grp = { { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' } };
cmp_nme = { { 'LTLE'                     'RTLE' } };

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cog_dta;
fcfg.sbj = cog_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

fcfg = [];

fcfg.sbj_nme = grp_sbj{1};

fcfg.dta     = grp_dta{1};
fcfg.dta_nme = cog_dta_col;

fcfg.grp     = grp_typ{1};
fcfg.grp_nme = cmp_out(1);

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'ttest' '/' 'Post_TLE_ttest' '/'];

ejk_ttest2_independent( fcfg );

%% Fischer's exact  test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_col_one = [ 3 4 ];
cog_col_two = [ 3 4 ];

grp_nme = { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

sbj_cat = [ cog_dta_sbj( grp.(grp_nme{1}), 1) ; cog_dta_sbj( grp.(grp_nme{2}), 1)];

for iT = 1:numel(cog_col_one)
    
    % Calculate impairments
    dta_one( cell2mat(cog_dta( grp.(grp_nme{1}), cog_col_one(iT) )) <= -1.5, 1) = {'Impaired'};
    dta_one( cell2mat(cog_dta( grp.(grp_nme{1}), cog_col_one(iT) )) >  -1.5, 1) = {'NotImpaired'};
    dta_one( isnan(cell2mat(cog_dta( grp.(grp_nme{1}), cog_col_one(iT) ))), 1) = {NaN};
    
    dta_two( cell2mat(cog_dta( grp.(grp_nme{2}), cog_col_two(iT) )) <= -1.5, 1) = {'Impaired'};
    dta_two( cell2mat(cog_dta( grp.(grp_nme{2}), cog_col_two(iT) )) >  -1.5, 1) = {'NotImpaired'};
    dta_two( isnan(cell2mat(cog_dta( grp.(grp_nme{2}), cog_col_two(iT) ))), 1)  = {NaN};
    
    grp_cat(:,iT) = [ dta_one ; dta_two ];
    
    % Make groups Groups
    typ_one = repmat( grp_nme(1),numel(grp.(grp_nme{1})),1 );
    typ_two = repmat( grp_nme(2),numel(grp.(grp_nme{2})),1 );
    
    typ_cat(:,iT) = [ typ_one ; typ_two ];
    
end

% Caluclate Fishers
fcfg = [];
            
fcfg.sbj = sbj_cat;

fcfg.dta_one = grp_cat;
fcfg.lbl_one = cog_dta_col(cog_col_one);

fcfg.dta_two = typ_cat;
fcfg.lbl_two = strcat( 'group_', cog_dta_col(cog_col_one));

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Fisher' '/'];

ejk_fisher_test( fcfg );

%% Correlations
% Correlation setups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_out = { 'TLE_Controls_pre_pre' ...
            'LTLE_Controls_pre_pre' ...
            'RTLE_Controls_pre_pre' ...
            'LTLE_post_post'  ...  
            'RTLE_post_post' ...
            'LTLE_pre_post'  ...  
            'RTLE_pre_post' };

cmp_nme = { 'tle_controls_pre_3T_allSurg_all' ... 
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_left'  ...  
            'tle_post_3T_ATLonly_right' ...
            'tle_post_3T_ATLonly_left'  ...  
            'tle_post_3T_ATLonly_right' };
        
cog_col_one = { [ 1 2 ] ...
                [ 1 2 ] ...
                [ 1 2 ] ...
                [ 3 4 ] ...
                [ 3 4 ] ...
                [ 1 2 ] ...
                [ 1 2 ] };

cog_col_two = { [ 1 2 ] ...
                [ 1 2 ] ...
                [ 1 2 ] ...
                [ 3 4 ] ...
                [ 3 4 ] ...
                [ 3 4 ] ...
                [ 3 4 ] };
            
% Run correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out)

    fcfg = [];
    
    fcfg.sbj_nme = cog_dta_sbj( grp.(cmp_nme{iG}), 1);
    
    fcfg.dta_one = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_one{iG}));
    fcfg.lbl_one = cog_dta_col(cog_col_one{iG});
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_two = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_two{iG}));
    fcfg.lbl_two = strcat('x', cog_dta_col(cog_col_two{iG}));
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.20;
    
    fcfg.force_plot = 1;
    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Correlation' '/' cmp_out{iG} '/' ];
    
    ejk_cross_cor( fcfg );
    
end

