out_put = [ prj_dir '/' prj_nme '/' 'SpecificCor'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% Pre-Surgical Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_grp = { { 'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right' } };
cmp_nme = { { 'HC'                          'LTLE'                    'RTLE' } };

%% ANOVA
cmp_out = { 'TLE_Controls_pre_anova' };

% Age, Education, WTAR
[~, use_dta_col ] = intersect( cln_dta_col, { 'AgeAtSurgery' 'Educ' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cln_dta(:,use_dta_col);
fcfg.sbj = cln_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

% Test
for iG = 1:numel(cmp_out)
    cfg.sbj_nme = grp_sbj{iG};
    
    cfg.dta     = grp_dta{iG};
    cfg.dta_nme = cln_dta_col(use_dta_col);
    
    cfg.grp     = grp_typ{iG};
    cfg.grp_nme = cmp_out(iG);
    
    cfg.out_dir = [ out_put '/' 'Clinical' '/' 'ANOVA' '/' 'Pre_Omnibus_anova' '/'];
    
    ejk_1way_anova( cfg )
end

%% Fisher's test
cmp_out = { 'TLE_Controls_pre_fishers' };

% Sex, Handedness, Race, Ethnicity
[~, use_dta_col ] = intersect( cln_dta_col, { 'Sex' 'Handedness' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cln_dta(:,use_dta_col);
fcfg.sbj = cln_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};

% Caluclate Fishers
fcfg = [];
            
fcfg.sbj = grp_sbj{1};

fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cln_dta_col(use_dta_col);

fcfg.dta_two = repmat(grp_typ{1},1,numel(use_dta_col));
fcfg.lbl_two = strcat( 'group_', cln_dta_col(use_dta_col));

fcfg.out_dir = [ out_put '/' 'Clinical' '/' 'Fisher' '/' cmp_out{1}];

ejk_fisher_test( fcfg );

%% Post-Surgical Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_grp = { { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' } };
cmp_nme = { { 'LTLE'                     'RTLE' } };

%% ttest
% Setup
cmp_out = { 'TLE_post_ttest' };

% Age, Education, WTAR, Age Onset, Duration, # ASMs, Seizure Frequency
[~, use_dta_col ] = intersect( cln_dta_col, { 'AgeAtSurgery' 'Educ' 'AgeOfSeizureOnset' 'NumAEDs' 'SeizureFreq' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cln_dta(:,use_dta_col);
fcfg.sbj = cln_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

% Test
fcfg = [];

fcfg.sbj_nme = grp_sbj{1};

fcfg.dta     = grp_dta{1};
fcfg.dta_nme = cln_dta_col(:,use_dta_col);

fcfg.grp     = grp_typ{1};
fcfg.grp_nme = cmp_out(1);

fcfg.out_dir = [ out_put '/' 'Clinical' '/' 'ttest' '/' 'Post_TLE_ttest' '/'];

ejk_ttest2_independent( fcfg );

%% Fisher's test
cmp_out = { 'TLE_post_fishers' };

% Sex, Handedness, Race, Ethnicity, MTS, Engel
[~, use_dta_col ] = intersect( cln_dta_col, { 'Sex' 'Handedness' 'MTS' 'EngelOutcome' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = cmp_grp;
fcfg.grp_nme = cmp_nme;
fcfg.dta = cln_dta(:,use_dta_col);
fcfg.sbj = cln_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

grp_dta{1}( ~strcmpi(grp_dta{1}(:,1),"I") & ~cellfun(@isempty,grp_dta{1}(:,1)), 1 ) = {'II+'};

grp_dta{1}( strcmpi(grp_dta{1}(:,3),'L') ,3) = {'Y'};
grp_dta{1}( strcmpi(grp_dta{1}(:,3),'R') ,3) = {'Y'};
grp_dta{1}( strcmpi(grp_dta{1}(:,3),'N/A') ,3) = {'N'};
grp_dta{1}( cellfun(@isempty,grp_dta{1}(:,3)) ,3) = {'N'};

grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};

% Caluclate Fishers
fcfg = [];
            
fcfg.sbj = grp_sbj{1};

fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cln_dta_col(use_dta_col);

fcfg.dta_two = repmat(grp_typ{1},1,numel(use_dta_col));
fcfg.lbl_two = strcat( 'group_', cln_dta_col(use_dta_col));

fcfg.out_dir = [ out_put '/' 'Clinical' '/' 'Fisher' '/' cmp_out{1}];

ejk_fisher_test( fcfg );

%% Correlations
% Correlation setups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_out = { 'TLE_Controls_pre_cln' ...
            'LTLE_Controls_pre_cln' ...
            'RTLE_Controls_pre_cln' ...
            'LTLE_post_cln'  ...  
            'RTLE_post_cln' ...
            'TLE_post_cln' };

cmp_nme = { 'tle_controls_pre_3T_allSurg_all' ... 
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_left'  ...  
            'tle_post_3T_ATLonly_right' ...
            'tle_post_3T_ATLonly_all' };
        
cog_col_one = { [ 1 2 ] ...
                [ 1 2 ] ...
                [ 1 2 ] ...
                [ 3 4 ] ...
                [ 3 4 ] ...
                [ 3 4 ] };

cln_col_two = { [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] ...
                [ 4 5 6 7 8 ] };
            
% Run correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out)

    fcfg = [];
    
    fcfg.sbj_nme = cln_dta_sbj( grp.(cmp_nme{iG}), 1);
    
    fcfg.dta_two = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_one{iG}));
    fcfg.lbl_two = cog_dta_col(cog_col_one{iG});
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_one = cell2mat(cln_dta( grp.(cmp_nme{iG}), cln_col_two{iG}));
    fcfg.lbl_one = cln_dta_col(cln_col_two{iG});
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.20;
    
    fcfg.force_plot = 1;
    
    fcfg.out_dir = [ out_put '/' 'Clinical' '/' 'Correlation' '/' cmp_out{iG} '/' ];
    
    ejk_cross_cor( fcfg );
    
end


