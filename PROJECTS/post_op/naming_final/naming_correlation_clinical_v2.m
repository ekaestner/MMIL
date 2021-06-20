out_put = [ prj_dir '/' prj_nme '/' 'SpecificCor'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

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
