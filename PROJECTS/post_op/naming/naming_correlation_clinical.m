cmp_out = { 'TLE_Controls_pre_pre' ...
            'LTLE_Controls_pre_pre' ...
            'RTLE_Controls_pre_pre' ...
            'LTLE_post_post'  ...  
            'RTLE_post_post' };

cmp_nme = { 'tle_controls_pre_3T_allSurg_all' ... 
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_left'  ...  
            'tle_post_3T_ATLonly_right' };
        
cln_col_one     =   [ 4 5 6 7 8 ] ;
cln_col_tle     =   [ 6 7 8 ] ;
cln_col_tle_con =   [ 4 5 ] ;

cog_col_two = { [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 4 5 6 ] ...
                [ 4 5 6 ] };


%%
cln_dta = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ]);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);
            
%% Correlations
for iG = 1:numel(cmp_out)

    fcfg = [];
    
    fcfg.sbj_nme = cln_dta_sbj( grp.(cmp_nme{iG}), 1);
    
    fcfg.dta_one = cell2mat(cln_dta( grp.(cmp_nme{iG}), cln_col_one));
    fcfg.lbl_one = cln_dta_col(cln_col_one);
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_two = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_two{iG}));
    fcfg.lbl_two = strcat('x', cog_dta_col(cog_col_two{iG}));
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.20;
    
    fcfg.force_plot = 1;
    
    fcfg.out_dir = [ out_put '/' 'Clinical' '/' cmp_out{iG} '/' ];
    
    ejk_cross_cor( fcfg );
    
end

%% Comparisons
% RTLE / LTLE ttest
fcfg = [];

fcfg.sbj_nme = [ cln_dta_sbj( grp.tle_post_3T_ATLonly_left, 1) ; cln_dta_sbj( grp.tle_post_3T_ATLonly_right, 1) ];

fcfg.dta     = [ cell2mat(cln_dta( grp.tle_post_3T_ATLonly_left, cln_col_tle)) ; cell2mat(cln_dta( grp.tle_post_3T_ATLonly_right, cln_col_tle)) ];
fcfg.dta_nme = cln_dta_col(cln_col_tle);

fcfg.grp     = [ repmat({'left'},numel(grp.tle_post_3T_ATLonly_left),1) ; repmat({'right'},numel(grp.tle_post_3T_ATLonly_right),1) ];
fcfg.grp_nme = { 'TLE_post' };

fcfg.out_dir = [ out_put '/' 'Clinical' '/' 'Post_TLE_ttest' '/'];

ejk_ttest2_independent( fcfg );

% Pre-operative ANOVA
cfg.sbj_nme = [ cln_dta_sbj( grp.controls_pre_3T_allSurg_all, 1) ; cln_dta_sbj( grp.tle_pre_3T_allSurg_left, 1) ; cln_dta_sbj( grp.tle_pre_3T_allSurg_right, 1) ];

cfg.dta     = [ cell2mat(cln_dta( grp.controls_pre_3T_allSurg_all, cln_col_tle_con)) ; cell2mat(cln_dta( grp.tle_pre_3T_allSurg_left, cln_col_tle_con)) ; cell2mat(cln_dta( grp.tle_pre_3T_allSurg_right, cln_col_tle_con)) ];
cfg.dta_nme = cln_dta_col(cln_col_tle_con);

cfg.grp     = [ repmat({'controls'},numel(grp.controls_pre_3T_allSurg_all),1) ; repmat({'left'},numel(grp.tle_pre_3T_allSurg_left),1) ; repmat({'right'},numel(grp.tle_pre_3T_allSurg_right),1) ];
cfg.grp_nme = {'HC_TLE'};

cfg.out_dir = [ out_put '/' 'Clinical' '/' 'Pre_Omnibus_anova' '/'];

ejk_1way_anova( cfg )

