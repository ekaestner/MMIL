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
        
cog_col_one = { [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 4 5 6 ] ...
                [ 4 5 6 ] ...
                [ 1 2 3 ] ...
                [ 1 2 3 ] };

cog_col_two = { [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 1 2 3 ] ...
                [ 4 5 6 ] ...
                [ 4 5 6 ] ...
                [ 4 5 6 ] ...
                [ 4 5 6 ] };


%% Correlations
for iG = 1:numel(cmp_out)

    fcfg = [];
    
    fcfg.sbj_nme = cln_dta_sbj( grp.(cmp_nme{iG}), 1);
    
    fcfg.dta_one = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_one{iG}));
    fcfg.lbl_one = cog_dta_col(cog_col_one{iG});
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_two = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col_two{iG}));
    fcfg.lbl_two = strcat('x', cog_dta_col(cog_col_two{iG}));
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.20;
    
    fcfg.force_plot = 1;
    
    fcfg.out_dir = [ out_put '/' 'Cognitive' '/' cmp_out{iG} '/' ];
    
    ejk_cross_cor( fcfg );
    
end

%% Comparisons
% RTLE / LTLE ttest
fcfg = [];

fcfg.sbj_nme = [ cln_dta_sbj( grp.tle_post_3T_ATLonly_left, 1) ; cln_dta_sbj( grp.tle_post_3T_ATLonly_right, 1) ];

fcfg.dta     = [ cell2mat(cog_dta( grp.tle_post_3T_ATLonly_left, 4:6)) ; cell2mat(cog_dta( grp.tle_post_3T_ATLonly_right, 4:6)) ];
fcfg.dta_nme = cog_dta_col(4:6);

fcfg.grp     = [ repmat({'left'},numel(grp.tle_post_3T_ATLonly_left),1) ; repmat({'right'},numel(grp.tle_post_3T_ATLonly_right),1) ];
fcfg.grp_nme = { 'TLE_post' };

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Post_TLE_ttest' '/'];

ejk_ttest2_independent( fcfg );

% Pre-operative ANOVA
cfg.sbj_nme = [ cln_dta_sbj( grp.controls_pre_3T_allSurg_all, 1) ; cln_dta_sbj( grp.tle_pre_3T_allSurg_left, 1) ; cln_dta_sbj( grp.tle_pre_3T_allSurg_right, 1) ];

cfg.dta     = [ cell2mat(cog_dta( grp.controls_pre_3T_allSurg_all, 1:3)) ; cell2mat(cog_dta( grp.tle_pre_3T_allSurg_left, 1:3)) ; cell2mat(cog_dta( grp.tle_pre_3T_allSurg_right, 1:3)) ];
cfg.dta_nme = cog_dta_col(1:3);

cfg.grp     = [ repmat({'controls'},numel(grp.controls_pre_3T_allSurg_all),1) ; repmat({'left'},numel(grp.tle_pre_3T_allSurg_left),1) ; repmat({'right'},numel(grp.tle_pre_3T_allSurg_right),1) ];
cfg.grp_nme = {'HC_TLE'};

cfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Pre_Omnibus_anova' '/'];

ejk_1way_anova( cfg )

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % POST-OPERATIVE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%
% % Correlate  within LTLE
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Within_LTLE' '/'];
% 
% ejk_cross_cor( fcfg );
%             
% % Correlate within RTLE
% iG = 4;
% 
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Within_RTLE' '/'];
% 
% ejk_cross_cor( fcfg );
% 
% % Correlate across LTLE/RTLE
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% iG = 3;
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% iG = 4;
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Across_TLE' '/'];
% 
% ejk_cross_cor( fcfg );
% 
% % PRE-OPERATIVE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%
% % Correlate  within LTLE
% iG = 1;
% 
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Within_LTLE' '/'];
% 
% ejk_cross_cor( fcfg );
%             
% % Correlate within RTLE
% iG = 2;
% 
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Within_RTLE' '/'];
% 
% ejk_cross_cor( fcfg );
% 
% % Correlate across LTLE/RTLE
% fcfg = [];
% 
% fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);
% 
% iG = 1;
% fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.cor_typ = 'spearman';
% 
% iG = 2;
% fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
% fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
% 
% fcfg.pvl_cut = 0.05;
% fcfg.pvl_lib = 0.20;
% 
% fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Across_TLE' '/'];
% 
% ejk_cross_cor( fcfg );