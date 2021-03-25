cmp_nme = { }


% POST-OPERATIVE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%
% Correlate  within LTLE
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Within_LTLE' '/'];

ejk_cross_cor( fcfg );
            
% Correlate within RTLE
iG = 4;

fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Within_RTLE' '/'];

ejk_cross_cor( fcfg );

% Correlate across LTLE/RTLE
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

iG = 3;
fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

iG = 4;
fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'PstOpr_Sample_Across_TLE' '/'];

ejk_cross_cor( fcfg );

% PRE-OPERATIVE DATASET %%%%%%%%%%%%%%%%%%%%%%%%%
% Correlate  within LTLE
iG = 1;

fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Within_LTLE' '/'];

ejk_cross_cor( fcfg );
            
% Correlate within RTLE
iG = 2;

fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Within_RTLE' '/'];

ejk_cross_cor( fcfg );

% Correlate across LTLE/RTLE
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.(nme_typ{iG}), 1);

iG = 1;
fcfg.dta_one = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_one = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.cor_typ = 'spearman';

iG = 2;
fcfg.dta_two = cell2mat(cog_dta( grp.(nme_typ{iG}), cog_col.(nme_typ{iG})));
fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Cognitive' '/' 'Total_Sample_Across_TLE' '/'];

ejk_cross_cor( fcfg );