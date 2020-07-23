
%% ejk_ttest2_independent
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:50]') );

fcfg.dta     = [ round(rand(50,1)*100) round(rand(50,1)*100) ];
fcfg.dta_nme = { 'Random1' 'Random2' };

fcfg.grp     = [ num2cell([ones(25,1) ; zeros(25,1)]) [ repmat({'ABC'},25,1) ; repmat({'XYZ'},25,1) ]];
fcfg.grp_nme = { 'Grouping1' 'Grouping2' };

fcfg.out_dir = '/home/ekaestner/Downloads/ttest_output';

ejk_ttest2_paired( fcfg );

%% ejk_ttest2_paired

%% non-parametric ttest

%% ejk_1way_anova
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:75]') );

fcfg.dta     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_nme = { 'Random1' 'Random2' };

fcfg.grp     = [ num2cell([ones(25,1)*2 ; ones(25,1) ; zeros(25,1)]) [ repmat({'ABC'},25,1) ; repmat({'XYZ'},25,1) ; repmat({'LMN'},25,1) ]];
fcfg.grp_nme = { 'Grouping1' 'Grouping2' };

fcfg.out_dir = '/home/ekaestner/Downloads/anova';

ejk_1way_anova( fcfg );

%% ejk_1way_ancova
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:75]') );

fcfg.dta     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_nme = { 'Random1' 'Random2' };

fcfg.grp     = [ num2cell([ones(25,1)*2 ; ones(25,1) ; zeros(25,1)]) [ repmat({'ABC'},25,1) ; repmat({'XYZ'},25,1) ; repmat({'LMN'},25,1) ]];
fcfg.grp_nme = { 'Grouping1' 'Grouping2' };

fcfg.cov     = [ num2cell( rand(75,1) ) repmat([ repmat({'Slow'},5,1) ; repmat({'Fast'},5,1) ; repmat({'Medium'},5,1) ],5,1) ];
fcfg.cov_nme = { 'Covariate1' 'Covariate2' };

fcfg.out_dir = '/home/ekaestner/Downloads/ancova';

ejk_1way_ancova( fcfg );

%% ejk_person_correlation
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:75]') );

fcfg.dta_one     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_one_nme = { 'Random1'             'Random2' };

fcfg.dta_two     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_two_nme = { 'Other1'              'Other2' };

fcfg.out_dir = '/home/ekaestner/Downloads/pearson';

ejk_pearson_correlation( fcfg );
