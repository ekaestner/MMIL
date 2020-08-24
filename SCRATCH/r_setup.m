
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

%% Residual 
% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

cov_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_updated.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

% Collate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_ind = 1;
for iS = 1:size(thk_sbj_nme, 1)
    
    cov_ind = find(strcmpi( cov_sbj_nme, thk_sbj_nme{iS,1} ));
    
    if ~isempty(cov_ind) %&& ~isempty(cov_dta{cov_ind, 5}) % 5 % 6 % 7
        
        thk_dta_sbj_use{num_ind,1} = thk_sbj_nme{iS,1};
        thk_dta_use(num_ind,:)     = thk_dta(iS,:);
               
        cov_dta_use{num_ind,1}     = cov_dta{cov_ind, 1};
        if cov_dta{cov_ind, 2}==1
            cov_dta_use{num_ind,2}     = 'Male';
        elseif cov_dta{cov_ind, 2}==2
            cov_dta_use{num_ind,2}     = 'Female';
        end
        cov_dta_use{num_ind,3}     = cov_dta{cov_ind, 3};
        if cov_dta{cov_ind, 4}==1
           cov_dta_use{num_ind,4}     = '1.5T';
        elseif cov_dta{cov_ind, 4}==2
           cov_dta_use{num_ind,4}     = '3T';
        end

        ste_dta_use{num_ind,1} = cov_dta{cov_ind, 5};
        
        num_ind = num_ind+1;
       
    else
        error('')
    end
end

% Residual Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.cov     = cov_dta_use(:,1);
fcfg.cov_nme = cov_roi_nme(1);

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/residuals'; % split % desikan

rsd_mtx = ejk_residual_matrix( fcfg );

scatter()














