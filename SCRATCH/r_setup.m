
%% ejk_ttest2_independent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:50]') );

fcfg.dta     = [ round(rand(50,1)*100) round(rand(50,1)*100) ];
fcfg.dta_nme = { 'Random1' 'Random2' };

fcfg.grp     = [ num2cell([ones(25,1) ; zeros(25,1)]) [ repmat({'ABC'},25,1) ; repmat({'XYZ'},25,1) ]];
fcfg.grp_nme = { 'Grouping1' 'Grouping2' };

fcfg.out_dir = '/home/ekaestner/Downloads/ttest_output';

ejk_ttest2_paired( fcfg );

%% ejk_ttest2_paired %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% non-parametric ttest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ejk_1way_anova %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:75]') );

fcfg.dta     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_nme = { 'Random1' 'Random2' };

fcfg.grp     = [ num2cell([ones(25,1)*2 ; ones(25,1) ; zeros(25,1)]) [ repmat({'ABC'},25,1) ; repmat({'XYZ'},25,1) ; repmat({'LMN'},25,1) ]];
fcfg.grp_nme = { 'Grouping1' 'Grouping2' };

fcfg.out_dir = '/home/ekaestner/Downloads/anova';

ejk_1way_anova( fcfg );

%% ejk_1way_ancova %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% ejk_surface_1way_ancova %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_chs = .05;
pvl_cls = .05;
smt_stp = 176;

%
dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA';

new_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/epd_age_covariates_no_nan_no_mass.csv');

grp_dta     = new_dta(:,[1 6]);
grp_sbj_nme = grp_dta(2:end,1);
grp_roi_nme = grp_dta(1,2:end);
grp_dta     = grp_dta(2:end,2:end);

cov_dta     = new_dta(:,[1 2 3 4 8]);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%
fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = grp_sbj_nme;

fcfg.dta_lhs = [ dta_dir '/' 'surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_aMRI_thickness_rhs_sm176_no_nan_no_mass.mat'];

fcfg.grp     = grp_dta;
fcfg.grp_cmp = { { { 'TLE' 'HC' } { 'MCI' 'HC' } {'TLE' 'MCI'} } };
fcfg.grp_nme = grp_roi_nme;

fcfg.cov     = cov_dta;
fcfg.cov_nme = cov_roi_nme;

fcfg.out_dir = '/home/ekaestner/Downloads/surface_ancova';

ejk_surface_1way_ancova( fcfg );

%% ejk_surface_correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_chs = .20;
pvl_cls = .20;
smt_stp = 176;

%
dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA';

new_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_aparc__Epilepsy_and_Aging_updated_withTransverse.csv');
org_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/epd_age_covariates_no_nan_no_mass.csv');

cor_dta     = new_dta(:,[1 7 24]);
cor_sbj_nme = cor_dta(2:end,1);
cor_roi_nme = cor_dta(1,2:end);
cor_dta     = cor_dta(2:end,2:end);

grp_dta     = org_dta(:,[1 6]);
grp_sbj_nme = grp_dta(2:end,1);
grp_roi_nme = grp_dta(1,2:end);
grp_dta     = grp_dta(2:end,2:end);

cov_dta     = org_dta(:,[1 2 3 4 8]);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%
fcfg = [];

fcfg.smt_stp = smt_stp;
fcfg.pvl_chs = pvl_chs;
fcfg.pvl_cls = pvl_cls;

fcfg.sbj_nme = cor_sbj_nme;

fcfg.dta_lhs = [ dta_dir '/' 'surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat'];
fcfg.dta_rhs = [ dta_dir '/' 'surf_aMRI_thickness_rhs_sm176_no_nan_no_mass.mat'];

fcfg.cor     = cor_dta;
fcfg.cor_nme = cor_roi_nme;

fcfg.grp     = grp_dta;
fcfg.grp_nme = grp_roi_nme;
fcfg.grp_cmp = { 'TLE' 'HC' 'MCI' };

fcfg.cov     = cov_dta;
fcfg.cov_nme = cov_roi_nme;

fcfg.out_dir = '/home/ekaestner/Downloads/surface_correlations';

ejk_surface_correlations( fcfg );


%% ejk_pearson_correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = strcat('s', num2str([1:75]') );

fcfg.dta_one     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_one_nme = { 'Random1'             'Random2' };

fcfg.dta_two     = [ round(rand(75,1)*100) round(rand(75,1)*100) ];
fcfg.dta_two_nme = { 'Other1'              'Other2' };

fcfg.out_dir = '/home/ekaestner/Downloads/pearson';

ejk_pearson_correlation( fcfg );

%% Residual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data %%%%
thk_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

cov_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_updated.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

% Collate Data %%%%
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

% Residual Matrix %%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.cov     = cov_dta_use(:,1);
fcfg.cov_nme = cov_roi_nme(1);

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/residuals'; % split % desikan

rsd_mtx = ejk_residual_matrix( fcfg );

scatter()

%% Cross-Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lod_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    lod_dta( cellfun(@isempty,lod_dta) ) = {NaN};

neu_bio_dta     = lod_dta(:,[1 14:36]);
neu_bio_sbj_nme = neu_bio_dta(2:end,1);
neu_bio_roi_nme = neu_bio_dta(1,2:end);
neu_bio_dta     = cell2mat(neu_bio_dta(2:end,2:end));

cog_scr_dta     = lod_dta(:,[1 5:13]);
cog_scr_sbj_nme = cog_scr_dta(2:end,1);
cog_scr_roi_nme = cog_scr_dta(1,2:end);
cog_scr_dta     = cell2mat(cog_scr_dta(2:end,2:end));

cov_dta     = lod_dta(:,1:4);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%
fcfg = [];

fcfg.sbj_nme = neu_bio_sbj_nme;

fcfg.dta_one = neu_bio_dta;
fcfg.lbl_one = neu_bio_roi_nme;

fcfg.dta_two = cog_scr_dta;
fcfg.lbl_two = cog_scr_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.15;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/scores_by_wmparc';

ejk_cross_cor( fcfg );






