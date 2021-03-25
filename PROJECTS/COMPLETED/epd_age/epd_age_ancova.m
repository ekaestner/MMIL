clear; clc;    

%% Data
thk_dta = mmil_readtext( '/home/ekaestne/PROJECTS//OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_aparc__Epilepsy_and_Aging_updated_withTransverse.csv' ); %MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv %
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

%% 
grp_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/sbj_nme_sm313_no_nan.csv');
grp_sbj_nme = grp_dta(2:end,1);
grp_roi_nme = grp_dta(1,2:end);
grp_dta     = grp_dta(2:end,2:end);

%% Covariates
cov_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_updated.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%% Collate data
num_ind = 1;
for iS = 1:size(thk_sbj_nme, 1)
    
    cov_ind = find(strcmpi( cov_sbj_nme, thk_sbj_nme{iS,1} ));
    grp_ind = find(strcmpi( grp_sbj_nme, thk_sbj_nme{iS,1} ));
    
    if ~isempty(cov_ind) %&& ~isempty(cov_dta{cov_ind, 5}) % 5 % 6 % 7
        
        thk_dta_sbj_use{num_ind,1} = thk_sbj_nme{iS,1};
        thk_dta_use(num_ind,:)     = thk_dta(iS,:);
        
        grp_var_use{num_ind,1}     = grp_dta{grp_ind, 1}; % 5 % 6 % 7
        
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
        
        num_ind = num_ind+1;
        
    end
end

%% ANCOVA
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use;
fcfg.grp_nme = { 'Diagnosis' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan'; % split % desikan

ejk_1way_ancova( fcfg );







