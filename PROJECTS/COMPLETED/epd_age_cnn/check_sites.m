clear; clc;    

%% Data
thk_dta = mmil_readtext( '/home/ekaestne/PROJECTS//OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' );
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

%% All Data
thk_dta = mmil_readtext( '/home/ekaestne/PROJECTS//OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_average_cleaned_no_wis.csv' );
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

%% Covariates
cov_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_nowisconsin_3T.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

%% ALL SITES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collate data
num_ind = 1;
for iS = 1:size(thk_sbj_nme, 1)
    
    cov_ind = find(strcmpi( cov_sbj_nme, thk_sbj_nme{iS,1} ));
    
    if ~isempty(cov_ind) && ~isempty(cov_dta{cov_ind, 5}) % 5 % 6 % 7
        
        thk_dta_sbj_use{num_ind,1} = thk_sbj_nme{iS,1};
        thk_dta_use(num_ind,:)     = thk_dta(iS,:);
        
        grp_var_use{num_ind,1}     = cov_dta{cov_ind, 5}; % 5 % 6 % 7
        
        cov_dta_use{num_ind,1}     = cov_dta{cov_ind, 1};
        if cov_dta{cov_ind, 2}==1
            cov_dta_use{num_ind,2}     = 'Male';
        elseif cov_dta{cov_ind, 2}==2
            cov_dta_use{num_ind,2}     = 'Female';
        end
        cov_dta_use{num_ind,3}     = cov_dta{cov_ind, 3};
        %if cov_dta{cov_ind, 4}==1
        %    cov_dta_use{num_ind,4}     = '1.5T';
        %elseif cov_dta{cov_ind, 4}==2
        %    cov_dta_use{num_ind,4}     = '3T';
        %end
        
        num_ind = num_ind+1;
        
    end
end

all_data = setxor( thk_dta_sbj_use, cov_sbj_nme(~cellfun(@isempty,cov_dta(:,5)),1) ); % 5 % 6 % 7

%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use;
fcfg.grp_nme = { 'Site' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/sitetest/ancova_no_wis_newdata'; % ancova_all_sites % ancova_no_wisc % ancova_3T

ejk_1way_ancova( fcfg );

%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = mean(thk_dta_use,2);
fcfg.dta_nme = {'tot'};

fcfg.grp     = grp_var_use;
fcfg.grp_nme = { 'Site' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/sitetest/ancova_all_sites_avg';  % ancova_all_sites_avg % ancova_no_wisc_avg % ancova_3T_avg

ejk_1way_ancova( fcfg );

%% 
clear cov_dta_use thk_dta_sbj_use thk_dta_use grp_dta_use

