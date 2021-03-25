clear; clc;

%% Load Data
thk_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_aparc__Epilepsy_and_Aging_updated_withTransverse.csv');
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

vol_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/ROI/MRI_volume_updated.csv');
vol_sbj_nme = vol_dta(2:end,1);
vol_roi_nme = strcat('x', cellfun(@(x) mmil_spec_char(x,{'-'}), vol_dta(1,2:end), 'uni', 0) );
vol_dta     = cell2mat(vol_dta(2:end,2:end));

cov_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_updated.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);
cov_dta{34,3} = nan;

grp_dta = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/projects/Epilepsy_and_Aging/epilepsy_aging_graph_sample.csv'); % mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/sbj_nme_sm313_no_nan.csv');
grp_sbj_nme = grp_dta(2:end,1);
grp_roi_nme = grp_dta(1,2:end);
grp_dta     = grp_dta(2:end,2:end);

%% Collate data
num_ind = 1;
for iS = 1:size(thk_sbj_nme, 1)
    
    vol_ind = find(strcmpi( vol_sbj_nme, thk_sbj_nme{iS,1} ));
    cov_ind = find(strcmpi( cov_sbj_nme, thk_sbj_nme{iS,1} ));
    grp_ind = find(strcmpi( grp_sbj_nme, thk_sbj_nme{iS,1} ));
    
    if ~isempty(cov_ind) %&& ~isempty(cov_dta{cov_ind, 5}) % 5 % 6 % 7
        
        thk_dta_sbj_use{num_ind,1} = thk_sbj_nme{iS,1};
        thk_dta_use(num_ind,:)     = thk_dta(iS,:);
        
        vol_dta_sbj_use{num_ind,1} = vol_sbj_nme{vol_ind,1};
        vol_dta_use(num_ind,:)     = vol_dta(vol_ind,:);
        
        grp_var_use{num_ind,1}     = grp_dta{grp_ind, 1}; % 5 % 6 % 7
        grp_var_use{num_ind,2}     = grp_dta{grp_ind, 2}; % 5 % 6 % 7
        
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

%% ComBat & Recombine
epd_ind = find(strcmpi( grp_var_use, 'EPD_Old' ));

epd_sbj_nme = thk_dta_sbj_use( epd_ind );
epd_thk_dta = thk_dta_use( epd_ind, : );
epd_vol_dta = vol_dta_use( epd_ind, : );
btc_dta     = ste_dta_use( epd_ind );
cov_dta     = cov_dta_use( epd_ind, 4);
% 
% % THICKNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcfg = [];
% 
% fcfg.sbj_nme = epd_sbj_nme;
% 
% fcfg.dta     = epd_thk_dta;
% fcfg.dta_nme = thk_roi_nme';
% 
% fcfg.btc     = btc_dta;
% fcfg.btc_nme = {'Site'};
% 
% %fcfg.cov     = cov_dta;
% %fcfg.cov_nme = cov_roi_nme(4);
% 
% fcfg.out_dir = '/home/ekaestner/Downloads/ComBat';
% 
% com_bat_epd_dta = ejk_ComBat(fcfg);
% 
% thk_dta_com_bat            = thk_dta_use;
% thk_dta_com_bat(epd_ind,:) = com_bat_epd_dta;
% 
% % VOLUME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcfg = [];
% 
% fcfg.sbj_nme = epd_sbj_nme;
% 
% fcfg.dta     = epd_vol_dta;
% fcfg.dta_nme = vol_roi_nme';
% 
% fcfg.btc     = btc_dta;
% fcfg.btc_nme = {'Site'};
% 
% %fcfg.cov     = cov_dta;
% %fcfg.cov_nme = cov_roi_nme(4);
% 
% fcfg.out_dir = '/home/ekaestner/Downloads/ComBat';
% 
% com_bat_epd_dta = ejk_ComBat(fcfg);
% 
% vol_dta_com_bat            = vol_dta_use;
% vol_dta_com_bat(epd_ind,:) = com_bat_epd_dta;
% for iC = 1:size(vol_dta_com_bat,2)-1
%     vol_dta_use(:,iC)     = vol_dta_use(:,iC) ./ vol_dta_use(:,size(vol_dta_use,2));
%     vol_dta_com_bat(:,iC) = vol_dta_com_bat(:,iC) ./ vol_dta_com_bat(:,size(vol_dta_com_bat,2));
% end
% vol_roi_nme = vol_roi_nme(:, 1:size(vol_roi_nme,2)-1);
% vol_dta_use = vol_dta_use(:, 1:size(vol_dta_use,2)-1);
% vol_dta_com_bat = vol_dta_com_bat(:, 1:size(vol_dta_com_bat,2)-1);
% 
% vol_roi_nme = vol_roi_nme(:,[1:32 34:end]);
% vol_dta_use = vol_dta_use(:,[1:32 34:end]);
% vol_dta_com_bat = vol_dta_com_bat(:,[1:32 34:end]);

% Combined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = epd_sbj_nme;

fcfg.dta     = [ epd_thk_dta  epd_vol_dta ];
fcfg.dta_nme = [ thk_roi_nme' ; vol_roi_nme'];

fcfg.btc     = btc_dta;
fcfg.btc_nme = {'Site'};

%fcfg.cov     = cov_dta;
%fcfg.cov_nme = cov_roi_nme(4);

fcfg.out_dir = '/home/ekaestner/Downloads/ComBat';

com_bat_epd_dta = ejk_ComBat(fcfg);

thk_dta_com_bat            = thk_dta_use;
thk_dta_com_bat(epd_ind,:) = com_bat_epd_dta(:,1:size(thk_dta_use,2));

vol_dta_com_bat            = vol_dta_use;
vol_dta_com_bat(epd_ind,:) = com_bat_epd_dta(:,size(thk_dta_use,2)+1:end);
for iC = 1:size(vol_dta_com_bat,2)-1
    vol_dta_use(:,iC)     = vol_dta_use(:,iC) ./ vol_dta_use(:,size(vol_dta_use,2));
    vol_dta_com_bat(:,iC) = vol_dta_com_bat(:,iC) ./ vol_dta_com_bat(:,size(vol_dta_com_bat,2));
end
vol_roi_nme = vol_roi_nme(:, 1:size(vol_roi_nme,2)-1);
vol_dta_use = vol_dta_use(:, 1:size(vol_dta_use,2)-1);
vol_dta_com_bat = vol_dta_com_bat(:, 1:size(vol_dta_com_bat,2)-1);

vol_roi_nme = vol_roi_nme(:,[1:32 34:end]);
vol_dta_use = vol_dta_use(:,[1:32 34:end]);
vol_dta_com_bat = vol_dta_com_bat(:,[1:32 34:end]);

% Save data out
% find(strcmpi(fcfg.dta_nme,'rhs_postcentral'))
ind_out = [ 80 95 5 39 32 66 15 49 6 40 8 42 14 48 29 63 23 57 21 55];

out_dta = [ {'sbj_nme'} fcfg.dta_nme(ind_out)' ; epd_sbj_nme num2cell(com_bat_epd_dta(:,ind_out)) ];
for iC = 2:3
    out_dta(2:end,iC)     = num2cell( cell2mat(out_dta(2:end,iC)) ./ com_bat_epd_dta(:,end) );
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age/submission/Brain/Reviews/ComBat/ComBat_ROIs.csv', out_dta)

out_grp = [ {'sbj_nme'} grp_roi_nme(1:end-1) ; grp_sbj_nme grp_dta(:,1:2) ];
cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age/submission/Brain/Reviews/ComBat/ComBat_ROIs_groups.csv', out_grp)

%% ANCOVA original data - Omnibus
% Thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use(:,1);
fcfg.grp_nme = { 'Diagnosis_Thickness' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/orig'; % split % desikan

ejk_1way_ancova( fcfg );

% Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = vol_dta_use;
fcfg.dta_nme = vol_roi_nme;

fcfg.grp     = grp_var_use(:,1);
fcfg.grp_nme = { 'Diagnosis_Volume' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/orig'; % split % desikan

ejk_1way_ancova( fcfg );

%% ANCOVA original data - Onset
% Thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use(:,2);
fcfg.grp_nme = { 'Onset_Thickness' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/orig'; % split % desikan

ejk_1way_ancova( fcfg );

% Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = vol_dta_use;
fcfg.dta_nme = vol_roi_nme;

fcfg.grp     = grp_var_use(:,2);
fcfg.grp_nme = { 'Onset_Volumes' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/orig'; % split % desikan

ejk_1way_ancova( fcfg );

%% ANCOVA no covariate ComBat data - Omnibus
% Thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_com_bat;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use(:,1);
fcfg.grp_nme = { 'Diagnosis_Thickness' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/ComBat_ini'; % split % desikan

ejk_1way_ancova( fcfg );

% Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = vol_dta_com_bat;
fcfg.dta_nme = vol_roi_nme;

fcfg.grp     = grp_var_use(:,1);
fcfg.grp_nme = { 'Diagnosis_Volumes' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/ComBat_ini'; % split % desikan

ejk_1way_ancova( fcfg );

%% ANCOVA no covariate ComBat data - Onset
% Thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_com_bat;
fcfg.dta_nme = thk_roi_nme;

fcfg.grp     = grp_var_use(:,2);
fcfg.grp_nme = { 'Onset_Thickness' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/ComBat_ini'; % split % desikan

ejk_1way_ancova( fcfg );

% Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = vol_dta_com_bat;
fcfg.dta_nme = vol_roi_nme;

fcfg.grp     = grp_var_use(:,2);
fcfg.grp_nme = { 'Onset_Volumes' };

fcfg.cov     = cov_dta_use;
fcfg.cov_nme = cov_roi_nme(1:size(cov_dta_use,2));

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/desikan/ComBat_ini'; % split % desikan

ejk_1way_ancova( fcfg );





