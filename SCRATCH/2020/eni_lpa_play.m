clear; clc;

%%
roi_int_org = { 'CGC_L'    'CGC_R'      'CGH_L'    'CGH_R'      'EC_L'    'EC_R'      'SS_L'    'SS_R'      'UNC_L'    'UNC_R' };
roi_int_sde = { 'CGC_ipsi' 'CGC_contra' 'CGH_ipsi' 'CGH_contra' 'EC_ipsi' 'EC_contra' 'SS_ipsi' 'SS_contra' 'UNC_ipsi' 'UNC_contra' };

%% Put Together Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_dti_dfa = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Data/FA.csv');

% Covariates
cov_dta     = dta_dti_dfa(:,1:15);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

% Data
dta_dta     = dta_dti_dfa(:,[1 16:end]);
dta_sbj_nme = dta_dta(2:end,1);
dta_roi_nme = dta_dta(1,2:end);
dta_dta     = cell2mat(dta_dta(2:end,2:end));

%% Data Cleaning 
% Convert 'NA' to nan
for iR = 1:size(cov_dta,1)
    for iC = 1:size(cov_dta,2)
        if strcmpi(cov_dta{iR,iC},'NA')
            cov_dta{iR,iC} = -999;
        end
    end
end

% Labeled Covariates
cov_dta_lbl_nme = { 'Dx' 'SDx' 'Sex' 'Handedness' };
cov_dta_lbl = enigma_relabel(cov_dta);

% Check data for missing values
sum(isnan(dta_dta(:)))

%% Data Explore
con_ind = find(cell2mat(cov_dta(:,2))==0);
epl_ind = find(cell2mat(cov_dta(:,2))==1);

% Diagnoses
tle_lft_tot = find( cell2mat(cov_dta(:,3))==3 | cell2mat(cov_dta(:,3))==5 );
tle_rgh_tot = find( cell2mat(cov_dta(:,3))==4 | cell2mat(cov_dta(:,3))==6 );
tle_ext     = find( cell2mat(cov_dta(:,3))==7 );

% Drug Responsitivity
drg_sze_gdd = find( cell2mat(cov_dta(:,10))==1 );
drg_sze_mid = find( cell2mat(cov_dta(:,10))==2 );
drg_sze_bdd = find( cell2mat(cov_dta(:,10))==3 );

[ tbl, ~, ~, tbl_lbl] = crosstab( cell2mat(cov_dta(:,10)), cov_dta_lbl(:,2) );
cell2csv('/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play/DrugResponse.csv', [ {''} tbl_lbl(:,2)' ;  tbl_lbl(1:end-1,1) num2cell(tbl) ])

% ENGEL
eng_one = find( cell2mat(cov_dta(:,11))==1 );
eng_two = find( cell2mat(cov_dta(:,11))==2 );
eng_thr = find( cell2mat(cov_dta(:,11))==3 );
eng_for = find( cell2mat(cov_dta(:,11))==4 );
eng_mss = find( cell2mat(cov_dta(:,11))==5 );

[ tbl, ~, ~, tbl_lbl] = crosstab( cell2mat(cov_dta(:,11)), cov_dta_lbl(:,2) );
cell2csv('/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play/EngelResponse.csv', [ {''} tbl_lbl(:,2)' ;  tbl_lbl(1:end,1) num2cell(tbl) ])

%% ComBat
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme;

fcfg.dta     = dta_dta;
fcfg.dta_nme = cellfun( @(x) mmil_spec_char(x,{'-' '/'}), dta_roi_nme, 'uni', 0);

fcfg.btc     = cov_dta(:,1);
fcfg.btc_nme = cov_roi_nme(:,1);

fcfg.cov     = cov_dta(:,2);
fcfg.cov_nme = cov_roi_nme(2);

fcfg.ylm     = [0 1];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play/';

com_bat_epd_dta = ejk_ComBat(fcfg);

% Create to ipsi/contra
[ lat_dta_sbj_nme, lat_dta_dta, lat_dta_roi_nme, lat_cov_sbj_nme, lat_cov_dta, lat_cov_dta_lbl] = enigma_ipsi_contra( dta_sbj_nme, com_bat_epd_dta, dta_roi_nme, cov_sbj_nme, cov_dta, cov_dta_lbl);

%% PCA

%% LPA - Original
% Select Data
dta_nme = cellfun( @(x) mmil_spec_char(x,{'-' '/'}), dta_roi_nme, 'uni', 0);

[ ~, roi_ind ] = intersect( dta_nme, roi_int_org );
dta_use = com_bat_epd_dta( sort([ tle_lft_tot ; tle_rgh_tot]), roi_ind );

% Fit 
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme(sort([ tle_lft_tot ; tle_rgh_tot]));

fcfg.dta     = dta_use;
fcfg.dta_nme = dta_nme( roi_ind );

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play/';

ejk_fit_lpa( fcfg );

% LPA Model
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme(sort([ tle_lft_tot ; tle_rgh_tot]));

fcfg.mdl_num = 8;

fcfg.dta     = dta_use;
fcfg.dta_nme = dta_nme( roi_ind );

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play/';

ejk_lpa( fcfg );

%% LPA - Ipsi/Contra
% Select Data
dta_nme = cellfun( @(x) mmil_spec_char(x,{'-' '/'}), lat_dta_roi_nme, 'uni', 0);

[ ~, roi_ind ] = intersect( dta_nme, roi_int_sde );
dta_use = lat_dta_dta( :, roi_ind );

% Fit 
fcfg = [];

fcfg.sbj_nme = lat_dta_sbj_nme;

fcfg.dta     = dta_use;
fcfg.dta_nme = dta_nme( roi_ind );

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_Side/';

ejk_fit_lpa( fcfg );

% LPA Model
fcfg = [];

fcfg.sbj_nme = dta_sbj_nme(sort([ tle_lft_tot ; tle_rgh_tot]));

fcfg.mdl_num = 7;

fcfg.dta     = dta_use;
fcfg.dta_nme = dta_nme( roi_ind );

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/eni_lpa/Play_Side/';

ejk_lpa( fcfg );














