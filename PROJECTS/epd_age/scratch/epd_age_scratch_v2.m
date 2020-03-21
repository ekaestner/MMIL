clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_final_sample_2.csv';
cov_nme = 'Epilepsy_Aging_Final_Sample_ASC_edit.csv';

epd_non_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
epd_usd_dir = '/home/mmilmcdRSI/data/fsurf';
adn_dir     = '/home/mmilmcdRSI/data_ADNI/fsurf/';

%% Category of Subjects - Updated to 80
sbj_nme = mmil_readtext( [prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' sbj_nme ]);
sbj_nme_lbl = sbj_nme;

dta_loc(string_find(sbj_nme(:,1) , {'epd'}),1) = repmat({epd_usd_dir},numel(string_find(sbj_nme(:,1) , {'epd'})),1);
dta_loc(setdiff(string_find(sbj_nme(:,3) , {'EPD_Old'}),string_find(sbj_nme(:,1) , {'epd'})),1) = repmat({epd_non_dir},numel(setdiff(string_find(sbj_nme(:,3) , {'EPD_Old'}),string_find(sbj_nme(:,1) , {'epd'}))),1);

dta_loc(string_find(sbj_nme(:,3) , {'MCI'}),1) = repmat({adn_dir},numel(string_find(sbj_nme(:,3) , {'MCI'})),1);

dta_loc(string_find(sbj_nme(:,1) , {'fc'}),1) = repmat({epd_usd_dir},numel(string_find(sbj_nme(:,1) , {'fc'})),1);
dta_loc(setdiff(string_find(sbj_nme(:,3) , {'HC'}),string_find(sbj_nme(:,1) , {'fc'})),1) = repmat({adn_dir},numel(setdiff(string_find(sbj_nme(:,3) , {'HC'}),string_find(sbj_nme(:,1) , {'fc'}))),1);

%% Organize Covariates
cov_tbl     = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' cov_nme ]);
cov_tbl_lbl = cov_tbl(1,:);
cov_tbl_sbj = cov_tbl(2:end,1);
cov_tbl     = cell2mat(cov_tbl(2:end,2:end));

cov_col_ind = 2:5;

cov_out_tbl = [sbj_nme(2:end,1) cell(size(sbj_nme(2:end,1),1),4)];
cov_out_lbl = [ sbj_nme_lbl{1} cov_tbl_lbl(cov_col_ind) ];
for iR = 1:size(cov_tbl,1)

    sbj_ind = find(strcmpi(cov_tbl_sbj,cov_out_tbl{iR,1}));
    
    if isempty(sbj_ind)
        error('')
    else
        cov_out_tbl(iR,2:5) = num2cell(cov_tbl(sbj_ind,:));
    end    
    
end
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates.csv'],[cov_out_lbl ; cov_out_tbl])

%% Load Surface
% Missing -> UCSD_18, UCSD_9

prj_dir = '/home/ekaestne/PROJECTS/';

hms     = {'lhs' 'rhs'};

tic;
fcfg = [];

fcfg.prj_dir = prj_dir;

fcfg.mes_typ = 'aMRI_thickness';
fcfg.smt_stp  = 2819;

fcfg.anl_dir = 'analysis';

for iH = 1:numel(hms)
    
    srf_dta  = nan(size(sbj_nme,1)-1,163842);
    fcfg.hms = hms{iH}(1:2);
    
    for iS = 2:size(sbj_nme,1)
        
        fcfg.prc_dir     = dta_loc{iS,1};
        fcfg.sbj_fsr_dir = sbj_nme{iS,2};
        
        srf_dta_hld = ejk_extract_vertices(fcfg);
        if ~isempty(srf_dta_hld); srf_dta(iS-1,:) = srf_dta_hld; end
        
    end
    srf_dta_sbj = sbj_nme(2:end,1);
    save([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.csv'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% Fix Nan
% lhs
load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs.mat'  ]);
mss_ind_lhs = find(isnan(srf_dta(:,100)));

srf_dta(mss_ind_lhs,:)      = [];
srf_dta_sbj(mss_ind_lhs,:)  = [];

save( ['/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_no_nan.mat'] , 'srf_dta' , 'srf_dta_sbj' )
clear srf_dta srf_dta_sbj

% rhs
load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs.mat'  ]);
mss_ind_rhs = find(isnan(srf_dta(:,100)));

srf_dta(mss_ind_rhs,:)      = [];
srf_dta_sbj(mss_ind_rhs,:)  = [];

save( ['/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_no_nan.mat'] , 'srf_dta' , 'srf_dta_sbj' )
clear srf_dta srf_dta_sbj

if ~all(mss_ind_lhs == mss_ind_rhs); error(''); end

% covariates
cov_tbl = mmil_readtext([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates.csv']);
cov_tbl(mss_ind_lhs+1,:) = [];
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates_no_nan.csv'],cov_tbl)

% sbj_nme
sbj_nme_out = sbj_nme;
sbj_nme_out(mss_ind_lhs+1,:) = [];
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'sbj_nme_no_nan.csv'],sbj_nme_out)

%% Surface Map Plots
%
lhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_no_nan.mat'  ]);
rhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_no_nan.mat'  ]);

%
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.sbj_grp     = sbj_nme;
fcfg.sbj_grp_col = { 'Diagnosis' };
fcfg.sbj_grp_nme = { { 'HC' 'EPD_Old' 'MCI'} };
fcfg.sbj_grp_cmp = { { [2 1] [3 1] [2 3] } };

fcfg.dta_lbl = { 'corticalthickness' };
fcfg.dta     = { lhs_dta rhs_dta }; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

fcfg.hms     = {'lhs' 'rhs'};

ejk_surf_group_avg(fcfg)

%% ROI load
prj_dir = '/home/ekaestne/PROJECTS/';
hms     = {'lhs' 'rhs'};
prc_nme = { '' '.a2009s' };

tic;

fcfg = [];

fcfg.prj_dir = prj_dir;

fcfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md

fcfg.anl_dir = 'analysis';

for iPR = 1:numel(prc_nme)
    
    fcfg.prc_nme = prc_nme{iPR};
    
    for iS = 2:size(sbj_nme,1)
        
        fcfg.ovr_dir     = dta_loc{iS,1};
        fcfg.sbj_fsr_dir = sbj_nme{iS,2};
        
        [ gry_thk_dta_hld , tot_lbl ]= ejk_extract_grey_thickness(fcfg);
        try
            if ~isempty(gry_thk_dta_hld)
                gry_thk_dta(iS-1,:) = gry_thk_dta_hld;
                fprintf( [ sbj_nme{iS,2} ' : ' 'Data Loaded\n' ] )
            else
                gry_thk_dta(iS-1,:) = nan(1,size(gry_thk_dta,2));
                fprintf( [ sbj_nme{iS,2} ' : ' 'Missing\n' ] )
            end
        catch
            gry_thk_dta(iS-1,:) = nan(1,size(gry_thk_dta,2));
            fprintf( [ sbj_nme{iS,2} ' : ' 'Missing\n' ] )
        end
    end
    
    sve_sbj_nme = sbj_nme(2:end,1);
    cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI' '_' 'thickness' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(gry_thk_dta,1))] [tot_lbl ; num2cell(gry_thk_dta)] ]);
    clear gry_thk_dta tot_lbl
    
end
toc

%% Plot ROIs
prj_dir = '/home/ekaestne/PROJECTS/';

%
dkn_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI_thickness_aparc_.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme(2:end,:);
fcfg.fle_nme = dkn_dta;
dkn_dta = ejk_load_mcd_data(fcfg);

vol_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Volumes_aparc.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme(2:end,:);
fcfg.fle_nme = vol_dta;
vol_dta = ejk_load_mcd_data(fcfg);

% BAR
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'volumes_new' 'desikan_new' };
fcfg.ydt     = {  vol_dta      dkn_dta  };

fcfg.grp     = sbj_nme;
fcfg.grp_col = [ 3 ];
fcfg.nme_col = [ sbj_nme(1,3) ];

fcfg.grp_nme = { { 'HC' 'EPD_Old' 'MCI' } };
fcfg.grp_clr = { { 'black' 'orange' 'blue' } ...
                 { 'black' 'bluish grey' 'light blue' 'dark blue' } };
fcfg.xdt     = { [ 1 2 3] };

fcfg.plt_cmp = { { [1 2] [1 3] [2 3] } };
fcfg.plt_anv = { [1 2 3] };

ejk_roi_bar(fcfg)

%% Load Volumes
prj_dir = '/home/ekaestne/PROJECTS/';

%
tic;
fcfg = [];

fcfg.prj_dir = prj_dir;

for iS = 2:size(sbj_nme,1)
    
    fcfg.ovr_dir     = dta_loc{iS,1};
    fcfg.sbj_fsr_dir = sbj_nme{iS,2};
    
    [ vol_dta_hld , tot_lbl] = ejk_extract_volumes(fcfg);
    if ~isempty(vol_dta_hld); vol_dta(iS-1,:) = vol_dta_hld; else vol_dta(iS-1,:) = nan(1,size(vol_dta,2)); end
    
end

sve_sbj_nme = sbj_nme(2:end,1);
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'Volumes' '_' 'aparc' '.csv'],[ ['SubjId' ; sve_sbj_nme] [tot_lbl ; num2cell(vol_dta)] ]);
clear vol_dta tot_lbl

toc

%% Plot Volumes

