%%
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

cov_out_tbl = [sbj_nme(:,1) cell(size(sbj_nme(2:end,1),1),4)];
cov_out_lbl = [ sbj_nme_lbl{1} cov_tbl_lbl(cov_col_ind) ];
for iR = 1:size(cov_tbl,1)

    sbj_ind = find(strcmpi(cov_tbl_sbj,cov_out_tbl{iR,1}));
    
    if isempty(sbj_ind)
        error('')
    else
        cov_out_tbl(iR,2:5) = num2cell(cov_tbl(sbj_ind,:));
    end    
    
end
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates.csv'],[cov_out_lbl ; cov_out_tbl]);

%% Load Surface
smt_stp = [ 176 313 705 2819 ];

prj_dir = '/home/ekaestne/PROJECTS/';

hms     = {'lhs' 'rhs'};

for iSM = [1 3 4]
    tic;
    fcfg = [];
    
    fcfg.prj_dir = prj_dir;
    
    fcfg.mes_typ = 'aMRI_thickness';
    fcfg.smt_stp  = smt_stp(iSM);
    
    fcfg.anl_dir = 'analysis';
    
    for iH = 1:numel(hms)
        
        srf_dta  = nan(size(sbj_nme,1),163842);
        fcfg.hms = hms{iH}(1:2);
        
        for iS = 1:size(sbj_nme,1)
            
            fcfg.prc_dir     = dta_loc{iS,1};
            fcfg.sbj_fsr_dir = sbj_nme{iS,2};
            
            srf_dta_hld = ejk_extract_vertices(fcfg);
            if ~isempty(srf_dta_hld); srf_dta(iS,:) = srf_dta_hld; end
            
        end
        srf_dta_sbj = sbj_nme(1:end,1);
        save([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'surf' '_' fcfg.mes_typ '_' fcfg.hms 's' '_' 'sm' num2str(fcfg.smt_stp) '.mat'],'srf_dta_sbj','srf_dta');
        clear srf_dta
        
    end
    toc
end

%% Quick Check Surface Map Plots
smt_stp = [ 176 313 705 2819 ];

for iSM = 1:numel(smt_stp)
    
    sbj_nme_org = [ {'SbjID'} {['Diagnosis_' num2str(smt_stp(iSM))]} ; sbj_nme(:,[1 3]) ];
    
    %
    lhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    rhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    
    %
    fcfg = [];
    
    fcfg.prj_dir = prj_dir;
    fcfg.prj_nme = prj_nme;
    
    fcfg.sbj_grp     = sbj_nme_org;
    fcfg.sbj_grp_col = { ['Diagnosis_' num2str(smt_stp(iSM))] };
    fcfg.sbj_grp_nme = { { 'HC' 'EPD_Old' 'MCI'} };
    fcfg.sbj_grp_cmp = { { [2 1] [3 1] [2 3] } };
    
    fcfg.dta_lbl = { 'corticalthickness' };
    fcfg.dta     = { lhs_dta rhs_dta }; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'
    
    fcfg.hms     = {'lhs' 'rhs'};
    
    ejk_surf_group_avg(fcfg)
    
end

%% Fix Nan
smt_stp = [ 176 313 705 2819 ];

for iSM = 1:numel(smt_stp)
    
    % Total missing indices
    sbj_nme_org = [ {'SbjID'} {['Diagnosis_' num2str(smt_stp(iSM))]} ; sbj_nme(:,[1 3]) ];
    
    load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    mss_ind_lhs = find(isnan(srf_dta(:,100)));
    
    load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    mss_ind_rhs = find(isnan(srf_dta(:,100)));
    
    tot_mss_ind = unique([ mss_ind_lhs ; mss_ind_rhs]);
    
    % lhs
    load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    
    srf_dta(tot_mss_ind,:)      = [];
    srf_dta_sbj(tot_mss_ind,:)  = [];
    
    save( ['/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs_sm' num2str(smt_stp(iSM)) '_no_nan.mat'] , 'srf_dta' , 'srf_dta_sbj' )
    clear srf_dta srf_dta_sbj
    
    % rhs
    load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_sm' num2str(smt_stp(iSM)) '.mat'  ]);
    
    srf_dta(tot_mss_ind,:)      = [];
    srf_dta_sbj(tot_mss_ind,:)  = [];
    
    save( ['/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs_sm' num2str(smt_stp(iSM)) '_no_nan.mat'] , 'srf_dta' , 'srf_dta_sbj' )
    clear srf_dta srf_dta_sbj
    
    % covariates
    cov_tbl = mmil_readtext([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates.csv']);
    cov_tbl(tot_mss_ind+1,:) = [];
    cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'epd_age_covariates_sm' num2str(smt_stp(iSM)) '_no_nan.csv'],cov_tbl)
    
    % sbj_nme
    sbj_nme_out = sbj_nme_org;
    sbj_nme_out(tot_mss_ind+1,:) = [];
    cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'sbj_nme_sm' num2str(smt_stp(iSM)) '_no_nan.csv'],sbj_nme_out)
    
end

%% Cluster threshold
pvl_chs = .01;

srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

% 3549 2150
fcfg = [];

fcfg.nverts = 163842-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(313);
fcfg.df     = 144;
fcfg.alpha  = .05;
fcfg.pval   = pvl_chs;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ... 
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );

pvl_nme = num2str(pvl_chs); pvl_nme = pvl_nme(3:end);
cls_nme = num2str(round(cls_thr));

%% Quick Check Surface Map Plots
plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/11_26_19/';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme = { 'TLE_HC' };
smt_stp = { '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean; men_dff_lhs = men_dff_lhs(1:end-1);
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean;
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = [ -0.01 0.01 ];
        pcfg.hgh_rng_num = [ -0.20 0.20 ];
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs,.05);
        pvl_rhs_fdr = FDR(pvl_rhs,.05);
        
        lhs_bdd_ind = find(pvl_lhs>pvl_lhs_fdr);
        rhs_bdd_ind = find(pvl_rhs>pvl_rhs_fdr);
        
        men_dff_lhs_fdr = men_dff_lhs;
        men_dff_rhs_fdr = men_dff_rhs;
        
        men_dff_lhs_fdr(lhs_bdd_ind) = 0;
        men_dff_rhs_fdr(rhs_bdd_ind) = 0;
          
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'fdr' ];
        
        pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = [ -0.01 0.01 ];
        pcfg.hgh_rng_num = [ -0.20 0.20 ];
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
        fcfg.hms     = 'lh';
        pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );
        
        fcfg.hms     = 'rh';
        pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );
        
        % 
        lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
        rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);
        
        men_dff_lhs_cls = men_dff_lhs;
        men_dff_rhs_cls = men_dff_rhs;
        
        men_dff_lhs_cls(lhs_bdd_ind) = 0;
        men_dff_rhs_cls(rhs_bdd_ind) = 0;

        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'cluster' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = [ -0.01 0.01 ];
        pcfg.hgh_rng_num = [ -0.20 0.20 ];
        
        mmil_anat_surf_plot(pcfg)
        
    end
end
















