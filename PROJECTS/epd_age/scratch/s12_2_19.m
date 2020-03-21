%%
clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_final_sample_2.csv';
cov_nme = 'Epilepsy_Aging_Final_Sample_ASC_edit.csv';

epd_non_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
epd_usd_dir = '/home/mmilmcdRSI/data/fsurf';
adn_dir     = '/home/mmilmcdRSI/data_ADNI/fsurf/';

%% Cluster threshold
pvl_chs = .001;

srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

%
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
plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/12_02_19/';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme = { 'TLE_HC' 'MCI_HC' 'TLE_MCI' 'Early_HC' 'Late_HC' 'Early_Late' 'LTLE_HC' 'RTLE_HC' };
smt_stp = { '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean*-1; men_dff_lhs = men_dff_lhs(1:end-1);
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean*-1;
        
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
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr = 0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr = 0; end
        
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

%% Check clustering
plt_nme = { 'TLE_HC' 'MCI_HC' 'TLE_MCI' 'Early_HC' 'Late_HC' 'Early_Late' 'LTLE_HC' 'RTLE_HC' };
smt_stp = { '313' };

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_chs = .05;

srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

% 
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

% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%
for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean*-1;
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean*-1;

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
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'cluster' '_p' pvl_nme '_' cls_nme 'mm'];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = [ -0.01 0.01 ];
        pcfg.hgh_rng_num = [ -0.20 0.20 ];
        
        mmil_anat_surf_plot(pcfg)
        
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)
hist(pvl_lhs,10000);
xlim([-.001 .01]);
subplot(2,1,2)
hist(pvl_rhs,10000);
xlim([-.001 .01]);

min(pvl_lhs)
    sum(pvl_lhs<.01)
min(pvl_rhs)
    sum(pvl_rhs<.01)

%% Education
plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/12_02_19/education';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme = { 'TLE_HC' };
smt_stp = { '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_ed.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_ed.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_ed.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean*-1; men_dff_lhs = men_dff_lhs(1:end-1);
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_ed.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean*-1;
        
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
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr = 0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr = 0; end
        
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

%% Outliers, Check Vertices
% Data
sbj_nme     = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/epd_age_covariates_no_nan.csv');
    sbj_nme(cellfun(@isempty,sbj_nme(:,10)),10) = {''};
ely_ind = find(strcmpi(sbj_nme(:,10),'Early'))-1;
lte_ind = find(strcmpi(sbj_nme(:,10),'Late'))-1;
con_ind = find(strcmpi(sbj_nme(:,6),'HC'))-1;
mci_ind = find(strcmpi(sbj_nme(:,6),'MCI'))-1;
epd_ind = find(strcmpi(sbj_nme(:,6),'EPD_Old'))-1;

org_thk_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_no_nan.mat');
org_thk_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_rhs_no_nan.mat');

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Late_HC.mat');
    pvl_lhs_lte = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Late_HC.mat');
    pvl_rhs_lte = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Late_HC.mat');
    men_dff_lhs_lte = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Late_HC.mat');
    men_dff_rhs_lte = men_dff_rhs.emmean;

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Early_HC.mat');
    pvl_lhs_ely = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Early_HC.mat');
    pvl_rhs_ely = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Early_HC.mat');
    men_dff_lhs_ely = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Early_HC.mat');
    men_dff_rhs_ely = men_dff_rhs.emmean;

%
tbl = [  sbj_nme(epd_ind+1,1) num2cell(org_thk_lhs.srf_dta(epd_ind,thk_vtx_lhs)) sbj_nme(epd_ind+1,4) sbj_nme(epd_ind+1,10)] ;
[ ~ , tbl_srt] = sort(cell2mat(tbl(:,2)));
tbl = tbl(tbl_srt,:);
tbl(cell2mat(tbl(:,2))>max(org_thk_lhs.srf_dta(con_ind,thk_vtx_lhs)),:)

figure();
scatter( cell2mat(tbl(:,2)) , cell2mat(tbl(:,3)) , 30 , 'MarkerFaceColor', rgb('black') ); hold on;
    cof = fit( cell2mat(tbl([1:33 35:end],2)) , cell2mat(tbl([1:33 35:end],3)) , 'poly1' );
    plot( cof , cell2mat(tbl(:,2)) , cell2mat(tbl(:,3)) );

%% Etiology
eto_csv = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/projects/Epilepsy_and_Aging/etiology.csv');

nme = { 'Normal' 'MTS' 'Tumor' 'Vascular' 'Other' 'Unknown' };

ely_nme = find(strcmpi( tbl(:,4) , 'Early' ));
    ely_min_hld = 1 - 0.15;
    ely_max_hld = 1 + 0.15;
lte_nme = find(strcmpi( tbl(:,4) , 'Late' ));
    lte_min_hld = 2 - 0.15;
    lte_max_hld = 2 + 0.15;

plt_col = distinguishable_colors(numel(nme));
    
figure(); hold on;
for iN = 1:numel(nme)
    
    ely_nme_ind = intersect( find(strcmpi( eto_csv(:,2) , nme{iN} )) , ely_nme );
        ely_xdt = (ely_max_hld-ely_min_hld) .* rand(numel(ely_nme_ind),1) + ely_min_hld;
        ely_ydt = cell2mat(tbl(ely_nme_ind,2));
    lte_nme_ind = intersect( find(strcmpi( eto_csv(:,2) , nme{iN} )) , lte_nme );
        lte_xdt = (lte_max_hld-lte_min_hld) .* rand(numel(lte_nme_ind),1) + lte_min_hld;
        lte_ydt = cell2mat(tbl(lte_nme_ind,2));
            
    scatter( ely_xdt , ely_ydt , 30 , 'MarkerFaceColor', plt_col(iN,:) );
    scatter( lte_xdt , lte_ydt , 30 , 'MarkerFaceColor', plt_col(iN,:) );
    
end

figure(); hold on;
num_hld = 1; text( 0 , num_hld , nme{num_hld} , 

%% Education
plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/12_02_19/education/education_thinning/';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme = { 'TLE_HC' };...'MCI_HC' 'Early_HC' 'Late_HC' 'Early_Late' 'TLE_MCI'  };
smt_stp = { '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_sm313_ed.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_sm313_ed.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_sm313_ed.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean*-1; men_dff_lhs = men_dff_lhs(1:end-1);
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_sm313_ed.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean*-1;
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = [ -0.01 nan ];
        pcfg.hgh_rng_num = [ -0.20 nan ];
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr = 0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr = 0; end
        
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
        
        pcfg.low_rng_num = [ -0.01 nan ];
        pcfg.hgh_rng_num = [ -0.20 nan ];
        
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
        
        pcfg.low_rng_num = [ -0.01 nan ];
        pcfg.hgh_rng_num = [ -0.20 nan ];
        
        mmil_anat_surf_plot(pcfg)
        
    end
end



