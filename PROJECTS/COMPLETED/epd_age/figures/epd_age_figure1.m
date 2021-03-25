%%
clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_final_sample_2.csv';
cov_nme = 'Epilepsy_Aging_Final_Sample_ASC_edit.csv';

epd_non_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
epd_usd_dir = '/home/mmilmcdRSI/data/fsurf';
adn_dir     = '/home/mmilmcdRSI/data_ADNI/fsurf/';

%% Constants
pvl_chs = .05;
smt_stp = 176;
deg_fre = 148;

pvl_cls = .05;

low_rng_num = [ -0.01 nan ];
hgh_rng_num = [ -0.25 nan ];

plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/figures/figure1/';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

%% Cluster threshold
srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

% 3549 2150
fcfg = [];

fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(smt_stp);
fcfg.df     = deg_fre;
fcfg.alpha  = pvl_cls;
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
plt_nme = { 'TLE_HC' 'MCI_HC' };
smt_stp = { '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'pValueL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'pValueR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'emmeanL' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean * -1;
        men_dff_rhs = load([ jhn_dta_loc '/' 'emmeanR' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean * -1;
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr=0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr=0; end
        
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
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
        fcfg.hms     = 'lh';
        pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs_176 );
        
        fcfg.hms     = 'rh';
        pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs_176 );
        
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
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end

%% Colorbar
top_pct = 1;
% Make Colormap
cfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
cfg.fmr_col_map = cellfun(@rgb ,cfg.fmr_col_map,'uni',0);

fmr_col_map = [];
for iC = 1:numel(cfg.fmr_col_map)-1
    fmr_col_map = [fmr_col_map ; [linspace(cfg.fmr_col_map{iC}(1),cfg.fmr_col_map{iC+1}(1),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))' linspace(cfg.fmr_col_map{iC}(2),cfg.fmr_col_map{iC+1}(2),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))' linspace(cfg.fmr_col_map{iC}(3),cfg.fmr_col_map{iC+1}(3),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))']; ];
end

ttt = fmr_col_map(1:500,:);

pcfg = [];
pcfg.col_map = ttt;
pcfg.col_bar = [0 0.5];
pcfg.out_dir = plt_out;
pcfg.sve_pre = 'figure1';
mmil_color_bar(pcfg)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quick Check Surface Map Plots
plt_nme = { 'TLE_HC' 'MCI_HC' 'Late_HC' 'Early_HC' };
smt_stp = { '176' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         pvl_lhs = load([ jhn_dta_loc '/' 'orig' '/' 'pValueL' '_' plt_nme{iPL} '_sm313.mat' ]);
%         pvl_lhs = pvl_lhs.pValue;
%         pvl_rhs = load([ jhn_dta_loc '/' 'orig' '/' 'pValueR' '_' plt_nme{iPL} '_sm313.mat' ]);
%         pvl_rhs = pvl_rhs.pValue;
%         
        pvl_lhs_176 = load([ jhn_dta_loc '/' 'sm176/par/' '/' 'pValue_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        pvl_lhs_176 = pvl_lhs_176.pValue;
        pvl_rhs_176 = load([ jhn_dta_loc '/' 'sm176/par/' '/' 'pValue_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        pvl_rhs_176 = pvl_rhs_176.pValue;
        
        pvl_lhs_313 = load([ jhn_dta_loc '/' 'sm313/par/' '/' 'pValue_lhs' '_' plt_nme{iPL} '.mat' ]);
        pvl_lhs_313 = pvl_lhs_313.pValue;
        pvl_rhs_313 = load([ jhn_dta_loc '/' 'sm313/par/' '/' 'pValue_rhs' '_' plt_nme{iPL} '.mat' ]);
        pvl_rhs_313 = pvl_rhs_313.pValue;
        
        scatter(pvl_lhs_176, pvl_lhs_313)
        
        [ ~ , pnt_176 ] = min(pvl_lhs_176)
        [ ~ , pnt_313 ] = min(pvl_lhs_313)
        
        % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         men_dff_lhs = load([ jhn_dta_loc '/' 'orig' '/' 'emmeanL' '_' plt_nme{iPL} '_sm313.mat' ]);
%         men_dff_lhs = men_dff_lhs.emmean;
%         men_dff_rhs = load([ jhn_dta_loc '/' 'orig' '/' 'emmeanR' '_' plt_nme{iPL} '_sm313.mat' ]);
%         men_dff_rhs = men_dff_rhs.emmean;
        
        men_dff_lhs_176 = load([ jhn_dta_loc '/' 'sm176/par/' '/' 'emmean_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs_176 = men_dff_lhs_176.emmean;
        men_dff_rhs_176 = load([ jhn_dta_loc '/' 'sm176/par/' '/' 'emmean_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs_176 = men_dff_rhs_176.emmean;

        men_dff_lhs_313 = load([ jhn_dta_loc '/' 'sm313/par/' '/' 'emmean_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs_313 = men_dff_lhs_313.emmean;
        men_dff_rhs_313 = load([ jhn_dta_loc '/' 'sm313/par/' '/' 'emmean_rhs' '_' plt_nme{iPL} '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs_313 = men_dff_rhs_313.emmean;
        
        scatter(men_dff_lhs_176, men_dff_lhs_313)
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMMEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        low_rng_num = [ -0.01 nan ];
        hgh_rng_num = [ -0.25 nan ];
        
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs_176' men_dff_rhs_176' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs_313' men_dff_rhs_313' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % PVALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        low_rng_num = [ -0.95 nan ];
        hgh_rng_num = [ -1.00 nan ];
        
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { (1-pvl_rhs_176)'*-1 (1-pvl_lhs_176)'*-1 };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { (1-pvl_rhs_313)'*-1 (1-pvl_lhs_313)'*-1 };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs_176,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr=0; end
        pvl_rhs_fdr = FDR(pvl_rhs_176,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr=0; end
        
        lhs_bdd_ind = find(pvl_lhs_176>pvl_lhs_fdr);
        rhs_bdd_ind = find(pvl_rhs_176>pvl_rhs_fdr);
        
        men_dff_lhs_fdr = men_dff_lhs;
        men_dff_rhs_fdr = men_dff_rhs;
        
        men_dff_lhs_fdr(lhs_bdd_ind) = 0;
        men_dff_rhs_fdr(rhs_bdd_ind) = 0;
          
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'fdr' ];
        
        pcfg.plt_dta = { men_dff_rhs_fdr' men_dff_lhs_fdr' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs_org' men_dff_rhs_org' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
        fcfg.hms     = 'lh';
        pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs_176 );
        
        fcfg.hms     = 'rh';
        pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs_176 );
        
        % 
        lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
        rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);
        
        men_dff_lhs_cls = men_dff_lhs_176;
        men_dff_rhs_cls = men_dff_rhs_176;
        
        men_dff_lhs_cls(lhs_bdd_ind) = 0;
        men_dff_rhs_cls(rhs_bdd_ind) = 0;

        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'cluster' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end





