clear; clc;

out_dir = '/home/ekaestne/PROJECTS/OUTPUT/ChristineHelp';

smt_stp = 10;
deg_fre = 44;

%%
lhs_sig = fs_load_mgh('/home/ekaestne/PROJECTS/OUTPUT/ChristineHelp/sigLH.mgh');
    lhs_pvl = 2*(1-normcdf(lhs_sig)); % Z-score -1.89
rhs_sig = fs_load_mgh('/home/ekaestne/PROJECTS/OUTPUT/ChristineHelp/sigRH.mgh');
    rhs_pvl = 2*(1-normcdf(rhs_sig));

figure()
subplot(2,1,1)
hist(lhs_sig(lhs_sig~=0),1000)
title('Left Hemisphere')
xlim([-3 6])
subplot(2,1,2)
hist(rhs_sig(rhs_sig~=0),1000)
xlim([-3 6])
title('Right Hemisphere')
print(gcf,'/home/ekaestne/PROJECTS/OUTPUT/ChristineHelp/distribution.png','-dpng')
close all

%% Plot pattern
% p05
low_rng_num = [ -1.96 1.96 ];
hgh_rng_num = [ -5 5 ];

pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction_p05' ];

pcfg.plt_dta = { lhs_sig' rhs_sig' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

% p01
low_rng_num = [ -2.57 2.57 ];
hgh_rng_num = [ -5 5 ];

pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction_p01' ];

pcfg.plt_dta = { lhs_sig' rhs_sig' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

%% Cluster Correct
% fs_calc_cluster_thresh %%%%%%%%%%%%%%
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/home/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/home/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

% p01
fcfg = [];
fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = smt_stp;
fcfg.df     = deg_fre;
fcfg.alpha  = .05;
fcfg.pval   = .01;

cls_thr_p01 = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ...
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );

pvl_p01_nme = num2str(roundsd(fcfg.pval,3));
pvl_p01_nme = pvl_p01_nme(3:end);
cls_p01_nme = num2str(round(cls_thr_p01));

% p05
fcfg = [];
fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = smt_stp;
fcfg.df     = deg_fre;
fcfg.alpha  = .05;
fcfg.pval   = .05;

cls_thr_p05 = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ...
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );

pvl_p05_nme = num2str(roundsd(fcfg.pval,3));
pvl_p05_nme = pvl_p05_nme(3:end);
cls_p05_nme = num2str(round(cls_thr_p05));

%% Plot Cluster Corrected p05
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = cls_thr_p05; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , lhs_pvl );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , rhs_pvl );

%
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

lhs_sig_lhs_cls = lhs_sig;
lhs_sig_rhs_cls = rhs_sig;

lhs_sig_lhs_cls(lhs_bdd_ind) = 0;
lhs_sig_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction' '_' 'cluster' '_pvl'  pvl_p05_nme '_cls' cls_p05_nme];

pcfg.plt_dta = { lhs_sig_lhs_cls' lhs_sig_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

%% Plot Cluster Corrected p01
fcfg = [];
fcfg.pvl_thr = .01;
fcfg.cls_thr = cls_thr_p01; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , lhs_pvl );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , rhs_pvl );

%
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

lhs_sig_lhs_cls = lhs_sig;
lhs_sig_rhs_cls = rhs_sig;

lhs_sig_lhs_cls(lhs_bdd_ind) = 0;
lhs_sig_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction' '_' 'cluster' '_pvl'  pvl_p01_nme '_cls' cls_p01_nme];

pcfg.plt_dta = { lhs_sig_lhs_cls' lhs_sig_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

%% Plot Cluster Corrected with 10
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 10; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , lhs_pvl );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , rhs_pvl );

%
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

lhs_sig_lhs_cls = lhs_sig;
lhs_sig_rhs_cls = rhs_sig;

lhs_sig_lhs_cls(lhs_bdd_ind) = 0;
lhs_sig_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction' '_' 'cluster' '_pvl'  pvl_p05_nme '_cls' num2str(10)];

pcfg.plt_dta = { lhs_sig_lhs_cls' lhs_sig_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)



