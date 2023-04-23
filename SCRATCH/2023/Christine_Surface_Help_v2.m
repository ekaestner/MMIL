clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Christine';

smt_stp = 10;
deg_fre = 36;

%%
lhs_sig = fs_load_mgh([out_dir '/' 'lhs_sig.mgh']);
    lhs_pvl = 2*(1-normcdf(lhs_sig)); % Z-score
    
rhs_sig = fs_load_mgh([out_dir '/' 'rhs_sig.mgh']);
    rhs_pvl = 2*(1-normcdf(rhs_sig)); % Z-score
    
figure()
subplot(2,1,1)
hist(lhs_sig(lhs_sig~=0),200)
title('Left Hemisphere')
xlim([-2.25 3.25])
subplot(2,1,2)
hist(rhs_sig(rhs_sig~=0),200)
xlim([-2.25 3.25])
title('Right Hemisphere')

print(gcf,[out_dir '/' 'statistic_distribution.png'],'-dpng')
close all

%% Plot pattern
% match Christine plot
low_rng_num = [ -1.7 1.7 ];
hgh_rng_num = [ -5 5 ];

pcfg = [];

pcfg.hme_wrk = 2;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '01_statistic_1.7' ];

pcfg.plt_dta = { lhs_sig' rhs_sig' };

pcfg.srf_typ = 'inflated';

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

% p15
fcfg = [];
fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = smt_stp;
fcfg.df     = deg_fre;
fcfg.alpha  = .05;
fcfg.pval   = .05;

cls_thr_p15 = fs_calc_cluster_thresh( fcfg.nverts , ...
                                      fcfg.area , ...
                                      fcfg.fwhm , ...
                                      fcfg.df , ...
                                      fcfg.alpha , ...
                                      fcfg.pval );

pvl_p15_nme = num2str(roundsd(fcfg.alpha,3));
pvl_p15_nme = pvl_p15_nme(3:end);
cls_p15_nme = num2str(round(cls_thr_p15));

fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = cls_thr_p15; % mm^2
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
rhs_sig_rhs_cls = rhs_sig;

lhs_sig_lhs_cls(lhs_bdd_ind) = 0;
rhs_sig_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 2;

pcfg.plt_dta = { lhs_sig_lhs_cls' rhs_sig_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.srf_typ = 'inflated';

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ 'test_statistic_no_correction' '_' 'cluster' '_pvl'  pvl_p15_nme '_cls' num2str(cls_p15_nme) '_' pcfg.srf_typ];

mmil_anat_surf_plot(pcfg)


pcfg.srf_typ = 'pial';
pcfg.out_pre_fix = [ 'test_statistic_no_correction' '_' 'cluster' '_pvl'  pvl_p15_nme '_cls' num2str(cls_p15_nme) '_' pcfg.srf_typ];
mmil_anat_surf_plot(pcfg)


