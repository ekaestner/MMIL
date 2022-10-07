%%
deg_fre = sum(strcmpi(grp_var.(cfg.grp_nme{iG}),fst_nme)) - 2;
        
% Cluster threshold
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/home/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/home/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

fcfg = [];

fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(cfg.smt_stp);
fcfg.df     = deg_fre;
fcfg.alpha  = cfg.pvl_cls;
fcfg.pval   = cfg.pvl_chs;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ...
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );
        
pvl_nme = num2str(roundsd(cfg.pvl_chs,3)); 
pvl_nme = pvl_nme(3:end);
cls_nme = num2str(round(cls_thr));

%%
pvl_lat = load([ out_hld '/' 'pvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
pvl_lat = pvl_lat.pvalues;

rvl_dff_lat = load([ out_hld '/' 'rvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
rvl_dff_lat = rvl_dff_lat.rvalues;

% pvl_lat = load([ out_hld '/' 'pvalues_ltle_atl_ltle_slah_lhs_dep_var.mat' ]); %  '_' 'sm' '313'
% pvl_lat = pvl_lat.pvalues; %pvl_lat(1) = [];
% 
% rvl_dff_lat = load([ out_hld '/' 'zvalues_ltle_atl_ltle_slah_lhs_dep_var.mat' ]); %  '_' 'sm' '313'
% rvl_dff_lat = rvl_dff_lat.zvalues; %rvl_dff_lat(1) = [];

%%
low_rng_num = [ -.30  .30 ];
hgh_rng_num = [ -.60  .60 ];

% low_rng_num = [ -1.6  1.6 ];
% hgh_rng_num = [ -3.25 3.25 ];

pcfg = [];

pcfg.hme_wrk = 1;
pcfg.srf_typ = 'pial';

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];

pcfg.plt_dta = { rvl_dff_lat' rvl_dff_lat' };


pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];
mmil_anat_surf_plot(pcfg)

%%
rvl_dff_lat_p05 = rvl_dff_lat;

rvl_dff_lat_p05( pvl_lat>.05) = 0;

pcfg = [];

pcfg.hme_wrk = 1;

pcfg.srf_typ = 'pial';

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p05' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];

pcfg.plt_dta = { rvl_dff_lat_p05' rvl_dff_lat_p05' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p05' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];
mmil_anat_surf_plot(pcfg)

%% 
rvl_dff_lat_p05 = rvl_dff_lat;

rvl_dff_lat_p05( pvl_lat>.01) = 0;

pcfg = [];

pcfg.hme_wrk = 1;

pcfg.srf_typ = 'pial';

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p01' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];

pcfg.plt_dta = { rvl_dff_lat_p05' rvl_dff_lat_p05' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p01' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) '_' pcfg.srf_typ ];
mmil_anat_surf_plot(pcfg)

%%
fcfg = [];
fcfg.pvl_thr = cfg.pvl_chs;
fcfg.cls_thr = 50;...cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lat_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lat );

%
lat_bdd_ind = find(pvl_lat_cls>fcfg.pvl_thr);

rvl_dff_lat_cls = rvl_dff_lat;

rvl_dff_lat_cls(lat_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 1;

pcfg.srf_typ = 'pial';

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rval' '_' fst_nme '_' 'cluster' '_cls' num2str(fcfg.cls_thr) '_' pcfg.srf_typ ];

pcfg.plt_dta = { rvl_dff_lat_cls' rvl_dff_lat_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'rval' '_' fst_nme '_' 'cluster' '_cls' num2str(fcfg.cls_thr) '_' pcfg.srf_typ ];
mmil_anat_surf_plot(pcfg)

%%
fcfg = [];
fcfg.pvl_thr = cfg.pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lat_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lat );

%
lat_bdd_ind = find(pvl_lat_cls>fcfg.pvl_thr);

rvl_dff_lat_cls = rvl_dff_lat;

rvl_dff_lat_cls(lat_bdd_ind) = 0;

%
pcfg = [];

pcfg.hme_wrk = 1;

pcfg.srf_typ = 'pial';

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rval' '_' fst_nme '_' 'cluster' '_cls' num2str(fcfg.cls_thr) '_' pcfg.srf_typ ];

pcfg.plt_dta = { rvl_dff_lat_cls' rvl_dff_lat_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'rval' '_' fst_nme '_' 'cluster' '_cls' num2str(fcfg.cls_thr) '_' pcfg.srf_typ ];
mmil_anat_surf_plot(pcfg)

