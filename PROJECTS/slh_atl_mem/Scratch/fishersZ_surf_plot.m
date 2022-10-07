out_hld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Analysis/FinalAnalysis/surf/fishersZ_256';

%%
pvl_lat = load([ out_hld '/' 'pvalues_ltle_atl_ltle_slah_lhs_dep_var.mat' ]); %  '_' 'sm' '313'
pvl_lat = pvl_lat.pvalues(2:end);

rvl_dff_lat = load([ out_hld '/' 'zvalues_ltle_atl_ltle_slah_lhs_dep_var.mat' ]); %  '_' 'sm' '313'
rvl_dff_lat = rvl_dff_lat.zvalues(2:end);

%%
low_rng_num = [ -1.65  1.65 ];
hgh_rng_num = [ -3 3 ];

%%
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 50;
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
pcfg.out_pre_fix = [ 'zval_ltle_cls50_pial' ];

pcfg.plt_dta = { rvl_dff_lat_cls' rvl_dff_lat_cls' };

pcfg.fmr_col_map = {'yellow green' 'green' 'dark green' 'grey' 'dark magenta' 'magenta' 'bright magenta'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

pcfg.srf_typ = 'inflated';
pcfg.out_pre_fix = [ 'zval_ltle_cls50_inflated' ];
mmil_anat_surf_plot(pcfg)



