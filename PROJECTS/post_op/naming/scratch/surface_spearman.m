% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/bnt_raw_scr_pst_spearman';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/bnt_raw_scr_pst';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst_spearman';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/bnt_raw_scr_spearman';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/bnt_raw_scr';
% out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/ant_mem_raw_scr_spearman';
out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/ant_mem_bnt_raw_scr';

% fst_nme = 'tle_post_3T_ATLonly_left';
fst_nme = 'tle_controls_pre_3T_allSurg_all';

% p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_lhs = load([ out_hld '/' 'pvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
pvl_lhs = pvl_lhs.pvalues;
pvl_rhs = load([ out_hld '/' 'pvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
pvl_rhs = pvl_rhs.pvalues;

% em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvl_dff_lhs = load([ out_hld '/' 'rvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
rvl_dff_lhs = rvl_dff_lhs.rvalues;
rvl_dff_rhs = load([ out_hld '/' 'rvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
rvl_dff_rhs = rvl_dff_rhs.rvalues;

% rvalue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low_rng_num = [ -0.20 0.20 ];
hgh_rng_num = [ -0.50 0.50 ];

pcfg = [];

pcfg.hme_wrk = 1;

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2)) ];

pcfg.plt_dta = { rvl_dff_lhs' rvl_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)

% rvalue, p<.05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvl_dff_lhs_p05 = rvl_dff_lhs;
rvl_dff_rhs_p05 = rvl_dff_rhs;

rvl_dff_lhs_p05( pvl_lhs>.05) = 0;
rvl_dff_rhs_p05( pvl_rhs>.05) = 0;

pcfg = [];

pcfg.hme_wrk = 1;

pcfg.out_dir     = out_hld;
pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p05' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2))];

pcfg.plt_dta = { rvl_dff_lhs_p05' rvl_dff_rhs_p05' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = low_rng_num;
pcfg.hgh_rng_num = hgh_rng_num;

mmil_anat_surf_plot(pcfg)