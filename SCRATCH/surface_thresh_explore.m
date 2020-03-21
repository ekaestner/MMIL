clear; clc;

dta_loc = '/home/mmilmcdRSI/MetaData/MCD_Phenotypes/SurfGroupAvgs/output_n116/';
dta_fle = 'CortMD_diff_tval_HC_Language & memory-sm313-lh.mgh';

dta_out = '/home/ekaestne/PROJECTS/OUTPUT/mri_srf_cls/';
dta_lhs = 'CortMD_diff_tval_HC_LanguageMemory-sm313-lh.mgh';
dta_rhs = 'CortMD_diff_tval_HC_LanguageMemory-sm313-rh.mgh';

pvl_dir = '';
pvl_fle = '';

pvl_chs = .005;

%% Check cluster threshold
srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

%3549
fcfg = [];

fcfg.nverts = 163842-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(2819);
fcfg.df     = 156;
fcfg.alpha  = .05;
fcfg.pval   = pvl_chs;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ... 
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );

%% Test ANCOVA_TLE_vs_HC
plt_nme = 'ANCOVA_TLE_vs_HC';
pvl_nme = num2str(pvl_chs); pvl_nme = pvl_nme(3:end);
cls_nme = num2str(round(cls_thr));

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% threshold
fcfg = [];
fcfg.pvl_thr = pvl_chs;
fcfg.cls_thr = cls_thr/2; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% Plot
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p' pvl_nme '_' cls_nme 'mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%% Test ANCOVA_MCI_vs_HC
plt_nme = 'ANCOVA_MCI_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_MCI_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_MCI_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_MCI_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean*-1;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_MCI_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean*-1;

% threshold
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% Plot
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p05_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%% Test ANCOVA_TLE_vs_MCI
plt_nme = 'ANCOVA_TLE_vs_MCI';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_MCI.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_MCI.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_MCI.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_MCI.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% threshold
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

% Plot
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p05_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%% Test ANCOVA_EARLY_vs_HC
plt_nme = 'ANCOVA_EARLY_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Early_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Early_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Early_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Early_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% Plot Diffs
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot p<.05 & >100mm^2
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p05_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot p<.01 & >100mm^2
fcfg = [];
fcfg.pvl_thr = .01;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p01_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%% Test ANCOVA_LATE_vs_HC
plt_nme = 'ANCOVA_LATE_vs_HC';

% Load
pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_Late_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_Late_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_Late_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_Late_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;

% Plot Diffs
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_diff'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot p<.05 & >100mm^2
fcfg = [];
fcfg.pvl_thr = .05;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p05_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

% Plot p<.01 & >100mm^2
fcfg = [];
fcfg.pvl_thr = .01;
fcfg.cls_thr = 100; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

men_dff_lhs_cls = men_dff_lhs;
men_dff_rhs_cls = men_dff_rhs;

men_dff_lhs_cls(lhs_bdd_ind) = 0;
men_dff_rhs_cls(rhs_bdd_ind) = 0;

pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/cluster';
pcfg.out_pre_fix = [ plt_nme '_cluster_p01_100mm2'];

pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_thr = .05;
dof     = 46 + 17 - 1;

tvl_thr = abs(tinv(pvl_thr,dof));
tvl_hld = mmil_rowvec(fs_load_mgh([ dta_loc '/' dta_fle]));

tcdf(tvl_thr,dof)

%% surf cluster
g.fname_lhs = [ dta_out '/' dta_lhs];
g.fname_rhs = [ dta_out '/' dta_rhs];
g.p_thresh  = .05;
g.dof       = 46 + 17 - 1;
g.tails     = 2;
g.subj      = 'fsaverage';
g.hemi      = 'lh';
g.subjdir   = '/home/mmilmcd/data/FSRECONS/';
g.clust_thresh = 2819; % cluster size threshold (mm^2)

% Defaults
if ~isfield(g,'val_thresh'), g.val_thresh = 0; end
if ~isfield(g,'thresh_abs_flag'), g.thresh_abs_flag = 1; end

% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_lhs = g.fname_lhs;
fname_rhs = g.fname_rhs;
p_thresh = g.p_thresh;
dof = g.dof;
tails = g.tails;
val_thresh = g.val_thresh;
thresh_abs_flag = g.thresh_abs_flag;
clust_thresh = g.clust_thresh;
subj = g.subj;
hemi = g.hemi;
subjdir = g.subjdir;
clear g;

% LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tails==1, p_thresh = 2*p_thresh; end
t_thresh = abs(tinv(p_thresh,dof));
tst_lhs = mmil_rowvec(fs_load_mgh(fname_lhs,[],1)); % first frame only
tst_rhs = mmil_rowvec(fs_load_mgh(fname_rhs,[],1)); % first frame only

tmp_val_lhs = tst_lhs;
tmp_val_rhs = tst_rhs;
tmp_thresh = t_thresh;

% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lhs
tmp_fname_out_lhs = [ dta_out '/' 'test' '_' 'clust' '_' 'lhs.w' ];

cmd = 'mri_surfcluster';
cmd = sprintf('%s --in %s',cmd,fname_lhs);
cmd = sprintf('%s --thmin %0.6f',cmd,tmp_thresh);
cmd = sprintf('%s --subject %s --hemi %s',cmd,subj,'lh');
cmd = sprintf('%s --sd %s',cmd,subjdir);
cmd = sprintf('%s --o %s',cmd,tmp_fname_out_lhs);
cmd = sprintf('%s --minarea %0.6f',cmd,clust_thresh);
[status,result] = unix(cmd);

% rhs
tmp_fname_out_rhs = [ dta_out '/' 'test' '_' 'clust' '_' 'rhs.w' ];

cmd = 'mri_surfcluster';
cmd = sprintf('%s --in %s',cmd,fname_rhs);
cmd = sprintf('%s --thmin %0.6f',cmd,tmp_thresh);
cmd = sprintf('%s --subject %s --hemi %s',cmd,subj,'rh');
cmd = sprintf('%s --sd %s',cmd,subjdir);
cmd = sprintf('%s --o %s',cmd,tmp_fname_out_rhs);
cmd = sprintf('%s --minarea %0.6f',cmd,clust_thresh);
[status,result] = unix(cmd);

% MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lhs
[w,v] = fs_read_wfile(tmp_fname_out_lhs);
v = v(find(w)); % exclude 0 vals (shouldn't be there anyway)
mask = ones(size(tmp_val_lhs));
mask(v) = 1;

tmp_vals = zeros(size(tmp_val_lhs));
tmp_vals(v) = tmp_val_lhs(v);
tmp_val_lhs = tmp_vals;

tmp_val_lhs(abs(tmp_val_lhs)<abs(val_thresh))=0;

% rhs
[w,v] = fs_read_wfile(tmp_fname_out_rhs);
v = v(find(w)); % exclude 0 vals (shouldn't be there anyway)
mask = ones(size(tmp_val_rhs));
mask(v) = 1;

tmp_vals = zeros(size(tmp_val_rhs));
tmp_vals(v) = tmp_val_rhs(v);
tmp_val_rhs = tmp_vals;

tmp_val_rhs(abs(tmp_val_rhs)<abs(val_thresh))=0;

% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/mri_srf_cls//';
pcfg.out_pre_fix = [ 'Anny_HC_LM_cluster_threshold'];

pcfg.plt_dta = { tmp_val_lhs tmp_val_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%% fdr cluster



%% 