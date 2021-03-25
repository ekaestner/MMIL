clear; clc;

jhn_dta_loc = '/home/ekaestne/PROJECTS/EXTERNAL/R_Output/';% '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_out = '/home/ekaestne/PROJECTS/OUTPUT/epd_age_phn';

low_rng_num = [ -0.01 nan ];
hgh_rng_num = [ -0.25 nan ];
pvl_chs = .05;
smt_stp = 176;
deg_fre = 79;
pvl_cls = .05;

%% Cluster threshold
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/rh.pial'); %fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/Misc/thickness-sphere-sm2819-lh.mgz'); % fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
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
cl_nme = num2str(round(cls_thr));


%% Quick Check Surface Map Plots
plt_nme = { 'TLE-MND_HCall_HC' 'TLE-MND_HC' 'TLE-MND_HC_wisconsin' };
smt_stp = { '176' '313' };

for iSM = 1:numel(smt_stp)
    for iPL = 1:numel(plt_nme)
        
        % Load
        pvl_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp{iSM} '/' 'par_new' '/' 'pValue' '_' 'lhs' '_' plt_nme{iPL} '.mat' ]);
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp{iSM} '/' 'par_new' '/'  'pValue' '_' 'rhs' '_' plt_nme{iPL} '.mat' ]);
        pvl_rhs = pvl_rhs.pValue;
        
        men_dff_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp{iSM} '/' 'par_new' '/' 'emmean' '_' 'lhs' '_' plt_nme{iPL} '.mat' ]);
        men_dff_lhs = men_dff_lhs.emmean;
        men_dff_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp{iSM} '/' 'par_new' '/' 'emmean' '_' 'rhs' '_' plt_nme{iPL} '.mat' ]);
        men_dff_rhs = men_dff_rhs.emmean;
        
        % Initial Plot
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_home'];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        %
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; %'/home/mmilmcd/data/FSRECONS/';
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
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'cluster_home' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end
    
    
    
    
    
    
    
    
    