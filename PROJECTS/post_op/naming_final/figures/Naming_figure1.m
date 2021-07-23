load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

out_nme = { 'BNT' ...
            'ANT' };

out_hld = { '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/bnt_raw_scr' ...
            '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/ant_mem_raw_scr' };

fst_nme = 'tle_controls_pre_3T_allSurg_all';

cfg.smt_stp = 313;
cfg.pvl_chs = .05;
cfg.pvl_cls = .05;

low_rng_num = [ -0.15 0.15 ];
hgh_rng_num = [ -0.30 0.30 ];

for iO = 1:numel(out_hld)
    
    %% Correct
    % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pvl_lhs = load([ out_hld{iO} '/' 'pvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
    pvl_lhs = pvl_lhs.pvalues;
    pvl_rhs = load([ out_hld{iO} '/' 'pvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
    pvl_rhs = pvl_rhs.pvalues;
    
    % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rvl_dff_lhs = load([ out_hld{iO} '/' 'rvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
    rvl_dff_lhs = rvl_dff_lhs.rvalues;
    rvl_dff_rhs = load([ out_hld{iO} '/' 'rvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
    rvl_dff_rhs = rvl_dff_rhs.rvalues;

    % Cluster threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deg_fre = numel(grp.(fst_nme)) - 2;
    
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
    
    % Cluster correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fcfg = [];
    fcfg.pvl_thr = cfg.pvl_chs;
    fcfg.cls_thr = cls_thr; % mm^2
    fcfg.fsr_sbj = 'fsaverage';
    fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
    fcfg.hms     = 'lh';
    pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );
    
    fcfg.hms     = 'rh';
    pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );
    
    %
    lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
    rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);
    
    rvl_dff_lhs_cls = rvl_dff_lhs;
    rvl_dff_rhs_cls = rvl_dff_rhs;
    
    rvl_dff_lhs_cls(lhs_bdd_ind) = 0;
    rvl_dff_rhs_cls(rhs_bdd_ind) = 0;
    
    %% Plot
    pcfg = [];
    
    pcfg.hme_wrk = 1;
    
    pcfg.out_dir     = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/figures/Figure1';
    pcfg.out_pre_fix = [ out_nme{iO} '_' 'rval' '_' fst_nme '_' 'cluster' '_pvl'  pvl_nme '_cls' cls_nme '_test'];
    
    pcfg.plt_dta = { rvl_dff_lhs_cls' rvl_dff_rhs_cls' };
    
    pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    
    pcfg.low_rng_num = low_rng_num;
    pcfg.hgh_rng_num = hgh_rng_num;
    
    mmil_anat_surf_plot(pcfg)
    
end