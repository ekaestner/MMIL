
function ejk_surface_1way_ancova(cfg)

%
%
%
%

%% Save Data
grp_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.grp_nme)
    grp_var.(cfg.grp_nme{iG}) = cfg.grp(:,iG);
end

cov_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.cov_nme)
    cov_var.(cfg.cov_nme{iG}) = cfg.cov(:,iG);
end

ejk_chk_dir(cfg.out_dir)
save( [ cfg.out_dir '/' 'grp_var.mat' ], 'grp_var' )
save( [ cfg.out_dir '/' 'cov_var.mat' ], 'cov_var' )

%% Command
for iG = 1:numel(cfg.grp_nme)
    
    out_hld = [ cfg.out_dir '/' cfg.grp_nme{iG}];
    ejk_chk_dir(out_hld)
    
    for iCM = 1:numel(cfg.grp_cmp{iG})
        
        fst_nme = cfg.grp_cmp{iG}{iCM}{1};
        scd_nme = cfg.grp_cmp{iG}{iCM}{2};
        
        sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
        sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/R/ejk_surface_1way_ancova.r'')\n',sve_cmd);
        sve_cmd = sprintf('%sejk_surface_1way_ancova( ''%s'', ''%s'', ''%s/grp_var.mat'', ''%s/cov_var.mat'', %i, ''%s'', ''%s'', ''%s'' )', sve_cmd, cfg.dta_lhs, cfg.dta_rhs, cfg.out_dir, cfg.out_dir, iG, fst_nme, scd_nme, out_hld);
        cell2csv( [out_hld '/example_R_script.r'], {sve_cmd} );
        unix( [ 'Rscript ' out_hld '/example_R_script.r' ] );
        
        %% Make Plots
        deg_fre = sum(strcmpi(grp_var.(cfg.grp_nme{iG}),fst_nme)) + sum(strcmpi(grp_var.(cfg.grp_nme{iG}),scd_nme)) - numel(cfg.cov_nme) - 1;
        
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
        
        pvl_nme = num2str(cfg.pvl_chs); 
        pvl_nme = pvl_nme(3:end);
        cls_nme = num2str(round(cls_thr));
        
        % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs = load([ out_hld '/' 'pValue_lhs_dep_var' '_' fst_nme '_' scd_nme '.mat' ]); %  '_' 'sm' '313'
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ out_hld '/' 'pValue_rhs_dep_var' '_' fst_nme '_' scd_nme '.mat' ]); % '_' 'sm' '313'
        pvl_rhs = pvl_rhs.pValue;
        
        % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        men_dff_lhs = load([ out_hld '/' 'emmean_lhs_dep_var' '_' fst_nme '_' scd_nme '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs = men_dff_lhs.emmean;
        men_dff_rhs = load([ out_hld '/' 'emmean_rhs_dep_var' '_' fst_nme '_' scd_nme '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs = men_dff_rhs.emmean;
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMMEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        low_rng_num = [ -0.01 nan ];
        hgh_rng_num = [ -0.25 nan ];
        
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = out_hld;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' fst_nme '_' scd_nme];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        men_dff_lhs_cls = men_dff_lhs;
        men_dff_rhs_cls = men_dff_rhs;
        
        men_dff_lhs_cls(lhs_bdd_ind) = 0;
        men_dff_rhs_cls(rhs_bdd_ind) = 0;
        
        %
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = out_hld;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' fst_nme '_' scd_nme '_' 'cluster' '_'  pvl_nme '_' cls_nme];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end

end