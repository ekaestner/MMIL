
function ejk_surface_correlations(cfg)

%
%
%
%

%% Save Data
cor_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.cor_nme)
    cor_var.(cfg.cor_nme{iG}) = cfg.cor(:,iG);
end

grp_var.sbj_nme   = cfg.sbj_nme;
for iG = 1:numel(cfg.grp_nme)
    grp_var.(cfg.grp_nme{iG}) = cfg.grp(:,iG);
end

if ~isempty(cfg.cov)
    cov_var.sbj_nme   = cfg.sbj_nme;
    for iG = 1:numel(cfg.cov_nme)
        cov_var.(cfg.cov_nme{iG}) = cfg.cov(:,iG);
    end
end

ejk_chk_dir(cfg.out_dir)
out_hld = [ cfg.out_dir '/' cfg.out_pre];
ejk_chk_dir(out_hld)

%% Command
for iG = 1:numel(cfg.cor_nme)
    
    save( [ out_hld '/' 'cor_var.mat' ], 'cor_var' )
    save( [ out_hld '/' 'grp_var.mat' ], 'grp_var' )
    if ~isempty(cfg.cov); save( [ out_hld '/' 'cov_var.mat' ], 'cov_var' ); end

    for iCM = 1:numel(cfg.grp_cmp)
        
        fst_nme = cfg.grp_cmp{iCM};
        
        if ~isempty(cfg.cov)
            cov_loc = [ out_hld '/' 'cov_var.mat'];
        else
            cov_loc = '';
        end
        
        sve_cmd = sprintf('.libPaths(.libPaths()[2:5])\n');
        sve_cmd = sprintf('%ssource(''/home/ekaestner/gitrep/MMIL/ejk/R/R/ejk_surface_1way_ancova.r'')\n',sve_cmd);
        sve_cmd = sprintf('%sejk_surface_correlation( ''%s'', ''%s'', ''%s/cor_var.mat'', ''%s/grp_var.mat'', ''%s'', %i, ''%s'', ''%s'' )', sve_cmd, cfg.dta_lhs, cfg.dta_rhs, out_hld, out_hld, cov_loc, iG, fst_nme, out_hld);
        cell2csv( [out_hld '/example_R_script.r'], {sve_cmd} );
        unix( [ 'Rscript ' out_hld '/example_R_script.r' ] );
        
        cmd_pst = sprintf('\n\nlhs_srf_var_loc = ''%s''\n', cfg.dta_lhs);
        cmd_pst = sprintf('%srhs_srf_var_loc = ''%s''\n', cmd_pst, cfg.dta_rhs);
        cmd_pst = sprintf('%scor_loc = ''%s/cor_var.mat''\n', cmd_pst, out_hld);
        cmd_pst = sprintf('%sgrp_loc = ''%s/grp_var.mat''\n', cmd_pst, out_hld);
        cmd_pst = sprintf('%scov_loc = ''%s''\n', cmd_pst, cov_loc);
        cmd_pst = sprintf('%siG = %i\n', cmd_pst, iG );
        cmd_pst = sprintf('%sout_put_loc = ''%s''\n\n\n', cmd_pst, out_hld);
               % cmd_pst = sprintf('%fst_nme = ''%s''\n', cmd_pst, fst_nme );
        cmd_pst
        
        %% Make Plots
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
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rvl_dff_lhs(rvl_dff_lhs>0) = 0;
%         rvl_dff_rhs(rvl_dff_rhs>0) = 0;
        % rvalue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        low_rng_num = [ -0.40 0.40 ];
        hgh_rng_num = [ -0.90 0.90 ];
        
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
        
        low_rng_num = [ -0.30 0.30 ];
        hgh_rng_num = [ -0.60 0.60 ];
        
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = out_hld;
        pcfg.out_pre_fix = [ 'rvalues_diff' '_' fst_nme '_' 'p05'];
        
        pcfg.plt_dta = { rvl_dff_lhs_p05' rvl_dff_rhs_p05' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        % rvalue histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('Visible','off','Position',[ 0 0 1440 1080],'Units','pixels')
        subplot(2,1,1)
        hist(rvl_dff_lhs,1000)
        xlim([-1 1])
        title('R-Values LHS');
        subplot(2,1,2)
        hist(rvl_dff_rhs,1000)
        xlim([-1 1])
        title('R-Values RHS');
        tightfig();
        print([ out_hld '/' 'rvalues_diff' '_' fst_nme '_histogram' ],'-dpng')
        close all
        
        % pvalue histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('Visible','off','Position',[ 0 0 1440 1080],'Units','pixels')
        subplot(2,1,1)
        hist(pvl_lhs,1000)
        xlim([.0 .1])
        title('P-Values LHS');
        subplot(2,1,2)
        hist(pvl_rhs,1000)
        xlim([.0 .1])
        title('P-Values RHS');
        tightfig();
        print([ out_hld '/' 'pvalues_diff' '_' fst_nme '_histogram' ],'-dpng')
        close all
        
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
        
        rvl_dff_lhs_cls = rvl_dff_lhs;
        rvl_dff_rhs_cls = rvl_dff_rhs;
        
        rvl_dff_lhs_cls(lhs_bdd_ind) = 0;
        rvl_dff_rhs_cls(rhs_bdd_ind) = 0;
        
        %
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = out_hld;
        pcfg.out_pre_fix = [ 'rval' '_' fst_nme '_' 'cluster' '_pvl'  pvl_nme '_cls' cls_nme];
        
        pcfg.plt_dta = { rvl_dff_lhs_cls' rvl_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
        %%
        clear pvl_lhs pvl_rhs rvl_dff_lhs rvl_dff_rhs
        
    end
end

end