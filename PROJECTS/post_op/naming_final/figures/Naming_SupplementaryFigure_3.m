load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

out_nme = { 'BNT' ...
            'ANT' };

out_hld = { '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/LTLE_post_post/bnt_raw_scr_pst' ...
            '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SurfaceCorrelation/LTLE_post_post/ant_mem_raw_scr_pst_novelpatients' };

fst_nme = 'tle_post_3T_ATLonly_left';

cfg.smt_stp = 313;
cfg.pvl_chs = .05;
cfg.pvl_cls = .05;

low_rng_num = [ -0.40 0.40 ];
hgh_rng_num = [ -0.70 0.70 ];

for iO = 1:numel(out_hld)
        
    % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rvl_dff_lhs = load([ out_hld{iO} '/' 'rvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
    rvl_dff_lhs = rvl_dff_lhs.rvalues;
    rvl_dff_rhs = load([ out_hld{iO} '/' 'rvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
    rvl_dff_rhs = rvl_dff_rhs.rvalues;
    
    % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pvl_lhs = load([ out_hld{iO} '/' 'pvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
    pvl_lhs = pvl_lhs.pvalues;
    pvl_rhs = load([ out_hld{iO} '/' 'pvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
    pvl_rhs = pvl_rhs.pvalues;
    
    %% Plot
    pcfg = [];
    
    pcfg.hme_wrk = 1;
    
    pcfg.out_dir     = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/figures/SupplementaryFigure_3';
    pcfg.out_pre_fix = [ out_nme{iO} '_' 'rval' '_' fst_nme ];
    
    pcfg.plt_dta = { rvl_dff_lhs' rvl_dff_rhs' };
    
    pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    
    pcfg.low_rng_num = low_rng_num;
    pcfg.hgh_rng_num = hgh_rng_num;
    
    mmil_anat_surf_plot(pcfg)
    
    %% Plot (p<.05)
    % rvalue, p<.05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rvl_dff_lhs_p05 = rvl_dff_lhs;
    rvl_dff_rhs_p05 = rvl_dff_rhs;
    
    rvl_dff_lhs_p05( pvl_lhs>.05) = 0;
    rvl_dff_rhs_p05( pvl_rhs>.05) = 0;
    
    pcfg = [];
    
    pcfg.hme_wrk = 1;
    
    pcfg.out_dir     = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/figures/SupplementaryFigure_3';
    pcfg.out_pre_fix = [ out_nme{iO} '_' 'rvalues_diff' '_' fst_nme '_' 'p05' '_' num2str(low_rng_num(2)) '_' num2str(hgh_rng_num(2))];
    
    pcfg.plt_dta = { rvl_dff_lhs_p05' rvl_dff_rhs_p05' };
    
    pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    
    pcfg.low_rng_num = low_rng_num;
    pcfg.hgh_rng_num = hgh_rng_num;
    
    mmil_anat_surf_plot(pcfg)
    
end