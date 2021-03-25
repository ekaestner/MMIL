%%
clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_final_sample_2.csv';
cov_nme = 'Epilepsy_Aging_Final_Sample_ASC_edit.csv';

epd_non_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
epd_usd_dir = '/home/mmilmcdRSI/data/fsurf';
adn_dir     = '/home/mmilmcdRSI/data_ADNI/fsurf/';

%% Constants
pvl_chs = .05;
smt_stp = 176;
deg_fre = 148;

pvl_cls = .25;

low_rng_num = [ -0.01 nan ];
hgh_rng_num = [ -0.25 nan ];

plt_out = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/recorticalsurfacemapq/';
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

%% Cluster threshold
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/home/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/home/thickness-sphere-sm2819-lh.mgz');
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
cls_nme = num2str(round(cls_thr));

%% Quick Check Surface Map Plots
plt_nme = { 'Bilateral_HC' 'Left_HC' 'Right_HC' };
smt_stp = { '' '176' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        pvl_lhs_176 = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/recorticalsurfacemapq/' '/' 'pValue_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        pvl_lhs_176 = pvl_lhs_176.pValue;
        pvl_rhs_176 = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/recorticalsurfacemapq/' '/' 'pValue_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        pvl_rhs_176 = pvl_rhs_176.pValue;
                
        % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        men_dff_lhs_176 = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/recorticalsurfacemapq/' '/' 'emmean_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs_176 = men_dff_lhs_176.emmean;
        men_dff_rhs_176 = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/recorticalsurfacemapq/' '/' 'emmean_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs_176 = men_dff_rhs_176.emmean;
        
        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMMEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        low_rng_num = [ -0.01 nan ];
        hgh_rng_num = [ -0.25 nan ];
        
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs_176' men_dff_rhs_176' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
                                
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
        fcfg.hms     = 'lh';
        pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs_176 );
        
        fcfg.hms     = 'rh';
        pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs_176 );
        
        % 
        lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
        rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);
        
        men_dff_lhs_cls = men_dff_lhs_176;
        men_dff_rhs_cls = men_dff_rhs_176;
        
        men_dff_lhs_cls(lhs_bdd_ind) = 0;
        men_dff_rhs_cls(rhs_bdd_ind) = 0;

        %
        pcfg = [];
        
        pcfg.hme_wrk = 1;
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp{iSM} '_' 'cluster' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end



