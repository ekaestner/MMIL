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

pvl_cls = .05;

low_rng_num = [ -0.01 nan ];
hgh_rng_num = [ -0.25 nan ];

%
deg_fre = [ 143 149 120 92 ];
smt_stp = [ 176 313 ];

%%
jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme     = { 'TLE_HC' 'MCI_HC' 'Early_HC' 'Late_HC' };
smt_stp_nme = { '176' '313' };

for iPL = 1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        plt_out = [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/figures/figure1/check/' smt_stp_nme{iSM} '_no_mass/' ];
        
        % Cluster threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
        srf_chr = fs_calc_triarea(srf_hld);
        
        ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
        bad_ind = find(ttt==0);
        gdd_ind = 1:numel(ttt);
        gdd_ind(bad_ind) = [];
        
        %
        fcfg = [];
        
        fcfg.nverts = numel(ttt)-numel(bad_ind);
        fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
        fcfg.fwhm   = 1.25 * sqrt(smt_stp(iSM));
        fcfg.df     = deg_fre(iPL);
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
        
        % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        pvl_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'pValue_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'pValue_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        pvl_rhs = pvl_rhs.pValue;
                
        % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        men_dff_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'emmean_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs = men_dff_lhs.emmean;
        men_dff_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'emmean_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs = men_dff_rhs.emmean;

        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMMEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
                
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr=0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr=0; end
        
        lhs_bdd_ind = find(pvl_lhs>pvl_lhs_fdr);
        rhs_bdd_ind = find(pvl_rhs>pvl_rhs_fdr);
        
        men_dff_lhs_fdr = men_dff_lhs;
        men_dff_rhs_fdr = men_dff_rhs;
        
        men_dff_lhs_fdr(lhs_bdd_ind) = 0;
        men_dff_rhs_fdr(rhs_bdd_ind) = 0;
          
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM} '_' 'fdr' ];
        
        pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
                
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
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
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM} '_' 'cluster' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end

%% More Liberal PValue 
pvl_chs = .05;
pvl_cls = .25;

jhn_dta_loc = '/home/jxrao/Desktop/Lab/surface_map/data/R_Output/';

plt_nme     = { 'TLE_HC' 'MCI_HC' 'Early_HC' 'Late_HC' };
smt_stp_nme = { '176' '313' };

for iPL = 3 %1:numel(plt_nme)
    for iSM = 1:numel(smt_stp)
        
        plt_out = [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/figures/figure1/check/' smt_stp_nme{iSM} '_no_mass_liberal/' ];
        ejk_chk_dir(plt_out);
        
        % Cluster threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        srf_hld = fs_read_surf('/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial');
        srf_chr = fs_calc_triarea(srf_hld);
        
        ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
        bad_ind = find(ttt==0);
        gdd_ind = 1:numel(ttt);
        gdd_ind(bad_ind) = [];
        
        %
        fcfg = [];
        
        fcfg.nverts = numel(ttt)-numel(bad_ind);
        fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
        fcfg.fwhm   = 1.25 * sqrt(smt_stp(iSM));
        fcfg.df     = deg_fre(iPL);
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
        
        % p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        pvl_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'pValue_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        pvl_lhs = pvl_lhs.pValue;
        pvl_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'pValue_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        pvl_rhs = pvl_rhs.pValue;
                
        % em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        men_dff_lhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'emmean_lhs' '_' plt_nme{iPL} '.mat' ]); %  '_' 'sm' '313'
        men_dff_lhs = men_dff_lhs.emmean;
        men_dff_rhs = load([ jhn_dta_loc '/' 'sm' smt_stp_nme{iSM} '/par/' '/' 'emmean_rhs' '_' plt_nme{iPL}  '.mat' ]); % '_' 'sm' '313'
        men_dff_rhs = men_dff_rhs.emmean;

        % plot threshold differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMMEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM}];
        
        pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
%         mmil_anat_surf_plot(pcfg)
                
        % pvalue-fdr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pvl_lhs_fdr = FDR(pvl_lhs,.05); if isempty(pvl_lhs_fdr); pvl_lhs_fdr=0; end
        pvl_rhs_fdr = FDR(pvl_rhs,.05); if isempty(pvl_rhs_fdr); pvl_rhs_fdr=0; end
        
        lhs_bdd_ind = find(pvl_lhs>pvl_lhs_fdr);
        rhs_bdd_ind = find(pvl_rhs>pvl_rhs_fdr);
        
        men_dff_lhs_fdr = men_dff_lhs;
        men_dff_rhs_fdr = men_dff_rhs;
        
        men_dff_lhs_fdr(lhs_bdd_ind) = 0;
        men_dff_rhs_fdr(rhs_bdd_ind) = 0;
          
        %
        pcfg = [];
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM} '_' 'fdr' ];
        
        pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };
        
        pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
%         mmil_anat_surf_plot(pcfg)
                
        % pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcfg = [];
        fcfg.pvl_thr = pvl_chs;
        fcfg.cls_thr = cls_thr; % mm^2
        fcfg.fsr_sbj = 'fsaverage';
        fcfg.fsr_dir = '/home/mmilmcd/data/FSRECONS/';
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
        
        pcfg.out_dir     = plt_out;
        pcfg.out_pre_fix = [ 'estmean_diff' '_' plt_nme{iPL} '_' 'sm' smt_stp_nme{iSM} '_' 'cluster' ];
        
        pcfg.plt_dta = { men_dff_lhs_cls' men_dff_rhs_cls' };
        
        pcfg.fmr_col_map = {'greenish grey' 'greenish grey' 'greenish grey' 'grey' 'dark red' 'red' 'bright yellow'};
        
        pcfg.low_rng_num = low_rng_num;
        pcfg.hgh_rng_num = hgh_rng_num;
        
        mmil_anat_surf_plot(pcfg)
        
    end
end
















































