clear; clc;

slh_atl_mem_constants

load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

roi_ind = [ 2 6 7 8 9 10 11 13 14 15 16 17 19 20 21 23 24 25 26 27 31 32 34 36 ];

%% Cluster size
grp_use = grp.surgery.pst_cog_dti.ltle_atl;

deg_fre = numel(grp_use) - 2;

% Cluster threshold
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/home/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/home/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

% Find acceptable ROIs
[col_loc,albl,actbl] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/lh.aparc.annot');
    [albl num2cell(unique(col_loc))];
    use_ind = ismember( col_loc, roi_ind);

bad_ind = unique( [ bad_ind  ; find(~use_ind) ]);
gdd_ind = intersect( gdd_ind', find(use_ind) );

% Calculate
fcfg = [];

fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = sqrt(256);
fcfg.df     = deg_fre;
fcfg.alpha  = .10;
fcfg.pval   = .05;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
                                  fcfg.area , ...
                                  fcfg.fwhm , ...
                                  fcfg.df , ...
                                  fcfg.alpha , ...
                                  fcfg.pval );



