clear; clc;

%% Initial Load
fsr_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/';

avg_sbj = 'fsaverage';
idv_sbj = 'FSURF_fc082_fmri_160923_20160923.152947_1'; 

% Groups
load( [ '/home/ekaestne/PROJECTS/OUTPUT' '/' 'PostOperative/Naming' '/' 'groups.mat' ] );
fst_nme = 'tle_controls_pre_3T_allSurg_all';
out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/ant_mem_raw_scr';

% Initial Parcellations
[ lhs_prc_loc, lhs_prc_lbl, ~]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'lh' '.aparc.annot' ] );
[ rhs_prc_loc, rhs_prc_lbl, ~]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'rh' '.aparc.annot' ] );

%% Find ROIs (fsaverage)
deg_fre = numel(grp.(fst_nme)) - 2;

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
fcfg.fwhm   = 1.25 * sqrt(313);
fcfg.df     = deg_fre;
fcfg.alpha  = .05;
fcfg.pval   = .05;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
    fcfg.area , ...
    fcfg.fwhm , ...
    fcfg.df , ...
    fcfg.alpha , ...
    fcfg.pval );

pvl_nme = num2str(roundsd(.05,3));
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

% pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.pvl_thr = .05;
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

% find indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_fus_roi = find(rvl_dff_lhs_cls>0);
rhs_new_fus_roi = find(rvl_dff_rhs_cls>0);

lhs_fus = find(lhs_prc_loc==8);
rhs_fus = find(rhs_prc_loc==8);

lhs_ovr_lap_fus = intersect(lhs_new_fus_roi,lhs_fus);
rhs_ovr_lap_fus = intersect(rhs_new_fus_roi,rhs_fus);

lhs_non_lap_fus = setxor(lhs_fus,lhs_ovr_lap_fus);
rhs_non_lap_fus = setxor(rhs_fus,rhs_ovr_lap_fus);

lhs_non_fus = setxor(lhs_new_fus_roi,lhs_ovr_lap_fus);
rhs_non_fus = setxor(rhs_new_fus_roi,rhs_ovr_lap_fus);

%% Create New ROIs (fsaverage)
kep_ind = [ 6 13 16 25 32 ];

[ lhs_prc_loc, lhs_prc_lbl, lhs_prc_ctb ]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'lh' '.aparc.annot' ] );
[ rhs_prc_loc, rhs_prc_lbl, rhs_prc_ctb ]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'rh' '.aparc.annot' ] );

% Setup Output structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_loc = zeros(size(lhs_prc_loc));
rhs_new_loc = zeros(size(rhs_prc_loc));

clear lhs_new_ctb rhs_new_ctb

% Save old ROIs for checking
for iK = 1:numel(kep_ind)
    lhs_new_loc( lhs_prc_loc==kep_ind(iK) ) = iK;
    rhs_new_loc( rhs_prc_loc==kep_ind(iK) ) = iK;
    
    lhs_new_ctb.struct_names{iK,1} = lhs_prc_ctb.struct_names{kep_ind(iK),1};
    rhs_new_ctb.struct_names{iK,1} = rhs_prc_ctb.struct_names{kep_ind(iK),1};
    
    lhs_new_ctb.table(iK,:) = lhs_prc_ctb.table(kep_ind(iK),:);
    rhs_new_ctb.table(iK,:) = rhs_prc_ctb.table(kep_ind(iK),:);
    
end

% New ROI - ?hs_ovr_lap_fus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_loc(lhs_ovr_lap_fus) = iK + 1;
rhs_new_loc(rhs_ovr_lap_fus) = iK + 1;

lhs_new_ctb.struct_names{iK+1,1} = 'new_fusiform';
rhs_new_ctb.struct_names{iK+1,1} = 'new_fusiform';

lhs_new_ctb.table(iK+1,:) = lhs_prc_ctb.table(end,:);
rhs_new_ctb.table(iK+1,:) = rhs_prc_ctb.table(end,:);

% New ROI - ?hs_non_lap_fus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_loc(lhs_non_lap_fus) = iK + 2;
rhs_new_loc(rhs_non_lap_fus) = iK + 2;

lhs_new_ctb.struct_names{iK+2,1} = 'non_overlapping_fusiform';
rhs_new_ctb.struct_names{iK+2,1} = 'non_overlapping_fusiform';

lhs_new_ctb.table(iK+2,:) = lhs_prc_ctb.table(end-1,:);
rhs_new_ctb.table(iK+2,:) = rhs_prc_ctb.table(end-1,:);

% New ROI - ?hs_non_fus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_loc(lhs_non_fus) = iK + 3;
rhs_new_loc(rhs_non_fus) = iK + 3;

lhs_new_ctb.struct_names{iK+3,1} = 'rest_of_new_roi';
rhs_new_ctb.struct_names{iK+3,1} = 'rest_of_new_roi';

lhs_new_ctb.table(iK+3,:) = lhs_prc_ctb.table(end-2,:);
rhs_new_ctb.table(iK+3,:) = rhs_prc_ctb.table(end-2,:);

% New ROI - everything else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load
[ lhs_old_loc, ~, lhs_old_ctb ]          = fs_read_annotation([ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'lh' '.aparc.annot' ]);
[ rhs_old_loc, ~, rhs_old_ctb ]          = fs_read_annotation([ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'rh' '.aparc.annot' ]);

% Label lhs
prc_ded_reg = find(lhs_new_loc==0);
org_ded_reg = find(lhs_old_loc==0);
prc_ded_reg = setxor(prc_ded_reg, org_ded_reg);

lhs_new_loc( prc_ded_reg ) = 9;
lhs_new_lbl{9} = 'deadzone';

lhs_new_ctb.numEntries = 9;
lhs_new_ctb.struct_names{9} = 'deadzone';
lhs_new_ctb.table = [ lhs_new_ctb.table ; lhs_old_ctb.table(20,:) ];

% Label rhs
prc_ded_reg = find(rhs_new_loc==0);
org_ded_reg = find(rhs_old_loc==0);
prc_ded_reg = setxor(prc_ded_reg, org_ded_reg);

rhs_new_loc( prc_ded_reg ) = 9;
rhs_new_lbl{9} = 'deadzone';

rhs_new_ctb.numEntries = 9;
rhs_new_ctb.struct_names{9} = 'deadzone';
rhs_new_ctb.table = [ rhs_new_ctb.table ; rhs_old_ctb.table(20,:) ];

% Ctab nonsense
lhs_sve_ctb = zeros(size(lhs_new_loc));
rhs_sve_ctb = zeros(size(rhs_new_loc));
for iV = 1:size(lhs_sve_ctb,1)
    
    if lhs_new_loc(iV)==0
        lhs_sve_ctb(iV)=0;
    else
        lhs_sve_ctb(iV) = lhs_new_ctb.table( lhs_new_loc(iV),1) + lhs_new_ctb.table( lhs_new_loc(iV),2)*2^8 + lhs_new_ctb.table( lhs_new_loc(iV),3)*2^16;
    end
    
    if rhs_new_loc(iV)==0
        rhs_sve_ctb(iV)=0;
    else
        rhs_sve_ctb(iV) = rhs_new_ctb.table( rhs_new_loc(iV),1) + rhs_new_ctb.table( rhs_new_loc(iV),2)*2^8 + rhs_new_ctb.table( rhs_new_loc(iV),3)*2^16;
    end
    
end

% Finish Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs_new_ctb.numEntries = numel(lhs_new_ctb.struct_names);
rhs_new_ctb.numEntries = numel(rhs_new_ctb.struct_names);

lhs_new_ctb.orig_tab = lhs_new_ctb.orig_tab;
rhs_new_ctb.orig_tab = rhs_new_ctb.orig_tab;

% Setup for saving
lhs_sve_ctb = zeros(size(lhs_new_loc));
rhs_sve_ctb = zeros(size(rhs_new_loc));
for iV = 1:size(lhs_sve_ctb,1)

    if lhs_new_loc(iV)==0
        lhs_sve_ctb(iV)=0;
    else
        lhs_sve_ctb(iV) = lhs_new_ctb.table( lhs_new_loc(iV),1) + lhs_new_ctb.table( lhs_new_loc(iV),2)*2^8 + lhs_new_ctb.table( lhs_new_loc(iV),3)*2^16;
    end
    
    if rhs_new_loc(iV)==0
        rhs_sve_ctb(iV)=0;
    else
        rhs_sve_ctb(iV) = rhs_new_ctb.table( rhs_new_loc(iV),1) + rhs_new_ctb.table( rhs_new_loc(iV),2)*2^8 + rhs_new_ctb.table( rhs_new_loc(iV),3)*2^16;
    end
    
end
    
% Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_write_annotation([ fsr_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus_cor' '.aparc.annot' ], ...
                     lhs_new_loc, lhs_sve_ctb, lhs_new_ctb)
fs_write_annotation([ fsr_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus_cor' '.aparc.annot' ], ...
                      rhs_new_loc, rhs_sve_ctb, rhs_new_ctb)
       
%% Figure of ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.fsr_dir = fsr_dir; 
fcfg.fsr_nme = avg_sbj;

fcfg.roi_loc = [ fsr_dir '/' avg_sbj '/' 'label' ];
fcfg.prc_nme = '.new_fus.aparc.annot';

fcfg.inc_reg = { 'new_fusiform' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = [ fsr_dir '/' 'roi_new_fus' '/'];
fcfg.out_nme  = 'roi_new_fus';

ejk_roi_plot(fcfg);
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.fsr_dir = fsr_dir; 
fcfg.fsr_nme = avg_sbj;

fcfg.roi_loc = [ fsr_dir '/' avg_sbj '/' 'label' ];
fcfg.prc_nme = '.new_fus_cor.aparc.annot';

fcfg.inc_reg = { 'new_fusiform'       'non_overlapping_fusiform' 'rest_of_new_roi' };
fcfg.roi_col = { rgb('light magenta') rgb('cerulean')            rgb('dark magenta') };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = [ fsr_dir '/' 'roi_new_fus_cor' '/'];
fcfg.out_nme  = 'roi_new_fus_cor';

ejk_roi_plot(fcfg);
                  
%% Resample (fsaverage to individual subject)
% LHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_parms = [];
tmp_parms.outdir      = [ fsr_dir '/' idv_sbj '/' 'label' '/' ];
tmp_parms.source_subj = 'fsaverage';
tmp_parms.subj        = idv_sbj;
tmp_parms.subjdir     = [ fsr_dir '/' ];
tmp_parms.verbose     = 1;
tmp_parms.forceflag   = 0;
args = mmil_parms2args(tmp_parms);
fname_out = fs_annot2annot( [ fsr_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus' '.aparc.annot' ], args{:});

% RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_parms = [];
tmp_parms.outdir      = [ fsr_dir '/' idv_sbj '/' 'label' '/' ];
tmp_parms.source_subj = 'fsaverage';
tmp_parms.subj        = idv_sbj;
tmp_parms.subjdir     = [ fsr_dir '/' ];
tmp_parms.verbose     = 1;
tmp_parms.forceflag   = 0;
args = mmil_parms2args(tmp_parms);
fname_out = fs_annot2annot( [ fsr_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus' '.aparc.annot' ], args{:});
   
%% Create Volume (individual subject)
fname_tmp = sprintf('%s/mri/%s+aseg.mgz', [ fsr_dir '/' idv_sbj '/' ], ['aparc.new_fus']);

cmd = [];
cmd = sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , fsr_dir);
cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj);
cmd = sprintf('%s   --annot %s', cmd, ['new_fus' '.aparc']);
cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ fsr_dir '/' idv_sbj], ['new_fus' '.aparc'] );
cmd = sprintf('%s   --volmask',cmd);
cmd = sprintf('%s   --o %s\n',cmd,fname_tmp);
[s,r] = unix(cmd);

%% Create wmparc (individual subject)
fname_out = sprintf('%s/%s.mgz',[ fsr_dir '/' idv_sbj '/' 'mri' '/' ],['wmparc' '_' 'new_fus']);

cmd = [];
cmd = sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , fsr_dir);
cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj );
cmd = sprintf('%s   --annot %s', cmd, ['new_fus' '.aparc'] );
cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ fsr_dir '/' idv_sbj], ['new_fus' '.aparc']);
cmd = sprintf('%s   --wmparc-dmax 5', cmd);
cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown', cmd);
cmd = sprintf('%s   --volmask', cmd);
cmd = sprintf('%s   --ctxseg %s+aseg.mgz', cmd, ['aparc.new_fus']);
cmd = sprintf('%s   --o %s', cmd, fname_out);
[s,r] = unix(cmd);

%% Extract new wmparc (individual subject)



