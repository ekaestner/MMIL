clear; clc;

fsa_dir = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi';

out_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/FinalAnalysis/surface/xhemi';

fsaverage_sym_brain.surf_brain =  fs_read_surf([fsa_dir '/' 'lrrev.register.dat']);

%% Understand input
% Load
lhs_surf_brain.surf_brain =  fs_read_surf([fsa_dir '/' 'surf' '/' 'lh.pial']);
lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
lhs_srf_lbl = fs_read_label([fsa_dir '/' 'label' '/' 'lh.cortex.label']);

rhs_surf_brain.surf_brain =  fs_read_surf([fsa_dir '/' 'surf' '/' 'rh.pial']);
rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
rhs_srf_lbl = fs_read_label([ fsa_dir '/' 'label' '/' 'rh.cortex.label']);

% Load xhemi
lhs_surf_brain_xhm.surf_brain =  fs_read_surf([fsa_dir '/' 'xhemi' '/' 'surf' '/' 'lh.pial']);
lhs_surf_brain_xhm.surf_brain.coords = lhs_surf_brain_xhm.surf_brain.vertices;
lhs_srf_lbl_xhm = fs_read_label([fsa_dir '/' 'xhemi' '/' 'label' '/' 'lh.cortex.label']);

rhs_surf_brain_xhm.surf_brain =  fs_read_surf([fsa_dir '/' 'xhemi' '/' 'surf' '/' 'rh.pial']);
rhs_surf_brain_xhm.surf_brain.coords = rhs_surf_brain_xhm.surf_brain.vertices;
rhs_srf_lbl_xhm = fs_read_label([ fsa_dir '/' 'xhemi' '/' 'label' '/' 'rh.cortex.label']);

% Load sphere
lhs_surf_brain_xhm.surf_brain =  fs_read_surf([fsa_dir '/' 'xhemi' '/' 'surf' '/' 'lh.fsaverage_sym.sphere.reg']);
lhs_surf_brain_xhm.surf_brain.coords = lhs_surf_brain_xhm.surf_brain.vertices;
lhs_srf_lbl_xhm = fs_read_label([fsa_dir '/' 'xhemi' '/' 'label' '/' 'lh.cortex.label']);

% Load gwcsurf
lhs_surf_brain_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'FA_gwcsurf_gm-sphere-sm256-lh.mgz']));

% Create fake pattern
plt_val_lhs = zeros(1,lhs_surf_brain.surf_brain.nverts);
plt_val_lhs(1000:5000) = 0.5;
plt_val_lhs(50000:54000) = 1;
plt_val_lhs(90000:94000) = -0.5;
plt_val_lhs(120000:124000) = -1;

% Create ROI pattern
lhs_col_loc = fs_read_annotation([ fsa_dir '/' 'label' '/' 'lh.aparc.annot' ]);
rhs_col_loc = fs_read_annotation([ fsa_dir '/' 'label' '/' 'rh.aparc.annot' ]);

plt_roi_val_lhs = zeros(1,lhs_surf_brain.surf_brain.nverts);
plt_roi_val_rhs = zeros(1,rhs_surf_brain.surf_brain.nverts);

plt_roi_val_lhs(lhs_col_loc==8) = 0.5;
plt_roi_val_lhs(lhs_col_loc==9) = 1;
plt_roi_val_lhs(lhs_col_loc==16) = -0.5;
plt_roi_val_lhs(lhs_col_loc==25) = -1;

plt_roi_val_rhs(rhs_col_loc==8) = 1;
plt_roi_val_rhs(rhs_col_loc==9) = 0.5;
plt_roi_val_rhs(rhs_col_loc==16) = -1;
plt_roi_val_rhs(rhs_col_loc==25) = -0.5;

% Save ROI pattern
[ lhs_vol, lhs_mmm, lhs_prm, ~ ] = fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'FA_gwcsurf_gm-sphere-sm256-lh.mgz']);
[ rhs_vol, rhs_mmm, rhs_prm, ~ ] = fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'FA_gwcsurf_gm-sphere-sm256-rh.mgz']);

fs_save_mgh( plt_roi_val_lhs', [out_dir '/' 'epd080' '/' 'orig' '/' 'lhs_roi_dummy.mgz'], lhs_mmm, lhs_prm )
fs_save_mgh( plt_roi_val_rhs', [out_dir '/' 'epd080' '/' 'orig' '/' 'rhs_roi_dummy.mgz'], rhs_mmm, rhs_prm )


% Calculate Laterality Indices
plt_roi_lat_ind_lhs = zeros(1,lhs_surf_brain.surf_brain.nverts);
plt_roi_lat_ind_rhs = zeros(1,rhs_surf_brain.surf_brain.nverts);

plt_roi_lat_ind_lhs(lhs_srf_lbl) = ( plt_roi_val_lhs(lhs_srf_lbl)     - plt_roi_val_rhs(rhs_srf_lbl_xhm)) ./ ( plt_roi_val_lhs(lhs_srf_lbl)     + plt_roi_val_rhs(rhs_srf_lbl_xhm));
plt_roi_lat_ind_rhs(rhs_srf_lbl) = ( plt_roi_val_lhs(lhs_srf_lbl_xhm) - plt_roi_val_rhs(rhs_srf_lbl))     ./ ( plt_roi_val_lhs(lhs_srf_lbl_xhm) + plt_roi_val_rhs(rhs_srf_lbl));

%% FSAVERAGE
% 01 First plot
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '01_initial_hemisphere' ];

pcfg.plt_dta = { plt_val_lhs plt_val_lhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

% 02 xhemi plot
pcfg = [];

pcfg.fsr_dir = [ fsa_dir '/' 'xhemi' '/'];

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '02_xhemi_hemisphere' ];

pcfg.plt_dta = { plt_val_lhs plt_val_lhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

%% ROIs
% 03 First plot - ROI
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '03_initial_hemisphere_ROI' ];

pcfg.plt_dta = { plt_roi_val_lhs plt_roi_val_lhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

% sphere
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '03_initial_hemisphere_ROI_sphere' ];

pcfg.plt_dta = { plt_roi_val_lhs plt_roi_val_lhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

% 04 xhemi plot - ROI
pcfg = [];

pcfg.fsr_dir = [ fsa_dir '/' 'xhemi' '/'];

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '04_xhemi_hemisphere_ROI' ];

pcfg.plt_dta = { plt_roi_val_lhs plt_roi_val_lhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

%% LI ROIs
% First plot - Same ROIs
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '05_initial_hemisphere_same_ROI' ];

pcfg.plt_dta = { plt_roi_val_lhs plt_roi_val_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.4  .4];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

% xhemi plot - LI ROIs
pcfg = [];

pcfg.fsr_dir = [ fsa_dir '/' ];

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '06_xhemi_lateralityindice_same_ROI_org' ];

pcfg.plt_dta = { plt_roi_lat_ind_lhs plt_roi_lat_ind_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [  -.2  .2];
pcfg.hgh_rng_num = [ -1.0 1.0];

mmil_anat_surf_plot(pcfg)

%% GWCSURF
% First plot - gwcsurf 
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '07_initial_hemisphere_gwcsurf' ];

pcfg.plt_dta = { lhs_surf_brain_gwc lhs_surf_brain_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.15 .15 ];
pcfg.hgh_rng_num = [ -.40 .40 ];

mmil_anat_surf_plot(pcfg)

% xhemi plot - gwcsurf
pcfg = [];

pcfg.fsr_dir = [ fsa_dir '/' 'xhemi' '/'];

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '08_xhemi_hemisphere_gwcsurf' ];

pcfg.plt_dta = { lhs_surf_brain_gwc lhs_surf_brain_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ .15 .28 ];
pcfg.hgh_rng_num = [ .27 .40 ];

mmil_anat_surf_plot(pcfg)

%% Investigate epd080
lhs_srf_org_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'FA_gwcsurf_gm-sphere-sm256-lh.mgz']));
rhs_srf_org_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'FA_gwcsurf_gm-sphere-sm256-rh.mgz']));
lhs_srf_xhm_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'xhemi' '/' 'FA_gwcsurf_gm-lh_xhemi_test_sphere_sm256.mgh']));
rhs_srf_xhm_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'xhemi' '/' 'FA_gwcsurf_gm-rh_xhemi_test_sphere_sm256.mgh']));

% 09 Plot original
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '09_initial_right_right' ];

pcfg.plt_dta = { rhs_srf_org_gwc rhs_srf_org_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.25 .25 ];
pcfg.hgh_rng_num = [ -.40 .40 ];

mmil_anat_surf_plot(pcfg)

% 10 Plot right & right
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '10_initial_rightxhemi_right' ];

pcfg.plt_dta = { rhs_srf_xhm_gwc rhs_srf_org_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.25 .25 ];
pcfg.hgh_rng_num = [ -.40 .40 ];

mmil_anat_surf_plot(pcfg)

%% Investigate dummy ROIs
lhs_srf_org_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'lhs_roi_dummy.mgz']));
rhs_srf_org_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'orig' '/' 'rhs_roi_dummy.mgz']));
lhs_srf_xhm_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'xhemi' '/' 'lh_to_lh_roi_dummy.mgz']));
rhs_srf_xhm_gwc = mmil_rowvec(fs_load_mgh([out_dir '/' 'epd080' '/' 'xhemi' '/' 'rh_to_lh_roi_dummy.mgz']));

plt_roi_lat_ind_lhs = ( lhs_srf_org_gwc - rhs_srf_xhm_gwc ) ./ ( lhs_srf_org_gwc + rhs_srf_xhm_gwc );

% 11 Plot original
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '11_initial_right_right' ];

pcfg.plt_dta = { rhs_srf_org_gwc rhs_srf_org_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.40  .40 ];
pcfg.hgh_rng_num = [ -1.00 1.00 ];

mmil_anat_surf_plot(pcfg)

% 12 Plot right & right
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '12_initial_rightxhemi_right' ];

pcfg.plt_dta = { rhs_srf_xhm_gwc rhs_srf_org_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.40  .40 ];
pcfg.hgh_rng_num = [ -1.00 1.00 ];

mmil_anat_surf_plot(pcfg)

% 13 Plot LI
pcfg = [];

pcfg.fsr_dir = fsa_dir;

pcfg.out_dir     = out_dir;
pcfg.out_pre_fix = [ '13_LI' ];

pcfg.plt_dta = { plt_roi_lat_ind_lhs rhs_srf_org_gwc };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -.20  .20 ];
pcfg.hgh_rng_num = [ -1.00 1.00 ];

mmil_anat_surf_plot(pcfg)



