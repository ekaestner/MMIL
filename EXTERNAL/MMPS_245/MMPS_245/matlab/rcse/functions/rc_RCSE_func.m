function [best_retmap,min_error]=rc_RCSE_func(avg_data,subjname,varargin)
%function [best_retmap,min_error]=rc_RCSE_func(avg_data,subjname,[options])
%
% Purpose: retinotopy-constrained source estimation
%
% Usage:
%  rc_RCSE_func(avg_data,subjname,'key1',value1,...);
%  [best_retmap,min_error]=rc_RCSE_func(avg_data,subjname,'key1',value1,...);
%
% Required Input:
%  avg_data: average data structure
%    (see ts_process_fif_data and ts_avg_fif_data)
%  subjname: freesurfer subject name (must be found in $SUBJECTS_DIR)
%
% Optional parameters of which one or the other is required:
%  'retmap': struct specifying vertices for visual areas and additional dipoles
%    visual area vertices can be defined explicitly with vertex numbers
%    or with wfiles (specifying a cluster of vertices)
%    See examp_define_retmap.m, examp_define_retmap_wfiles.m
%      or use rc_define_retmap_from_csv.m or rc_define_retmap_from_retfit.m
%     {default = []}
%  'retfit_results': struct containing the following fields:
%    lh_data, rh_data, lh_fit_data, lh_area_masks, rh_fit_data, rh_area_masks
%    These are the outputs fo retfit and define the fitted retinotopy data
%      and ROIs for V1, V2, and V3
%    If 'retfit_results' is supplied, input 'retmap' is ignored
%    Also, 'lh_verts' and 'rh_verts' will be derived from vertices in retfit
%    See rc_load_retfit_results
%    {default = []}
%
% Optional parameters:
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'prefix': prefix of all output files
%     {default = 'RCSE'}
%  'rootoutdir': root output directory (subdirectories will be created)
%    {default = pwd} (current working directory)
%  'verbose': [0|1] display status messages
%     {default = 1}
%  'forceflag': [0|1|2] whether to overwrite files
%     0: do nothing (except maybe plots) if output exists
%     1: overwrite results and other output
%     2: overwrite all output, including forward solution
%     {default = 0}
%  'use_areas': vector of area numbers defining subset of visual areas
%     defined in retmap to use
%     {default = []} (if empty, use all areas in retmap)
%  'conditions': vector of condition numbers (index to avg_data.averages)
%     defining subset of stimulus locations in retmap to use
%     If empty, use all conditions in retmap.cond_order
%       or conditions in cond_info with non-zero contrast
%     {default = []}
%  'event_codes': vector of event codes to use
%     defining subset of stimulus conditions in cond_info to use
%     If empty, use all conditions in cond_info with non-zero contrast
%     If specified, 'conditions' will be ignored
%     {default = []}
%  'wfiledir': directory (relative to running dir) containing FreeSurfer wfiles
%     used if retmap is defined using wfiles
%     {default = './weights'}
%  'calc_dipinfo_flag': [0|1] whether to calculate dipole info from surfaces
%     otherwise load from .dip files
%    {default = 0}
%  'lh_dip_file': name of left hemi freesurfer dipole file (from tksurfer)
%     {default = 'bem/lh_white.dip'} (relative to subject dir)
%  'rh_dip_file': name of right hemi freesurfer dipole file (from tksurfer)
%     {default = 'bem/rh_white.dip'} (relative to subject dir)
%  'lh_dip_info': matrix of left hemi dipole information (6 x ndips)
%     6 columns for x,y,z coordinates and nx,ny,nz normal vector
%     if specified, will ignore lh_dip_file
%     {default = []}
%  'rh_dip_info': matrix of right hemi dipole information (6 x ndips)
%     6 columns for x,y,z coordinates and nx,ny,nz normal vector
%     if specified, will ignore rh_dip_file
%     {default = []}
%  'lh_verts': vector of left hemisphere vertex numbers to use in
%     raw forward matrix (be sure these vertex numbers are 1-based)
%     if empty, will use only the vertices specified in retmap
%     {default = []}
%  'rh_verts': vector of right hemisphere vertex numbers to use in
%     raw forward matrix (be sure these vertex numbers are 1-based)
%     if empty, will use only the vertices specified in retmap
%     {default = []}
%  'baseline_start: start time of noise period (msec)
%     if ncov_type = 1, use this time window for noise covariance
%     and scaling factor calculation from average data
%     relative to trigger onset; negative times occur before trigger
%     { default = -Inf }
%  'inverse_type': [0|1|2] type of inverse calculations
%    0: unregularized pseudo-inverse (fast, very little memory)
%    1: regularized psuedo-inverse with identity matrix for noise covariance
%       no source covariance matrix (fast, very little memory)
%    2: regularized psuedo-inverse with noise covariance matrix from data
%       depending on ncov_type, uses source covariance matrix
%       if indy_locs_flag=1, inverse_type will be set to 2
%    {default = 1}
%  'ncov_type': specify what type of noise covariance matrix to use
%     if 0, use identity matrix
%       (assume uniform white noise, independently scaled for each sensor type)
%     if 1, calculate and use noise covariance matrix from average prestim
%     if 2, use noise covariance matrix calculated from single trials
%       (stored in avg_data.noise.covar)
%    {default = 1}
%  'raw_ncov_lambda': regularization factor used if ncov_type = 2
%    {default = 0.1}
%  'SNR': estimated signal-to-noise-ratio (for regularization parameter)
%      {default = 10}
%  'baseline_end': end time of noise period (msec)
%     { default = 0 }
%  'calc_scalefacts_flag': [1|0] Toggle calculation of scaling factors
%     for different channel types (standard deviation of avg baseline)
%     if 1, calculate these scaling factors and apply to data and forward matrix
%     if 0, use default or user specified scaling factors
%     {default = 0}
%  'grad_scalefact': scaling factor applied to gradiometer data
%     purpose of scaling factors is to get data from different channel types
%     into roughly the same scale
%     {default = 10^13}
%  'mag_scalefact': scaling factor applied to magnitometer data
%     {default = 10^15}
%  'EEG_scalefact': scaling factor applied to EEG data
%     {default = 10^6}
%  'bem_flag': [1|0] use boundary element method (BEM)
%      if 0, use spherical shell model instead
%      {default = 1}
%  'openmeeg_flag': [0|1] use OpenMEEG to calculate BEM solution
%    otherwise use M.X. Huang's BEM solution function
%    {default = 0}
%  'conductivities': vector of three values specifying conductivity of:
%      (1) brain (2) skull  (3) scalp (S/m)
%      {default = [0.3 0.01 0.3]}
%  'EEG_gain_scalefact': conductivity scaling factor
%     an optimized value is necessary for integration of MEG and EEG data
%     {default = 1}
%  'nlayers': number of layers for BEM or spherical shell model
%      For MEG only, can be 1 or 3; with EEG (alone or with MEG), must be 3
%      {default = 3}
%  'badchans': vector of bad channel indices
%      { default = [] }
%  'badchanfile': name of text file containing bad channel labels
%     {default = []}
%  'usegrad_flag': [1|0] Toggle use of gradiometer data, if available
%     {default = 1}
%  'usemag_flag': [1|0] Toggle use of magiometer data, if available
%     {default = 0}
%  'useEEG_flag': [1|0] Toggle use of EEG data, if available
%     {default = 0}
%  'refEEG_coords': reference EEG electrode coordinates (e.g. [0.004,0.020,0.001])
%     if specified, forward solution for these coordinates is subtracted 
%       from forward solution for all other EEG electrodes
%     for data originally in fif format, these coordinates are already
%       stored in the "loc" matrix for each EEG electrode
%     {default = []}
%  'visual_area_weight': source covariance matrix weighting for visual
%     area dipoles (prior estimate of relative strength)
%     {default = 1} (can range from 0-1)
%  'ret_dips_weight': weighting for extra retinotopic dipoles
%     {default = 1} (can range from 0-1)
%  'nonret_dips_weight': weighting for extra non-retinotopic dipoles
%     {default = 1} (can range from 0-1)
%  'ret_dips_lh_dec_file': name of left hemi FreeSurfer dec file
%     (from tksurfer) for extra retinotopic dipoles --
%     these can be instead of or in addition to extra dipoles defined in retmap
%     dec file contains 0's and 1's indicating subset of dipoles to use
%     {default = []} (relative to subject dir)
%  'ret_dips_rh_dec_file': name of right hemi dec file for retinotopic dipoles
%     {default = []} (relative to subject dir)
%  'nonret_dips_lh_dec_file': name of left hemi FreeSurfer dec file (from 
%     tksurfer) for extra non-retinotopic dipoles -- these can be instead of or
%     in addition to extra dipoles defined in retmap
%     {default = []} (relative to subject dir)
%  'nonret_dips_rh_dec_file': name of right hemi dec file for non-ret dipoles
%     {default = []} (relative to subject dir)
%  'ret_dip_qfield_flag': [1|0] Toggle quarterfield constraint
%     if 1, extra retinotopic dipoles are constrained to have the
%     same source strength and orientation within a visual quarter field
%     (reduces number of unknowns and prevents these dipoles from absorbing too
%      much variance)
%     {default = 1}
%  'norm_weights_flag': [0|1|2|3] whether and how to normalize vertex weights
%     0: do not normalize
%     1: normalize so sum of weights for each stimulus location equals 1
%     2: normalize so average sum of weights (within a visual area) equals 1
%     3: multiply forward matrix by vertex areas (changes source units to nA m/mm^2)
%     {default = 2}
%  'surfname': surface name to be loaded to
%     calculate area if norm_weights_flag = 3
%     {default = 'white'}
%  'w_thresh_patch' - threshold applied to cortical patch weights
%       relative to max val for each patch
%     this reduces the number of dipoles that included in the forward matrix
%       cutting off tails of smooth, gaussian-like cluster
%     {default: 0}
%  'near_nbr_weight': if non-zero, sum visual area dipoles across nearest
%     neighbors and apply this weighting to the nearest neighbors (central
%     dipole will get weight=1)
%     {default = 0}
%  'indy_locs_flag': [1|0] Toggle calculation of independent source estimates
%     for each stimulus location
%     {default = 0}
%  'indy_smfact': smoothness factor for independent location estimates
%     (use of neighbor covariance matrix) related to proximity of polar angle
%     Only applicable when indy_locs_flag=1
%     {default = 0.999}
%  'ecc_smfact': smoothness factor for independent location estimates
%     (neighbor covariance) related to proximity of eccentricity
%     Only applicable when indy_locs_flag=1
%     {default = 0.999} (can range from 0-1)
%  'upperlower_smfact': smoothness factor across upper and lower visual fields.
%     e.g. smfact=0 for independent sources
%          smfact=1 for identical sources
%    {default = 1}
%  'hemi_smfact': smoothness factor across left and right hemifields
%     e.g. smfact=0 for independent sources
%          smfact=1 for identical sources
%    {default = 1}
%  'loose_flag': [1|0] toggle loose orientation constraint
%     Automatically sets indy_locs_flag=1
%     {default = 0}
%  'loose_tang_weight': source covariance weight of tangential components
%     value of 1 means equal weighting with normal vector
%     Only applicable when loose_flag=1
%     {default = 0.1}
%  'condF_thresh': exclude stimulus locations for which cond(F'*F)
%     is greater than this threshold
%     {default = 0}
%  'write_err_flag': [1|0] Whether to save residual error sensor waveforms
%     {default = 0}
%  'write_fit_flag': [1|0] Whether to save fitted sensor waveforms
%     {default = 0}
%  'write_areafit_flag': [1|0] Whether to save fitted sensor waveforms
%     separately for each visual area
%     {default = 0}
%  'write_fif_flag': [1|0] Whether to save fit and residual error
%     sensor waveforms as fif files (displayable by Neuromag's xplotter)
%     {default = 0}
%  'forward_matfile': mat file containing gain matrix
%     This allows user to supply precalculated gain matrix from an earlier
%       run of this program -- for example, calculate gain matrix with
%       all channels (grad, mag, and EEG) and then rerun excluding certain
%       channels or channel types
%     The gain matrix contained in this file is not retinotopy constrained
%       but does only contain the dipoles specified in retmap
%       or lh_verts and rh_verts
%     Do not use this option if you change the number of areas or dipoles
%       or locations used (unless lh_verts and rh_verts include all possible dipoles)
%     {default = []}
%
% Optional parameters for defining retmap from retfit_results
%  'vf2ctx_flag': [0|1] precompute mapping from visual field to cortex
%    {default = 1}
%  'stim_type': type of stimulus to model
%    0: point: generation of retmap is very fast, but less accurate -- correction possible?
%    1: circle: slower, but more accurately takes into account receptive field size
%    2: wedge: slower, but more accurately takes into account receptive field size
%    {default = 2}
%  'ecc_width': eccentricity width of stimuli (deg. vis. ang.)
%    will be ignored if fname_conds has ecc_width column
%    {default = 1}
%  'theta_width': polar angle width of stimuli (degrees)
%    will be ignored if fname_conds has theta_width column or if stim_type=1
%    {default = 10}
%  'w_thresh': threshold applied to weights (relative to max)
%    {default = 0.01}
%  'vfnorm_flag': [0|1] for applying w_thresh for vf2ctx
%     0: normalize to global max
%     1: normalize to max for each cortical location
%     {default = 1}
%  'single_vertex_flag': [0|1] select one vertex for each stimulus location
%    0: use all vertices
%    1: vertex with maximum weight
%    2: vertex closest to center of mass
%    {default = 0}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'rf_sizes': vector of receptive field sizes for each visual area
%    Will be used as initial estimate if rf_niters>0
%    {default = [1,1,1]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Will be used as initial estimate if rf_niters>0
%    {default = [0.1,0.1,0.1]}
%  'rf_r_inter': radius value used as intercept at which rf size is rf_sizes
%    {default = 6}
%  'restrict_hemi_flag': [0|1] whether to restrict dipole clusters
%    so no ipsilateral
%    {default = 0}
%  'restrict_uplow_flag': [0|1] whether to restrict dipole clusters
%    so no upper-lower cross-over
%    {default = 0}
%  'retfit_data_flag': [0|1] whether to use refit area ROIs but original
%    retinotopy data instead of template values for selection of dipole clusters
%    {default = 0}
%
% Optional parameters for non-linear fitting
%  'reweight_flag': [0|1] whether to use iteratively re-weighted least squares
%     {default = 0}
%  'reweight_leverage_flag': [0|1] for IRLS, use "leverage"
%     {default = 1}
%  'reweight_leverage_max_flag': [0|1] for IRLS, take max leverage across sensors
%     otherwise take mean
%     {default = 1}
%  'hybrid_flag': [0|1] whether to use hybrid correlation offset search
%     {default = 0}
%  'offset_niters': number of iterations for random offset search to fit data
%     This uses retfit with different ecc or theta offsets
%     to see if they give a better fit to the data
%     If 0, do exhaustive search of offsets
%       specified by all combinations of r_offset and th_offset vectors
%     If not 0, retfit_results must be supplied
%     {default = 0}
%  'offset_mstarts': number of multi-start iterations
%       for random offset search (ignored if offset_niters==0)
%     If 0, use 0 offsets for starting point
%     If >=1, randomly choose offsets for starting point
%     {default = 0}
%  'prior_prefix': prefix of output files from a previous run of RCSE
%     If not full path, assumed to be relative to pwd/matfiles
%     The prior source waveforms will be used as a reference.
%     For offset or nbrhd search, cost function will be correlation with reference
%       between corr_time0 and corr_time1.
%     Note: only applicable if offset_niters>0 or nbrhd_niters>0
%     {default = []}
%  'prior_weight': weighting applied to prior relative to residual error
%     {default = 0.5}
%  'prior_avg_flag': [0|1|2] whether to use averaged prior
%     0: no averaging -- requires that prior and fitted waveform have
%          same number of visual areas and contrast levels
%     1: average prior across contrast levels, but not visual areas
%     2: average prior across visual areas, but not contrast levels
%     3: average prior across contrsat levels and visual areas
%     {default = 0}
%  'prior_mindiff_flag': [0|1|2] whether to minimize absolute difference
%     0: maximize correlation
%     1: minimize absolute difference
%     2: minimize absolute difference after optimally scaling
%     {default = 0}
%  'prior_zscore_flag': [0|1] if prior contains S_sem, use to calculate
%     z-score for prior-constrained optimization cost function
%     only applies if prior_mindiff_flag = 1 or 2
%     {default = 0}
%  'prior_sem_min': set all S_sem values to at least this value
%     prevent very small SEM values from magnifying small errors
%     {default = 0.5}
%  'prior_indy_wform_flag': [0|1] for prior, scale prior waveform
%     indepdently for each area and contrast
%     {default = 0}
%  'corr_time0': start of time range (msec) to calculate
%     correlation or difference with reference
%     Note: only applicable if prior_prefix not empty
%     {default = 0}
%  'corr_time1': end of time range (msec) to calculate
%     correlation or difference with reference
%     Note: only applicable if prior_prefix not empty
%     {default = 170}
%  'fit_range_flag': [0|1] for nonlinear fits, minimize error between
%     fit_time0 and fit_time1 (otherwise use entire time range)
%     {default = 0}
%  'fit_time0': start of time range (msec) to fit data
%     {default = 60}
%  'fit_time1': end of time range (msec) to fit data
%     {default = 80}
%  'polarity_penalty': cost multiplier to penalize waveforms with positive
%     or negative polarity
%     If positive, source waveforms with positive mean amplitude within
%       fit time window will be penalized according to this formula:
%         total_error = total_error*(1 + mean_amp*polarity_penalty)
%     If negative, source waveforms with negative mean amplitude within
%       fit time window will be penalized
%     If zero, no polarity penalty
%     {default = 0}
%  'mstart_rotsearch_flag': [0|1] whether to do rotation search for each
%     offset search multistart iteration
%     If 0, do rotation search (if rot_niters>0) once at end
%     If 1, do rotation search (if rot_niters>0) at end of each multistart
%     {default = 0}
%  'mstart_rfsearch_flag': [0|1] whether to do rf sizes search for each
%     offset search multistart iteration
%     If 0, do rf sizes search (if rf_niters>0) once at end
%     If 1, do rf sizes search (if rf_niters>0) at end of each multistart
%     {default = 0}
%  'mstart_randstart_flag': [0|1] whether to use random offsets to start
%     for multistart iteration (except first)
%     otherwise, use 0 offsets for each multistart iteration
%     {default = 1}
%  'mstart_average_flag': [0|1] whether to average waveforms across
%     multi-start iterations and calculation standard deviation
%     {default = 1}
%  'offset_niters_last': number of iterations for extra multi-start iteration of
%     random offset search using best offsets as starting point
%     Ignored if offset_niters==0 or offset_mstarts==0
%     {default = 0}
%  'offset_group_flag': [0|1|2|3|4] force conditions to have the same offsets
%     0: each stimulus locations can have a different offset
%     1: all within an eccentricity band and quadrant must have same offset
%     2: all within a visual field quadrant must have same offset
%     3: all within a hemifield must have same offset
%     4: all stimulus locations must have same offset
%     {default = 1}
%  'offset_group_patches_flag': [0|1] force patches for a stimulus location
%     to have the same offsets (i.e. ipsi, upper-lower cross-over)
%     {default = 0}
%  'offset_group_areas_flag': [0|1] force areas to have the same offsets
%     {default = 0}
%  'offset_const_areas': vector of area indices for which offsets should be
%     held constant (e.g. [1,2] to keep V1 and V2 constant and fit only for V3)
%     {default = []}
%  'r_step': standard deviation of gaussian step size for rand offset search
%    {default = 0.1}
%  'th_step': standard deviation of gaussian step size for rand offset search
%    {default = 2}
%  'r_offset_range': vector of min and max r_offsets for rand offset search
%    {default = [-3,3]}
%  'th_offset_range': vector of min and max th_offsets for rand offset search
%    {default = [-20,20]}
%  'grid_offset_flag': instead of r and th, make offets to grid
%    coordinates u and v (unit grid)
%    {default = 0}
%  'grid_offset_smooth': number of smoothing steps for expanding number of
%     vertices included in forward solution (only applies if grid_offset_flag=1)
%     value corresponds roughly to number of vertices by which to expand border
%    {default = 0}
%  'rot_niters': number of iterations for random rotation search to fit data
%     This tests different registrations between MEG and MRI
%     to see if they give a better fit to the data
%     {default = 0}
%  'rot_step': standard deviation of rotation step size (degrees)
%     {default = 1}
%  'rot_max': maximum rotation (degrees)
%     If zero, allow unconstrained random search for best rotation
%     {default = 10}
%  'best_rot': best fitting rotation (degrees) relative to mri2head 'trans'
%     vector of [x y z] rotations
%     {default = [0 0 0]}
%  'rscale_niters': number of iterations for eccentricity (r) scaling
%     This adjusts the mapping between MEG and MRI stimuli (scales r_max)
%     {default = 0}
%  'rscale_step': standard deviation of rscale step size
%     {default = 0.01}
%  'rscale_range': vector of min and max rscale
%     {default = [0.8,1.2]}
%  'best_rscale': best fitting scaling factor for r_max
%     {default = [1s]}
%  'nbrhd_niters': number of iterations for random dipole-search to fit data
%     This tests neighboring vertices, defined by retmap from retfit or w files,
%     to see if they give a better fit to the data
%     {default = 0}
%  'nbrhd_full_flag': [0|1] for nbrhd search, test one stimulus location
%     at a time and do exhaustive search of all neighboring vertices
%     {default = 0}
%  'best_retmap': "best-fit" retmap with subset of dipoles from previous fit
%     {default = []}
%  'min_error': minimum error from previous fits
%     {default = []}
%  'r_offset': eccentricity offset for generating retmap from retfit
%     If a vector, will test each r_offset specified for best fit
%     Ignored if offset_niters>0
%     {default = 0}
%  'th_offset': polar angle offset for generating retmap from retfit
%     If a vector, will test each th_offset specified for best fit
%     Ignored if offset_niters>0
%     {default = 0}
%  'rf_niters': number of iterations for random search for best-fitting
%     rreceptive field sizes for each visual area
%     {default = 0}
%  'rf_min_sizes': vector of minimum receptive field sizes (deg vis angle)
%     one value for each visual area
%     {default = [0.1,0.1,0.1]}
%  'rf_max_sizes': vector of minimum receptive field sizes (deg vis angle)
%     one value for each visual area
%     {default = [3,3,3]}
%  'rf_min_slopes': vector of minimum values for fitting rf_slopes
%     {default = [0 0 0]}
%  'rf_max_slopes': vector of maximum values for fitting rf_slopes
%     {default = [0.2 0.2 0.2]}
%  'rf_size_step': standard deviation of rf size step
%     {default = 0.1}
%  'rf_slope_step': standard deviation of rf slope step
%     {default = 0.01}
%  'best_rf_sizes': vector of best fitting receptive field sizes
%     {default = []}
%  'best_rf_slopes': vector of best fitting rf ecc slopes
%     {default = []}
%
% Optional parameters that may be required (depending on other parameters):
%  'cond_info': struct array containing condition information
%    see rc_read_cond_info
%    Required if retfit_results is supplied
%    {default = []}
%  'template_fif': full path name of online average fif file
%     Required if write_fif_flag=1
%     {default = []}
%  'bem_surf_files': cell array of file names for FreeSurfer tri files
%    (text files specifying triangles and faces) with MRI-derived skull
%     and scalp surfaces
%    This should contain three filenames (each relative to subject dir):
%      (1) inner skull
%      (2) outer skull
%      (3) outer scalp
%    Optionally, for single-shell MEG-only gain matrix, a single filename
%     for the inner skull surface can be specified here.
%    Files must exist if bem_flag=1
%    {default = {'bem/inner_skull.tri'...
%               'bem/outer_skull.tri'...
%               'bem/outer_scalp.tri'}     }
%  'cen_sph': center of sphere (mm)
%     Required and must have three elements (x,y,z) if bem_flag=0
%     {default = []}
%  'radii': vector of three values specifying radii (mm) of spherical shell
%     layers: (1) inner skull (2) outer skull (3) outer scalp
%     Required if useEEG_flag=1 and bem_flag=0 (spherical shell model)
%     {default = []}
%  'trans': 4x4 matrix specifying mri2head transformation
%     {default = []}
%  'transfile': text file containing 4x4 mri2head transformation
%     {default = []}
%  'alignment_fif' fif file containing 4x4 mri2head transformation
%     note on reading from fif's:
%      this requires Uutela's fiff access toolbox (32 bit only)
%     and core dump will result if mri2head transform is missing from
%      specified file
%     {default = []}
%
% Notes on retmap:
%   retmap should contain wfiles or vertex numbers defining the set of possible
%     dipole locations for each stimulus location and visual area
%   for each iteration, a vertex within this neighborhood of possibility,
%     will be randomly chosen
%   if the neighborhood is too large, it will take a lot of memory and time to
%     calculate the forward matrix (done once for all possible locations)
%   if the neighborhood has only one vertex, then the same estimates will be
%     generated for each iteration
%
% Notes on ret_dips and nonret_dips:
%   ret_dips are additional dipoles without orientation constraints that are
%     also allowed to have a different source strength for each stimulus location
%   nonret_dips are additional dipoles without orientation constraints that are
%     constrained to have the same source strength (and orientation) for each
%     stimulus location
%
% Notes on transform files:
%   trans, trans_file, and alignment_fif are different, mutually exclusive
%     ways to specify the mri2head transformation matrix.
%   This transformation matrix is a 4x4 afine transformation that
%     maps locations in the MRI (from dip files) to locations in the
%     MEG/EEG head space defined by the digitized cardinal points
%     (and HPI coils).
%   This transformation matrix can also be found in the avg_data structure.
%     The avg_data structure contains a field called coor_trans.
%     If the subfield "mri2head" exists and is not empty, this will be used.
%   If more than one transformation marix is supplied, the order of dominance is:
%     trans, trans_file, alignment_fif, avg_data.coor_trans.mri2head
%     i.e., if trans is supplied directly, all other methods of
%     supplying the trans will be ignored.
%
% Notes on iterative fitting:
%  If nbrhd_niters>0, this program can be used to iteratively search for the
%   best fitting dipoles within a defined neighborhood for each stimulus
%   location for each visual area.  This program can be called multiple times,
%   seeking improvements of the estimates, by supplying the output
%   (retmap,best_retmap,min_error) as input on the next round.
%  If nbrhd_niters=0, this program will use the vertices defined in retmap,
%   without doing the random search.
%
% Created:  06/26/06 by Don Hagler
% Last Mod: 03/09/15 by Don Hagler
%

%% todo: change defaults so that grid_offset_flag = 1
%%       and r_step, th_step, r_offset_range, and th_offset_range
%%       are appropriate for grid offset

%% todo: remove grid_offset_smooth option?
%%       or use transform_coords function to get u,v coords for any vertices?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output
best_retmap = [];
min_error = [];

% parse input parameters
if ~mmil_check_nargs(nargin,2), return; end;
parms = parse_input(subjname,varargin);

% check input parameters, set dependent values
parms = check_input(parms,avg_data);

% create output directories
create_output_dirs(parms);

% prepare dipole information
parms = prep_dipinfo(parms);

% prepare data
[parms,avg_data] = prep_data(parms,avg_data);

% calculate scaling factors
parms = calc_scalefacts(parms,avg_data);

% prepare noise covariance matrix
parms = prep_noise(parms,avg_data);

% prepare vertex information
parms = prep_verts(parms);

% save parms to mat file for future reference
save_parms(parms);

fname_results = sprintf('%s/matfiles/%s_results.mat',...
  parms.rootoutdir,parms.prefix);
results = [];
if ~exist(fname_results,'file') || parms.forceflag
  % prepare forward solution
  if parms.verbose
    fprintf('%s: preparing forward solution...\n',mfilename);
  end;
  forward = prep_forward(parms,avg_data);

  if parms.rscale_niters
    [parms,results] = rc_RCSE_rand_rscale_search(parms,avg_data,forward);
  end;

  if ~isempty(parms.retfit_results)
    if parms.hybrid_flag
      [parms,results] = rc_RCSE_hybrid_offset_search(parms,avg_data,forward);    
    elseif parms.offset_niters>0  % multi-start offset search
      [parms,results] = rc_RCSE_mstart_offset_search(parms,avg_data,forward);    
    elseif parms.num_r_offsets>1 || parms.num_th_offsets>1
      % if r_offset or th_offset are vectors, do exhaustive search
      parms = rc_RCSE_full_offset_search(parms,avg_data,forward);
    else
      parms.retmap = ...
        rc_RCSE_retmap_from_retfit(parms,parms.r_offset,parms.th_offset,[],forward);
    end;
  else
    if parms.offset_niters>0
      error('with offset_niters>0, retfit_results must be supplied');
    end;
  end;
  if parms.rot_niters && ...
     parms.offset_niters>0 && parms.offset_mstarts>1 &&...
     parms.mstart_rotsearch_flag
    forward = parms.best_rot_forward;
  end;

  % random neighborhoods search for best fitting dipoles
  if parms.nbrhd_niters>0
    if parms.nbrhd_full_flag
      parms = rc_RCSE_full_nbrhood_search(parms,avg_data,forward);
    else
      parms = rc_RCSE_rand_nbrhood_search(parms,avg_data,forward);
    end;
  end;

  % search for receptive field sizes (unless we did it as part of mstart)
  if parms.rf_niters && ...
    ~(parms.offset_niters && parms.offset_mstarts>1 && parms.mstart_rfsearch_flag)
    [parms,results] = rc_RCSE_rand_rf_search(parms,avg_data,forward);
  end;

  % search for best overall rotation to account for registration errors
  if parms.rot_niters && ...
    ~(parms.offset_niters && parms.offset_mstarts>1 && parms.mstart_rotsearch_flag)
    [parms,results] = rc_RCSE_rand_rot_search(parms,avg_data,forward);
    forward = parms.best_rot_forward;
  end;

  if isempty(parms.best_retmap), parms.best_retmap = parms.retmap; end;
  if parms.reweight_flag
    [parms,tmp_results,retforward,inverse] =...
      rc_RCSE_reweight_lsq(parms,parms.best_retmap,avg_data,forward);
  else
    [parms,tmp_results,retforward,inverse]=...
      rc_RCSE_fit_waveforms(parms,parms.best_retmap,avg_data,forward);
  end;
  parms.best_retmap = retforward.retmap;
  if parms.verbose
    fprintf('%s: error = %f\n',mfilename,parms.total_error);
  end;
  if parms.mstart_average_flag && parms.offset_mstarts>1 && parms.offset_niters>0
    if parms.verbose
      fprintf('%s: using multi-start average waveforms...\n',mfilename);
    end;
  else
    results = tmp_results;
  end;

  % save results
  if parms.verbose
    fprintf('%s: saving results...\n',mfilename);
  end;
  parms = rc_RCSE_save_results(parms,avg_data,results,retforward,inverse,forward);

  best_retmap = parms.best_retmap;
  min_error = parms.min_error;
elseif parms.verbose
  fprintf('%s: results file %s already exists\n',mfilename,fname_results);
end;

if parms.plotflag
  if parms.verbose
    fprintf('%s: plotting...\n',mfilename);
  end;
  if isempty(results), load(fname_results); end;
  figure(1);
  rc_plot_sources('parms',parms,'results',results,...
    'area_names',parms.area_names,'area_colors',parms.area_colors,...
    'ylim',parms.plot_ylim,'units',parms.source_units);
  figure(2);
  rc_plot_fitvar_results(parms,results);
  drawnow;
end;

if parms.verbose
  fprintf('%s: finished.\n',mfilename);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = parse_input(subjname,options)
  parms = mmil_args2parms(options,{...
    'subjname',subjname,[],...
  ...
    'retmap',[],[],...
    'retfit_results',[],[],...
    'cond_info',[],[],...
    'subjdir',[],[],...
    'prefix','RCSE',[],...
    'rootoutdir',pwd,[],...
    'verbose',true,[false true],...
    'forceflag',0,[0:2],...
    'conditions',[],[],...
    'event_codes',[],[],...
    'use_areas',[],[],...
    'wfiledir','./weights',[],...
    'calc_dipinfo_flag',false,[false true],...
    'dip_matfile',[],[],...
    'lh_dip_file','bem/lh_white.dip',[],...
    'rh_dip_file','bem/rh_white.dip',[],...
    'lh_dip_info',[],[],...
    'rh_dip_info',[],[],...
    'lh_verts',[],[],...
    'rh_verts',[],[],...
    'inverse_type',1,[0:2],...
    'ncov_type',1,[0,1,2],...
    'raw_ncov_lambda',0.1,[0,1],...
    'SNR',10,[eps Inf],...
    'calc_scalefacts_flag',false,[false true],...
    'baseline_start',-Inf,[-Inf Inf],...
    'baseline_end',0,[-Inf Inf],...
    'bem_flag',true,[false true],...
    'openmeeg_flag',false,[false true],...
    'radii',[],[],...
    'cen_sph',[],[],...
    'conductivities',[0.3 0.012 0.3],[],...
    'EEG_gain_scalefact',1,[-Inf Inf],...
    'nlayers',3,[1:4],...
    'bem_surf_files',...
      {'bem/inner_skull.tri','bem/outer_skull.tri','bem/outer_scalp.tri'},[],...
    'grad_scalefact',10^13,[-Inf Inf],...
    'mag_scalefact',10^15,[-Inf Inf],...
    'EEG_scalefact',10^6,[-Inf Inf],...
    'badchans',[],[],...
    'badchanfile',[],[],...
    'usegrad_flag',true,[false true],...
    'usemag_flag',true,[false true],...
    'useEEG_flag',true,[false true],...
    'visual_area_weight',1,[0,1],...
    'ret_dips_weight',1,[0,1],...
    'nonret_dips_weight',1,[0,1],...
    'ret_dips_lh_dec_file',[],[],...
    'ret_dips_rh_dec_file',[],[],...
    'ret_dip_qfield_flag',true,[false true],...
    'nonret_dips_lh_dec_file',[],[],...
    'nonret_dips_rh_dec_file',[],[],...
    'norm_weights_flag',2,[0 1 2 3],...
    'surfname','white',[],...
    'w_thresh_patch',0,[0,1],...
    'near_nbr_weight',0,[0,1],...
    'indy_locs_flag',false,[false true],...
    'indy_smfact',0.999,[0,1],...
    'ecc_smfact',0.999,[0,1],...
    'upperlower_smfact',1,[0,1],...
    'hemi_smfact',1,[0,1],...
    'loose_flag',false,[false true],...
    'loose_tang_weight',0.1,[0,1],...
    'condF_thresh',0,[],...
    'write_err_flag',false,[false true],...
    'write_fit_flag',false,[false true],...
    'write_areafit_flag',false,[false true],...
    'write_fif_flag',false,[false true],...
    'forward_matfile',[],[],...
    'template_fif',[],[],...
    'trans',[],[],...
    'transfile',[],[],...
    'alignment_fif',[],[],...
    'reweight_flag',false,[false true],...
    'reweight_leverage_flag',true,[false true],...
    'reweight_leverage_max_flag',true,[false true],...
    'hybrid_flag',false,[false true],...
    'hybrid_pre_niters',0,[0,Inf],...
    'fmincon_flag',false,[false true],...
    'offset_niters',0,[0,Inf],...
    'offset_mstarts',0,[0,Inf],...
    'mstart_rotsearch_flag',false,[false true],...
    'mstart_rfsearch_flag',false,[false true],...
    'mstart_randstart_flag',true,[false true],...
    'mstart_average_flag',true,[false true],...
    'offset_niters_last',0,[0,Inf],...
    'offset_group_flag',1,[0 1 2 3 4],...
    'offset_group_patches_flag',false,[false true],...
    'offset_group_areas_flag',false,[false true],...
    'offset_const_areas',[],[],...
    'total_offset_niters',0,[0,Inf],...
    'total_offset_mstarts',0,[0,Inf],...
    'r_step',0.1,[1e-4,10],...
    'th_step',2,[1e-4,90],...
    'r_offset_range',[-3,3],[-100,100],...
    'th_offset_range',[-20,20],[-180,180],...
    'grid_offset_flag',false,[false true],...
    'grid_offset_smooth',0,[0,1000],...
    'total_rot_niters',0,[0,Inf],...
    'rot_niters',0,[0,Inf],...
    'rot_step',1,[0,180],...
    'rot_max',10,[0,180],...
    'best_rot',[0 0 0],[-180 180],...
    'total_rscale_niters',0,[0,Inf],...
    'rscale_niters',0,[0,Inf],...
    'rscale_step',0.01,[0,1],...
    'rscale_range',[0.8,1.2],[0.1,10],...
    'best_rscale',1,[0.1,10],...
    'nbrhd_niters',0,[0,Inf],...
    'total_nbrhd_niters',0,[0,Inf],...
    'nbrhd_full_flag',false,[false true],...
    'best_retmap',[],[],...
    'min_error',[],[],...
    'refEEG_coords',[],[],...
    'prior_prefix',[],[],...
    'prior_weight',0.5,[0,1],...
    'prior_avg_flag',0,[0:3],...
    'prior_mindiff_flag',0,[0:2],...
    'prior_zscore_flag',false,[false true],...
    'prior_sem_min',0.5,[eps,100],...
    'prior_indy_wform_flag',false,[false true],...
    'corr_time0',0,[-Inf,Inf],...
    'corr_time1',170,[-Inf,Inf],...
    'r_offset',0,[-100,100],...
    'th_offset',0,[-180,180],...
    'vf2ctx_flag',true,[false true],...
    'stim_type',2,[0:2],...
    'ecc_width',1,[0,100],...
    'theta_width',10,[0,360],...
    'w_thresh',0.01,[0,100],...
    'vfnorm_flag',true,[false true],...
    'single_vertex_flag',0,[0 1 2],...
    'r_max',12.5,[0,Inf],...
    'rf_sizes',[1,1,1],[0.01,10],...
    'rf_slopes',[0.1,0.1,0.1],[0,10],...
    'rf_r_inter',6,[0,Inf],...
    'restrict_hemi_flag',false,[false true],...
    'restrict_uplow_flag',false,[false true],...
    'retfit_data_flag',false,[false true],...
    'fit_range_flag',false,[false true],...
    'fit_time0',60,[-Inf,Inf],...
    'fit_time1',80,[-Inf,Inf],...
    'polarity_penalty',0,[-Inf,Inf],...
    'rf_niters',0,[0,Inf],...
    'total_rf_niters',0,[0,Inf],...
    'rf_min_sizes',[0.1,0.1,0.1],[0.01,10],...
    'rf_max_sizes',[3,3,3],[0.01,10],...
    'rf_min_slopes',[0,0,0],[0,10],...
    'rf_max_slopes',[0.2,0.2,0.2],[0,10],...
    'rf_size_step',0.1,[0.01,10],...
    'rf_slope_step',0.01,[0.001,1],...
    'best_rf_sizes',[],[0.01,10],...
    'best_rf_slopes',[],[0,10],...
  ... % scaling factors to get forward matrix into same scale as fif data
    'grad_unitsfact',1e-13,[0 1],... % fT/cm to T/m
    'mag_unitsfact',1e-15,[0 1],... % fT to T
    'EEG_unitsfact',1e-6,[0 1],... % uV to V
  ... % undocumented
    'plotflag',false,[false true],...
    'plot_ylim',[-30,30],[-Inf,Inf],...
    'area_names',{'v1','v2','v3'},[],...
    'area_colors',{'b','g','r'},[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'stimres',100,[50,1000],...
    'vf2ctx_flag',true,[false true],...
    'vf2ctx',[],[],...
    'save_avg_flag',true,[false true],...
    'evcode_flag',true,[false true],...
    'hybrid_cutoff',0.5,[0 1],...
    'reweight_factor',2,[],... % 4.685
    'reweight_maxiter',100,[],...
    'reweight_tol',1e-7,[],...
    'cond_weights',[],[],...
    'condF',[],[],...
    'nnbrs',10,[1,1000],...
  ...
    'chan_tags',{'badchans','badchanfile','usegrad_flag','usemag_flag',...
                 'useEEG_flag','channames','verbose',},[],...
    'gainmat_tags',{'fname_lh_surf','fname_rh_surf',...
                    'lh_dip_info','rh_dip_info','trans','bem_surf_files',...
                    'bem_matfile','lh_dec_dips','rh_dec_dips',...
                    'bem_flag','openmeeg_flag','useMEG_flag','useEEG_flag',...
                    'conductivities','cen_sph','radii',...
                    'rootoutdir','forceflag','verbose','batchsize'},[],...
    'vf2ctx_tags',{'w_thresh' 'vfnorm_flag' 'r_max'...
                   'rf_sizes' 'rf_slopes' 'rf_r_inter'...
                   'retfit_data_flag' 'area_names' 'hemilist' 'stimres'},[],...
  });
  parms.nhemi = length(parms.hemilist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(parms,avg_data)
  % initialize random number generation
  rand('twister',sum(100*clock));
  randn('state',sum(100*clock));

  % check subjdir
  if isempty(parms.subjdir)
    parms.subjdir = deblank(getenv('SUBJECTS_DIR'));
    if(isempty(parms.subjdir))
      error('cannot find SUBJECTS_DIR environment variable');
    end
  end;
  parms.subjpath = [parms.subjdir '/' parms.subjname];
  if ~exist(parms.subjpath,'dir')
    error('subject dir %s not found',parms.subjpath);
  end

  % check retfit_results vs retmap and areas to use
  if ~isempty(parms.retfit_results)
    parms.nareas = length(parms.area_names); % should be in retfit_results
    if isempty(parms.cond_info)
      error('if supplying retfit_results, must specify ''cond_info''');
    end;
    parms.orig_cond_info = parms.cond_info;
    if length(parms.cond_info) ~= length(avg_data.averages)
      error('length of cond_info must match number of conditions in avg_data.averages');
    end;
  elseif ~isempty(parms.retmap)
    parms.nareas = length(parms.retmap.areas);
  else
    error('must supply either retmap or retfit_results');
  end;

  % check condition info
  parms = rc_RCSE_check_conds(parms,avg_data);
  if isempty(parms.conditions)
    error('no valid conditions specified');
  end;

  % mri2head transformation matrix
  parms = set_trans(parms);

  % determine which channels (sensors) to use
  parms = set_chans(parms,avg_data);

  % check forward solution parameters
  parms = check_forward(parms);

  % check which areas will be used
  if isempty(parms.use_areas)
    parms.use_areas = 1:parms.nareas;
  end;
  parms.use_areas = find(ismember(1:parms.nareas,parms.use_areas));
  if isempty(parms.use_areas)
    error('use_areas is empty');
  end;
  % check that area indices are in bounds
  if ~isempty(parms.offset_const_areas)
    if length(intersect(parms.offset_const_areas,parms.use_areas)) ~= ...
       length(parms.offset_const_areas)
      error('invalid area indices in offset_const_areas');
    end;
  end;

  % check dipole info files
  if isempty(parms.dip_matfile) && ~parms.calc_dipinfo_flag
    if isempty(parms.lh_dip_info)
      if mmil_isrelative(parms.lh_dip_file)
        parms.lh_dip_file = [parms.subjpath '/' parms.lh_dip_file];
      end;
      if ~exist(parms.lh_dip_file)
        error('lh_dip_file %s not found',parms.lh_dip_file);
      end
    end;
    if isempty(parms.rh_dip_info)
      if mmil_isrelative(parms.rh_dip_file)
        parms.rh_dip_file = [parms.subjpath '/' parms.rh_dip_file];
      end;
      if ~exist(parms.rh_dip_file)
        error('rh_dip_file %s not found',parms.rh_dip_file);
      end
    end;
  end;

  % ret dip dec files
  if ~isempty(parms.ret_dips_lh_dec_file)
    if mmil_isrelative(parms.ret_dips_lh_dec_file)
      parms.ret_dips_lh_dec_file = [parms.subjpath '/' parms.ret_dips_lh_dec_file];
    end;
    if ~exist(parms.ret_dips_lh_dec_file)
      error('ret_dips_lh_dec_file %s not found',parms.ret_dips_lh_dec_file);
    end;
  end;
  if ~isempty(parms.ret_dips_rh_dec_file)
    if mmil_isrelative(parms.ret_dips_rh_dec_file)
      parms.ret_dips_rh_dec_file = [parms.subjpath '/' parms.ret_dips_rh_dec_file];
    end;
    if ~exist(parms.ret_dips_rh_dec_file)
      error('ret_dips_rh_dec_file %s not found',parms.ret_dips_rh_dec_file);
    end;
  end;

  % nonret dip dec files
  if ~isempty(parms.nonret_dips_lh_dec_file)
    if mmil_isrelative(parms.ret_dips_rh_dec_file)
      parms.nonret_dips_lh_dec_file = [parms.subjpath '/' parms.nonret_dips_lh_dec_file];
    end;
    if ~exist(parms.nonret_dips_lh_dec_file)
      error('nonret_dips_lh_dec_file %s not found',parms.nonret_dips_lh_dec_file);
    end;
  end;
  if ~isempty(parms.nonret_dips_rh_dec_file)
    if mmil_isrelative(parms.ret_dips_rh_dec_file)
      parms.nonret_dips_rh_dec_file = [parms.subjpath '/' parms.nonret_dips_rh_dec_file];
    end;
    if ~exist(parms.nonret_dips_rh_dec_file)
      error('nonret_dips_rh_dec_file %s not found',parms.nonret_dips_rh_dec_file);
    end;
  end;

  % check which type of constraint will be used
  if parms.loose_flag
    parms.indy_locs_flag=1;
  end;
  if parms.indy_locs_flag
    parms.inverse_type = 2;
  end;

  % options for constructing ret mapping
  tags = {'wfiledir','w_thresh_patch','norm_weights_flag',...
    'ret_dips_lh_dec_file','ret_dips_rh_dec_file',...
    'nonret_dips_lh_dec_file','nonret_dips_rh_dec_file',...
    'ret_dip_qfield_flag'};
  parms.ret_mapping_opts = mmil_parms2args(parms,tags);

  % intitialize fit error
  if isempty(parms.min_error)
    parms.min_error = 10^29;
  end;

  % for receptive field size optimization
  if parms.rf_niters>1 & isempty(parms.retfit_results)
    error('must supply either retfit_results to search for receptive field sizes');
  end;
  parms.init_rf_sizes = parms.rf_sizes;
  parms.init_rf_slopes = parms.rf_slopes;

  % for rotation optimization
  if parms.rot_max==0, parms.rot_max = 360; end;

  % prior prefix
  if ~isempty(parms.prior_prefix)
    if isempty(regexp(parms.prior_prefix,'^/')) % relative
      fname = sprintf('%s/matfiles/%s_results.mat',...
        parms.rootoutdir,parms.prior_prefix);
    else
      fname = sprintf('%s_results.mat',parms.prior_prefix);
    end;
    if ~exist(fname,'file')
      error('prior results file %s not found',fname);
    end;
    parms.fname_prior = fname;
  else
    parms.fname_prior = [];
  end;  

  % fit start and end samples
  sfreq = avg_data.sfreq;
  t_trigger = avg_data.averages(1).time(1)*1000;
  parms.fit_t0 = max(1,round((parms.fit_time0 - t_trigger)*sfreq/1000));
  parms.fit_t1 = min(round((parms.fit_time1 - t_trigger)*sfreq/1000),...
    length(avg_data.averages(1).time));
  nfit_timepoints = parms.fit_t1 - parms.fit_t0 + 1;
  if parms.nbrhd_niters>0 |...
     length(parms.r_offset)>0 | length(parms.th_offset)>0 |...
     parms.offset_niters>0
    if parms.fit_range_flag
      fprintf('%s: fitting times %d - %d msec post-stim\n',...
        mfilename,parms.fit_time0,parms.fit_time1);
      fprintf('%s: fitting samples %d - %d\n',...
        mfilename,parms.fit_t0,parms.fit_t1);
    elseif parms.polarity_penalty
      fprintf('%s: polarity penalty times %d - %d msec post-stim\n',...
        mfilename,parms.fit_time0,parms.fit_time1);
      fprintf('%s: polarity penalty samples %d - %d\n',...
        mfilename,parms.fit_t0,parms.fit_t1);
    end;
  end;

  if ~isempty(parms.prior_prefix) || parms.hybrid_flag
    % correlation start and end samples
    parms.corr_t0 = max(1,round((parms.corr_time0 - t_trigger)*sfreq/1000));
    parms.corr_t1 = min(round((parms.corr_time1 - t_trigger)*sfreq/1000),...
      length(avg_data.averages(1).time));
    fprintf('%s: reference correlation samples %d - %d\n',...
      mfilename,parms.corr_t0,parms.corr_t1);
  else
    parms.corr_t0 = [];
    parms.corr_t1 = [];
  end;

  % set units for plotting
  if parms.norm_weights_flag==3
    parms.source_units = 'nA m/mm^2';
    parms.plot_ylim = parms.plot_ylim;
  else
    parms.source_units = 'nA m';
  end;

  % check template for writing fif files
  if parms.write_fif_flag && ~exist(parms.template_fif,'file')
    error('template_fif %s not found',parms.template_fif);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_trans(parms)
  T_mri2head = [];
  if ~isempty(parms.trans)
    fprintf('%s: using transformation matrix supplied...\n',mfilename);
    T_mri2head = parms.trans;
  elseif ~isempty(parms.transfile)
    if ~exist(parms.transfile)
      error('transfile %s not found',parms.transfile);
    end;
    fprintf('%s: reading transformation matrix from %s...\n',...
      mfilename,parms.transfile);
    T_mri2head = ts_read_transfile(parms.transfile);
  elseif ~isempty(parms.alignment_fif)
    if ~exist(parms.alignment_fif)
      error('alignment_fif %s not found',parms.alignment_fif);
    end
    fprintf('%s: reading transformation matrix from %s...\n',...
      mfilename,parms.alignment_fif);
    T_mri2head=loadtrans(parms.alignment_fif,'MRI','HEAD');
  elseif isfield(avg_data.coor_trans,'mri2head')
    fprintf('%s: using transformation matrix in avg_data structure...\n',mfilename);
    if ~isempty(avg_data.coor_trans.mri2head)
      T_mri2head=avg_data.coor_trans.mri2head;
    end;
  end;
  if isempty(T_mri2head)
    error('no mri2head transformation was supplied');
  end;
  parms.trans = T_mri2head;
  fprintf('%s: mri2head transformation matrix:\n',mfilename);
  disp(parms.trans);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_chans(parms,avg_data)
  parms.nchans = length(avg_data.sensor_info);
  args = mmil_parms2args(parms,parms.chan_tags);
  [ind_chans,ind_grad,ind_mag,ind_EEG,ind_bad] = ts_set_chans(avg_data,args{:});
  parms.goodchans = ind_chans;
  parms.grad_chans = ind_grad;
  parms.mag_chans = ind_mag;
  parms.EEG_chans = ind_EEG;
  parms.badchans = ind_bad;
  if isempty(parms.goodchans)
    error('no valid channels available');
  end;
  if isempty(parms.EEG_chans)
    parms.useEEG_flag = 0;
  end;
  if isempty(parms.grad_chans) && isempty(parms.mag_chans)
    parms.useMEG_flag = 0;
  else
    parms.useMEG_flag = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_forward(parms)
  % check that nlayers is 1 or 3 if openmeeg ~= 1
  if ~parms.openmeeg_flag && ~ismember(parms.nlayers,[1,3])
    error('openmeeg_flag=1 required if number of layers not 1 or 3');    
  end;
  % check number of layers are consistent with MEG/EEG
  % EEG must have at least 3 layers
  if parms.nlayers<3 & parms.useEEG_flag
    error('must have 3 layers (BEM surfacse or spherical shells) for EEG');
  end;
  % bem vs. sphere
  if parms.bem_flag
    if iscell(parms.bem_surf_files)
      num_bem_surfs = length(parms.bem_surf_files);
    else
      parms.bem_surf_files = {parms.bem_surf_files};
      num_bem_surfs=1;
    end;
    if num_bem_surfs > parms.nlayers % i.e. nlayers = 1
      parms.bem_surf_files = {parms.bem_surf_files{1:parms.nlayers}};
      num_bem_surfs = length(parms.bem_surf_files);
    elseif num_bem_surfs < parms.nlayers
      error('number of bem surf files must match nlayers');
    end;
    % check that files exist and change from relative to absolute path
    for s=1:num_bem_surfs
      fname = parms.bem_surf_files{s};
      if mmil_isrelative(fname)
        fname = [parms.subjpath '/' fname];
      end;
      if ~exist(fname,'file')
        error('bem surf %d file %s not found',s,fname);
      end
      parms.bem_surf_files{s} = fname;
    end;
    if parms.verbose
      fprintf('%s: will use BEM forward with surface files:\n',mfilename);
      for s=1:num_bem_surfs
        fprintf('   %d: %s\n',s,parms.bem_surf_files{s});
      end;
    end;
  else
    if parms.useEEG_flag
      if length(parms.radii) > parms.nlayers
        parms.radii = parms.radii(1:parms.nlayers);
      elseif length(parms.radii) < parms.nlayers
        error('number of radii must match number of layers');
      end;
    end;
    if isempty(parms.cen_sph)
      error('cen_sph must be specified if bem_flag=0');
    end
    if(length(parms.cen_sph)~=3)
      error('cen_sph (center of sphere) incorrectly specified');
    end
    if parms.verbose
      fprintf('%s: spherical shell forward with center of sphere: [%0.1f,%0.1f,%0.1f]\n',...
        mfilename,parms.cen_sph);
      if parms.nlayers==1
        fprintf('%s: with radius: %0.1f',...
          mfilename,parms.radii);
      else
        fprintf('%s: with radii: [%0.1f,%0.1f,%0.1f]\n',...
          mfilename,parms.radii);
      end;
    end;
  end
  if length(parms.conductivities) > parms.nlayers
    parms.conductivities = parms.conductivities(1:parms.nlayers);
  elseif length(parms.conductivities) < parms.nlayers
    error('number of conductivities must match number of layers');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_output_dirs(parms)
  mmil_mkdir([parms.rootoutdir '/matfiles']);
  if parms.write_fif_flag
    mmil_mkdir([parms.rootoutdir '/fifs']);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_parms(parms)
  fname = sprintf('%s/matfiles/%s_parms.mat',parms.rootoutdir,parms.prefix);
  if parms.verbose
    fprintf('%s: saving parameters to mat file %s...\n',mfilename,fname);
  end;
  save(fname,'parms');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_dipinfo(parms)
  fname = sprintf('%s/matfiles/%s_dipinfo.mat',parms.rootoutdir,parms.prefix);
  if ~exist(fname,'file') || parms.forceflag
    % load or calculate dipole info
    if ~isempty(parms.dip_matfile)
      if ~exist(parms.dip_matfile,'file')
        error('dip matfile %s not found',parms.dip_matfile);
      end;
      if parms.verbose
        fprintf('%s: loading dipole info from mat file %s...\n',...
          mfilename,parms.dip_matfile);
      end;
      load(parms.dip_matfile);
    elseif parms.calc_dipinfo_flag
      if parms.verbose
        fprintf('%s: calculating dipole information...\n',mfilename);
      end;
      [lh_dip_info,rh_dip_info] = ts_dip_info(parms.subjname,parms.subjdir);  
    else
      if isempty(parms.lh_dip_info)
        if parms.verbose
          fprintf('%s: loading left hemisphere dipole information from %s...\n',...
            mfilename,parms.lh_dip_file);
        end;
        lh_dip_info = ts_read_dip_file(parms.lh_dip_file);
      end;
      if isempty(parms.rh_dip_info)
        if parms.verbose
          fprintf('%s: loading right hemisphere dipole information from %s...\n',...
            mfilename,parms.rh_dip_file);
        end;
        rh_dip_info = ts_read_dip_file(parms.rh_dip_file);
      end;
    end;

    % ret dip dec files
    if ~isempty(parms.ret_dips_lh_dec_file)
      ret_dips_lh_dec_dips = ts_read_dec_file(parms.ret_dips_lh_dec_file);
    else
      ret_dips_lh_dec_dips = [];
    end;
    if ~isempty(parms.ret_dips_rh_dec_file)
      ret_dips_rh_dec_dips = ts_read_dec_file(parms.ret_dips_rh_dec_file);
    else
      ret_dips_rh_dec_dips = [];
    end;

    % nonret dip dec files
    if ~isempty(parms.nonret_dips_lh_dec_file)
      nonret_dips_lh_dec_dips = ts_read_dec_file(parms.nonret_dips_lh_dec_file);
    else
      nonret_dips_lh_dec_dips = [];
    end;
    if ~isempty(parms.nonret_dips_rh_dec_file)
      nonret_dips_rh_dec_dips = ts_read_dec_file(parms.nonret_dips_rh_dec_file);
    else
      nonret_dips_rh_dec_dips = [];
    end;

    % save dipole and decimation info
    if parms.verbose
      fprintf('%s: saving dipinfo to mat file %s...\n',mfilename,fname);
    end;
    save(fname,'lh_dip_info','rh_dip_info',...
               'ret_dips_lh_dec_dips','ret_dips_rh_dec_dips',...
               'nonret_dips_lh_dec_dips','nonret_dips_rh_dec_dips');
  else
    if parms.verbose
      fprintf('%s: loading dipinfo from mat file %s...\n',mfilename,fname);
    end;
    load(fname);
  end;
  parms.lh_dip_info = lh_dip_info;
  parms.rh_dip_info = rh_dip_info;
  parms.ret_dips_lh_dec_dips = ret_dips_lh_dec_dips;
  parms.ret_dips_rh_dec_dips = ret_dips_rh_dec_dips;
  parms.nonret_dips_lh_dec_dips = nonret_dips_lh_dec_dips;
  parms.nonret_dips_rh_dec_dips = nonret_dips_rh_dec_dips;
  parms.num_dips_lh = length(parms.lh_dip_info);
  parms.num_dips_rh = length(parms.rh_dip_info);
  if(parms.num_dips_lh==0 & parms.num_dips_rh==0)
    error('no dipole info found');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,avg_data] = prep_data(parms,avg_data)
  fname = sprintf('%s/matfiles/%s_avg_data.mat',...
    parms.rootoutdir,parms.prefix);
  if ~parms.save_avg_flag || ~exist(fname,'file') || parms.forceflag
    % replace/add reference coords to loc for each electrode
    if parms.useEEG_flag & ~isempty(parms.refEEG_coords)
      if length(parms.refEEG_coords)~=3
        error('regEEG_coords must contain 3 values (x,y,z)');
      end;
      for i=1:parms.EEG_chans
        avg_data.sensor_info(i).loc(1:3,1) = parms.refEEG_coords;
      end;
    end;
    if parms.save_avg_flag
      if parms.verbose
        fprintf('%s: saving avg_data to %s...\n',mfilename,fname);
      end;
      save(fname,'avg_data');
    end;
  else
    if parms.verbose
      fprintf('%s: loading avg_data from %s...\n',mfilename,fname);
    end;
    load(fname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_scalefacts(parms,avg_data)
  % calculate scaling factors for each channel (if necessary)
  if parms.calc_scalefacts_flag || parms.ncov_type==1
    parms = calc_channel_scalefacts(parms,avg_data);
  end;
  % calculate scaling factors for each channel type
  %   otherwise use default or user supplied scaling factors
  if parms.calc_scalefacts_flag
    parms = calc_chantype_scalefacts(parms);
  end;
  % identify additional bad chans based on flat baseline
  parms = detect_flat_chans(parms,avg_data);
  if parms.verbose
    fprintf('%s: scaling factors: grad=%0.2g, mag=%0.2g, EEG=%0.2g\n',...
      mfilename,parms.grad_scalefact,parms.mag_scalefact,parms.EEG_scalefact);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_channel_scalefacts(parms,avg_data)
  if parms.verbose
    fprintf('%s: calculating scaling factors from data...\n',...
      mfilename);
  end;
  parms.scalefacts = zeros(parms.nchans,1);
  total_num_trials = 0;
  nconds = length(parms.conditions);
  for c=1:nconds
    cond = parms.conditions(c);
    total_num_trials = total_num_trials + avg_data.averages(cond).num_trials;
    time = avg_data.averages(cond).time;
    [tmp,baseline_start_samp] = min(abs(time - parms.baseline_start/1000));
    [tmp,baseline_end_samp] = min(abs(time - parms.baseline_end/1000));
    data = avg_data.averages(cond).data;
    baseline = data(:,baseline_start_samp:baseline_end_samp);
    num_trials = avg_data.averages(cond).num_trials;
    if isempty(num_trials) | num_trials<2, num_trials = 2; end;
    parms.scalefacts = parms.scalefacts + std(baseline,0,2)*sqrt(num_trials-1);
  end;
  if isempty(total_num_trials) | total_num_trials<nconds+1
    if parms.ncov_type==2
      fprintf('%s: WARNING: invalid number of trials will result in invalid noise normalization\n',...
        mfilename);
    end;
    total_num_trials = nconds+1; % prevent division by 0
  end;
  % pooled standard deviation across conditions
  parms.scalefacts = parms.scalefacts/sqrt(total_num_trials-nconds);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_chantype_scalefacts(parms)
  if isempty(parms.grad_chans)
    parms.grad_scalefact = 0;
  else
    parms.grad_scalefact = 1.0/mean(parms.scalefacts(parms.grad_chans));
  end;
  if isempty(parms.mag_chans)
    parms.mag_scalefact = 0;
  else
    parms.mag_scalefact = 1.0/mean(parms.scalefacts(parms.mag_chans));
  end;
  if isempty(parms.EEG_chans)
    parms.EEG_scalefact = 0;
  else
    parms.EEG_scalefact = 1.0/mean(parms.scalefacts(parms.EEG_chans));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = detect_flat_chans(parms,avg_data)
  % multiply channel scale facts by channel type scale facts
  scale_vector = calc_scale_vector(parms);
  scalefacts = scale_vector .* parms.scalefacts;
  scalefacts = scalefacts(parms.goodchans);
  i_bad = find(scalefacts<eps);
  if ~isempty(i_bad)
    badchans = parms.goodchans(i_bad);
    sensor_labels = {avg_data.sensor_info.label};
    badchan_labels = sensor_labels(badchans);
    fprintf('%s: WARNING: autodetected these channels as flat: %s\n',...
      mfilename,sprintf('%s ',badchan_labels{:}));
    parms.badchans = [parms.badchans,badchans];
    parms.goodchans = setdiff(parms.goodchans,badchans);
    if isempty(parms.goodchans)
      error('all flat channels');
    end;
    parms.grad_chans = setdiff(parms.grad_chans,parms.badchans);
    parms.mag_chans = setdiff(parms.mag_chans,parms.badchans);
    parms.EEG_chans = setdiff(parms.EEG_chans,parms.badchans);
    % recalculate scaling factors for each channel type
    if parms.calc_scalefacts_flag && ~parms.prewhiten_flag
      parms = calc_chantype_scalefacts(parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_noise(parms,avg_data)
  % use calculated noise covariance matrix or get from avg_data.noise.covar
  if parms.ncov_type==0
    if parms.verbose
      fprintf('%s: using scaled identity matrix for noise covariance\n',...
        mfilename);
    end;
    C = create_identity_noisecovar(parms);
    navg = 1;
  elseif parms.ncov_type==1
    if parms.verbose
      fprintf('%s: using noise covariance calculated from average baseline\n',...
        mfilename);
    end;
    C = diag(parms.scalefacts.^2);
    navg = length(parms.conditions);
  elseif isempty(avg_data.noise.covar)
    if parms.verbose
      fprintf('%s: using scaled identity matrix for noise covariance\n',...
        mfilename);
    end;
    C = create_identity_noisecovar(parms);
    navg = 1;
  else
    if parms.verbose
      fprintf('%s: using noise covariance matrix from avg_data.noise.covar\n',...
        mfilename);
    end;
    C = avg_data.noise.covar;
    L = parms.raw_ncov_lambda; % regularization
    C = (1-L)*C + L*diag(diag(C));
    navg = avg_data.noise.num_trials;
  end;
  parms.noisecovar = C;
  parms.noisecovar_navg = navg;

  % create matrix for scaling noise covariance matrix
  scale_matrix = calc_noise_scale_matrix(parms);

  % scale noise covariance matrix
  parms.noisecovar = parms.noisecovar.*scale_matrix;

  % exclude bad channels from noise covariance matrix
  parms.noisecovar = ...
    parms.noisecovar(parms.goodchans,parms.goodchans);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = create_identity_noisecovar(parms)
  scalefacts = zeros(parms.nchans,1);
  scalefacts(parms.grad_chans) = 1.0/(parms.grad_scalefact^2);
  scalefacts(parms.mag_chans) = 1.0/(parms.mag_scalefact^2);
  scalefacts(parms.EEG_chans) = 1.0/(parms.EEG_scalefact^2);
  C = diag(scalefacts);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scale_vector = calc_scale_vector(parms)
  scale_vector = zeros(parms.nchans,1);
  scale_vector(parms.grad_chans) = parms.grad_scalefact;
  scale_vector(parms.mag_chans) = parms.mag_scalefact;
  scale_vector(parms.EEG_chans) = parms.EEG_scalefact;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scale_matrix = calc_noise_scale_matrix(parms)
  scale_vector = calc_scale_vector(parms);
  scale_matrix = scale_vector*scale_vector';  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scale_matrix = calc_data_scale_matrix(parms,avg_data)
  [num_sensors,num_tpoints] = size(avg_data.averages(1).data);
  scale_vector = calc_scale_vector(parms);
  scale_matrix = scale_vector*ones(1,num_tpoints);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scale_matrix = calc_forward_scale_matrix(parms,G_xyz)
  [num_sensors,num_sources] = size(G_xyz);
  % scale to get to units of data, then scale with scalefacts
  %   also scale EEG part with EEG_gain_scalefact
  scale_vector = zeros(parms.nchans,1);
  scale_vector(parms.grad_chans) = parms.grad_unitsfact*parms.grad_scalefact;
  scale_vector(parms.mag_chans) = parms.mag_unitsfact*parms.mag_scalefact;
  scale_vector(parms.EEG_chans) = parms.EEG_unitsfact*parms.EEG_scalefact*...
                                            parms.EEG_gain_scalefact;
  scale_matrix = scale_vector*ones(1,num_sources);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_verts(parms)
  % calculate vertex area
  parms = calc_vertex_areas(parms);

  % pre-compute mapping between visual field and cortex
  parms = calc_vf2ctx(parms);

  % pre-compute neighbors for each vertex in roi
  parms = calc_nbrhds(parms);

  % define set of vertices to be used in forward solution
  parms = set_verts(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_vertex_areas(parms)
  if parms.norm_weights_flag==3
    fname = sprintf('%s/matfiles/%s_vertex_area.mat',...
      parms.rootoutdir,parms.prefix);
    if ~exist(fname,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: calculating vertex areas...\n',mfilename);
        tic;
      end;
      vertex_areas = [];
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_surf = [parms.subjpath '/surf/' hemi '.' parms.surfname];
        surf = fs_read_surf(fname_surf);
        msurf = preprocessQ(surf);
        vertex_areas.(hemi) = areaPerVertex(msurf);
      end;
      if parms.verbose, toc; end;
      save(fname,'vertex_areas');
    else
      load(fname);
    end;
    parms.lh_vertex_areas = vertex_areas.lh;
    parms.rh_vertex_areas = vertex_areas.rh;
  else
    parms.lh_vertex_areas = [];
    parms.rh_vertex_areas = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_vf2ctx(parms)
  if parms.vf2ctx_flag && ~isempty(parms.retfit_results)
    fname = sprintf('%s/matfiles/%s_vf2ctx.mat',...
      parms.rootoutdir,parms.prefix);
    if ~exist(fname,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: calculating mapping between visual field and cortex...\n',...
          mfilename);
        tic
      end;
      args = mmil_parms2args(parms,parms.vf2ctx_tags);
      [vf2ctx,retfit_results] = rc_calc_vf2ctx(parms.retfit_results,args{:});
      if parms.verbose, toc; end;
      save(fname,'vf2ctx','retfit_results');
    else
      load(fname);
    end;
    parms.vf2ctx = vf2ctx;
    parms.retfit_results = retfit_results;
  end;    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_nbrhds(parms)
  nbrhds = [];
  if parms.grid_offset_flag && ~isempty(parms.retfit_results)
    fname = sprintf('%s/matfiles/%s_nbrhds.mat',...
      parms.rootoutdir,parms.prefix);
    if ~exist(fname,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: calculating neighborhoods...\n',mfilename);
        tic;
      end;
      clear surfs;
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_surf = sprintf('%s/%s/surf/%s.%s',...
          parms.subjdir,parms.subjname,hemi,parms.surfname);
        surf = fs_read_surf(fname_surf);
        surfQ = preprocessQ(surf);
        fit_data = parms.retfit_results.(sprintf('%s_fit_data',hemi));
        roi = find(fit_data.u | fit_data.v);
        nbrhds(h).n = cell(1,length(roi));
        nbrhds(h).verts = [];
        for k=1:length(roi)
          n0 = roi(k);
          v = neighborsQ(surfQ,n0,parms.nnbrs);
          nbrhds(h).n{k} = v;
          nbrhds(h).verts = union(nbrhds(h).verts,v);
        end;
        surfs(h) = surfQ;
      end;
      if parms.verbose, toc; end;
      save(fname,'nbrhds','surfs');
    else
      load(fname);
    end;
  end;
  parms.nbrhds = nbrhds;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_verts(parms)
  if ~isempty(parms.retfit_results)
    parms.lh_verts = find(parms.retfit_results.lh_fit_data.pol_r~=0 |...
                          parms.retfit_results.lh_fit_data.pol_i~=0);
    parms.rh_verts = find(parms.retfit_results.rh_fit_data.pol_r~=0 |...
                          parms.retfit_results.rh_fit_data.pol_i~=0);
  elseif isempty(parms.lh_verts) | isempty(parms.rh_verts)
    % construct ret mapping structure for all areas
    use_areas = [1:length(parms.retmap.areas)];
    % use unique_location_conds only (e.g. multiple contrast levels)
    tmp_event_codes = []; tmp_conditions = [];
    if ~isempty(parms.event_codes)
      tmp_event_codes = parms.event_codes(parms.unique_location_conds);
    end;
    if ~isempty(parms.conditions)
      tmp_conditions = parms.conditions(parms.unique_location_conds);
    end;
    retforward.retmap = rc_construct_ret_mapping(retmap,...
      'use_areas',parms.use_areas,'conditions',tmp_conditions,...
      'event_codes',tmp_event_codes,parms.ret_mapping_opts{:});
    if isempty(retmap)
      error('failed to create retinotopy mapping matrix');
    end
    if ~isempty(parms.lh_verts) % empty if no verts in this hemi
      parms.lh_verts = retmap.uniq_verts_lh;
    end;  
    if ~isempty(parms.rh_verts)
      parms.rh_verts = retmap.uniq_verts_rh;
    end;  
  end;
  % include all verts in nbrhds that are in ROI to forward solution
  if ~isempty(parms.nbrhds)
    lh_roi = find(parms.retfit_results.lh_fit_data.u |...
                  parms.retfit_results.lh_fit_data.v);
    rh_roi = find(parms.retfit_results.rh_fit_data.u |...
                  parms.retfit_results.rh_fit_data.v);
    parms.lh_verts = intersect(union(parms.lh_verts,parms.nbrhds(1).verts),lh_roi);
    parms.rh_verts = intersect(union(parms.rh_verts,parms.nbrhds(2).verts),rh_roi);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function forward = prep_forward(parms,avg_data)
  % calculate forward solution for all possible dipoles
  forward = [];
  forward.lh_verts = parms.lh_verts;
  forward.rh_verts = parms.rh_verts;
  if ~isempty(parms.forward_matfile)
    if mmil_isrelative(parms.forward_matfile)
      fname_forward = sprintf('%s/matfiles/%s',...
        parms.rootoutdir,parms.forward_matfile);
    else
      fname_forward = parms.forward_matfile;
    end;
  else
    fname_forward = sprintf('%s/matfiles/%s_forward.mat',...
      parms.rootoutdir,parms.prefix);
  end;
  fname_forward_prep=sprintf('%s/matfiles/%s_forward_prep.mat',...
    parms.rootoutdir,parms.prefix);
  if ~exist(fname_forward_prep,'file') || parms.forceflag
    if ~exist(fname_forward,'file') || parms.forceflag>1
      % expand lh_verts and rh_verts to include buffer zone
      %
      %   NOTE: this is not necessary if retfit used roi_dilate_niters
      %           to create buffer zone
      %         also, just because vertices are included in forward
      %           does not mean they will be used with grid_offset
      %         that requires that they have u,v values in retfit,
      %           which requires that they are in the retfit ROI
      %
      %         a potential solution is to transform the coordinates
      %           of additional vertices using roi_rotation, etc. from retfit
      %
      if parms.grid_offset_flag && parms.grid_offset_smooth>0
        for h=1:parms.nhemi
          hemi = parms.hemilist{h};
          verts_tag = [hemi '_verts'];
          surf = fs_load_subj(parms.subjname,hemi,parms.surfname,0,parms.subjdir);
          vals = zeros(surf.nverts,1);
          vals(parms.(verts_tag)) = 1;
          vals = fs_smooth_sparse(surf,vals,parms.grid_offset_smooth);
          forward.(verts_tag) = find(vals);
        end;
      end;
      forward.lh_dip_info = parms.lh_dip_info;
      forward.rh_dip_info = parms.rh_dip_info;
      forward.lh_dec_dips = zeros(parms.num_dips_lh,1);
      forward.lh_dec_dips(forward.lh_verts) = 1;
      forward.num_dec_dips_lh=length(find(forward.lh_dec_dips));
      forward.rh_dec_dips = zeros(parms.num_dips_rh,1);
      forward.rh_dec_dips(forward.rh_verts) = 1;
      forward.num_dec_dips_rh=length(find(forward.rh_dec_dips));
      parms.lh_dec_dips = forward.lh_dec_dips;
      parms.rh_dec_dips = forward.rh_dec_dips;    
      if parms.verbose
        fprintf('%s: number of dipoles in forward solution: lh=%d rh=%d\n',...
          mfilename,forward.num_dec_dips_lh,forward.num_dec_dips_rh);
      end;

      if parms.verbose
        fprintf('%s: calculating forward solution gain matrix...\n',mfilename);
      end;
      if parms.bem_flag
        parms.bem_matfile = sprintf('%s/matfiles/%s_bem.mat',...
          parms.rootoutdir,parms.prefix);
      end;
      args = mmil_parms2args(parms,parms.gainmat_tags);
      forward.G_xyz = ts_calc_gainmat(avg_data,args{:});
      if isempty(forward.G_xyz)
        error('failed to calculate gain matrix');
      end;
      if parms.verbose
        fprintf('%s: saving forward solution to matfile %s...\n',...
          mfilename,fname_forward);
      end;
      save(fname_forward,'forward');
    else
      if parms.verbose
        fprintf('%s: loading pre-calculated forward solution from %s...\n',...
          mfilename,fname_forward);
      end;
      load(fname_forward);
    end

    % scale gain matrix
    if parms.verbose
      fprintf('%s: scaling gain matrix...\n',mfilename);
    end;
    scale_matrix = calc_forward_scale_matrix(parms,forward.G_xyz);
    forward.G_xyz = forward.G_xyz.*scale_matrix;

    % scale gain matrix by vertex area
    if parms.norm_weights_flag==3
      forward.lh_vertex_areas = parms.lh_vertex_areas(forward.lh_verts);
      forward.rh_vertex_areas = parms.rh_vertex_areas(forward.rh_verts);
      scalefacts = [forward.lh_vertex_areas;forward.rh_vertex_areas];
      scalefacts = reshape([scalefacts,scalefacts,scalefacts]',...
        [size(forward.G_xyz,2),1])';
      scale_matrix = ones(size(forward.G_xyz,1),1)*scalefacts;
      forward.G_xyz = forward.G_xyz.*scale_matrix;
    end;

    % convert xyz to normal and two orthogonal tangential components
    if parms.verbose
      fprintf('%s: converting xyz components to normal and tangentials\n',...
        mfilename);
    end;
    [forward.G_norm,forward.G_tang1,forward.G_tang2]=ts_gain_xyz2norm(...
      forward.G_xyz,...
      forward.lh_dip_info,forward.rh_dip_info,...
      forward.lh_dec_dips,forward.rh_dec_dips,parms.trans);

    % resize to remove unwanted and bad chans
    if parms.verbose
      fprintf('%s: resizing gain matrix to remove bad channels...\n',mfilename);
    end;
    forward.G_xyz = forward.G_xyz(parms.goodchans,:);
    forward.G_norm = forward.G_norm(parms.goodchans,:);
    forward.G_tang1 = forward.G_tang1(parms.goodchans,:);
    forward.G_tang2 = forward.G_tang2(parms.goodchans,:);

    % create scale matrix for data
    forward.scale_matrix = calc_data_scale_matrix(parms,avg_data);
    forward.inv_scale_matrix = ones(size(forward.scale_matrix));
    forward.inv_scale_matrix(find(forward.scale_matrix)) = ...
      1./forward.scale_matrix(find(forward.scale_matrix));

    save(fname_forward_prep,'forward');
  else
    if parms.verbose
      fprintf('%s: loading pre-prepared forward solution from %s...\n',...
        mfilename,fname_forward_prep);
    end;
    load(fname_forward_prep);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

