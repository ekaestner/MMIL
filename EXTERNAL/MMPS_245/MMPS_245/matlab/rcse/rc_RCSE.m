function prefix = rc_RCSE(subj,varargin)
%function prefix = rc_RCSE(subj,[options])
%
% Required Parameters:
%   subj: freesurfer subject (directory name)
%
% Optional Parameters:
%  'rootoutdir': root output directory
%    (subdirectories will be created and/or assumed to exist)
%    {default = pwd} (current working directory)
%  'prefix': prefix of all output files
%    (relative to rootoutdir/matfiles)
%     {default = 'RCSE'}
%  'infix': additional string attached to output prefix
%     e.g. to indicate values of key parameters
%     {default = []}
%  'infix_flag': [0|1] attach automatically generated string
%     indicating values of certain parameters
%     ignored if 'infix' is not empty
%     {default = 0}
%  'datamatfile': name of mat file containing averages sensor data
%    If not absolute path, assumed to be relative to rootoutdir/matfiles
%    {default = 'proc_avg_data_subnull.mat'}
%  'subjdir': freesurfer root subject dir
%    {default = SUBJECTS_DIR} (environment variable)
%  'plotflag': [0|1|2|3] whether to plot source waveforms and residual error
%     0=no plots
%     1=save plots at end (do not display)
%     2=display plots during offset and nbrhd fitting (do not save)
%     3=do not run RCSE, only save plots from previous run (do not display)
%     {default = 0}
%  'transfile': name of .trans file with mri to head registration matrix
%     {default = 'mri2head.trans'}
%  'badchanfile': name of text file containing bad channel labels
%     {default = []}
%  'forward_matfile': mat file containing gain matrix
%     If does not exist, will calculate
%     {default = []}
%  'graymid_flag': [0|1] use graymid dip files (must exist) instead of white
%     {default = 0}
%  'calc_dipinfo_flag': [0|1] whether to calculate dipole info from surfaces
%     otherwise load from .dip files
%    {default = 0}
%  'inverse_type': [0|1|2] type of inverse calculations
%    0: unregularized pseudo-inverse (fast, very little memory)
%    1: regularized psuedo-inverse with identity matrix for noise covariance
%       no source covariance matrix (fast, very little memory)
%    2: regularized psuedo-inverse with noise covariance matrix from data
%       depending on ncov_type, uses source covariance matrix
%    {default = 1}
%  'ncov_type': specify what type of noise covariance matrix to use
%     if 0, use identity matrix
%       (assume uniform white noise, independently scaled for each sensor type)
%     if 1, calculate and use noise covariance matrix from average prestim
%     if 2, use noise covariance matrix calculated from single trials
%       (stored in avg_data.noise.covar)
%    {default = 1}
%  'SNR': estimated signal-to-noise-ratio (for regularization parameter)
%      {default = 10}
%  'verbose': [0|1] display status messages
%     {default = 1}
%  'forceflag': [0|1|2] whether to overwrite files
%     0: do nothing (except maybe plots) if output exists
%     1: overwrite results and other output
%     2: overwrite all output, including forward solution
%     {default = 0}
%
% Optional Parameters about retfit:
%  'retfit_dir': directory containing retfit results
%    {default = 'retfit'}
%  'retfit_stem': file stem for refit results
%    {default = 'retfit'}
%  'vf2ctx_flag': [0|1] precompute mapping from visual field to cortex
%    {default = 1}
%  'stim_type': type of stimulus to model
%    0: point: generation of retmap is very fast, but less accurate
%    1: circle: slower, but more accurately takes into account receptive field size
%    2: wedge: slower, but more accurately takes into account receptive field size
%    {default = 2}
%  'norm_weights_flag': [0|1|2|3] whether and how to normalize vertex weights
%     0: do not normalize
%     1: normalize so sum of weights for each stimulus location equals 1
%     2: normalize so average sum of weights (within a visual area) equals 1
%     3: multiply forward matrix by vertex areas (changes source units to nA m/mm^2)
%     {default = 2}
%  'use_areas': vector of visual area indices to include in forward model
%    If empty, use all areas specified in retfit results
%    {default = []}
%  'r_offset': vector of r_offsets
%    {default = 0}
%  'th_offset': vector of th_offsets
%    {default = 0}
%  'ecc_width': eccentricity width of stimuli (deg. vis. ang.)
%    will be ignored if fname_conds has ecc_width column
%    {default = 1}
%  'theta_width': polar angle width of stimuli (degrees)
%    will be ignored if fname_conds has theta_width column or if stim_type=1
%    {default = 10}
%  'w_thresh': threshold applied to weights
%    {default = 0.01}
%  'vfnorm_flag': [0|1] for applying w_thresh for vf2ctx
%     0: normalize to global max
%     1: normalize to max for each cortical location
%     {default = 1}
%  'w_thresh_patch' - threshold applied to cortical patch weights
%       relative to max val for each patch
%     this reduces the number of dipoles that included in the forward matrix
%       cutting off tails of smooth, gaussian-like cluster
%     {default: 0}
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
% Optional Parameters for nonlinear fitting:
%  'fit_range_flag': [0|1] for nonlinear fits, minimize error between
%     fit_time0 and fit_time1 (otherwise use entire time range)
%     {default = 0}
%  'fit_time0': start of time range (msec) to fit data
%     {default = 70}
%  'fit_time1': end of time range (msec) to fit data
%     {default = 100}
%  'polarity_penalty': cost multiplier to penalize waveforms with positive
%     or negative polarity
%     If positive, source waveforms with positive mean amplitude within
%       fit time window will be penalized according to this formula:
%         total_error = total_error*(1 + mean_amp*polarity_penalty)
%     If negative, source waveforms with negative mean amplitude within
%       fit time window will be penalized
%     If zero, no polarity penalty
%     {default = 0}
%  'best_retmap_prefix': prefix of output files from a previous run of RCSE
%     Example usage: find best fitting dipoles in neighborhood using
%       indy_locs_flag = 0, then rerun using this option with indy_locs_flag = 1
%     Note: if supplied, no fitting will be allowed: prior_prefix, nbrhd_niters,
%       offset_niters, r_offset, and th_offset will be ignored
%     {default = []}
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
%  'corr_time0': start of time range (msec) to calculate correlation with reference
%     Note: only applicable if prior_prefix not empty
%     {default = 0}
%  'corr_time1': end of time range (msec) to calculate correlation with reference
%     Note: only applicable if prior_prefix not empty
%     {default = 170}
%  'err_prefix': prefix of output files from a previous run of RCSE
%     If supplied, will fit residual error from previous run
%     {default = []}
%  'offset_niters': number of iterations for random offset search to fit data
%     This tests neighboring vertices, defined by retmap from retfit or w files,
%     to see if they give a better fit to the data
%     If 0, do exhaustive search of offsets
%       specified by all combinations of r_offset and th_offset vectors
%     {default = 0}
%  'max_offset_niters': maximum number of iterations for random offset search
%     {default = Inf}
%  'offset_mstarts': number of multi-start iterations
%       for random offset search (ignored if offset_niters==0)
%     If 0, use 0 offsets for starting point
%     If >=1, randomly choose offsets for starting point
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
%    If 0, no bounds
%    {default = 0}
%  'th_offset_range': vector of min and max th_offsets for rand offset search
%    If 0, no bounds
%    {default = 0}
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
%  'rscale_niters': number of iterations for eccentricity (r) scaling
%     This adjusts the mapping between MEG and MRI stimuli (scales r_max)
%     {default = 0}
%  'rscale_step': standard deviation of rscale step size
%     {default = 0.01}
%  'rscale_range': vector of min and max rscale
%     {default = [0.8,1.2]}
%  'nbrhd_niters': number of iterations for random dipole-search to fit data
%     This tests neighboring vertices, defined by retmap from retfit or w files,
%     to see if they give a better fit to the data
%     {default = 0}
%  'nbrhd_full_flag': [0|1] for nbrhd search, test one stimulus location
%     at a time and do exhaustive search of all neighboring vertices
%     {default = 0}
%  'rf_niters': number of iterations for random search for best-fitting
%     receptive field sizes for each visual area
%     {default = 0}
%  'rf_min_sizes': vector of minimum receptive field sizes
%     one value for each visual area
%     {default = [0.1,0.1,0.1]}
%  'rf_max_sizes': vector of minimum receptive field sizes
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
%  'reweight_flag': [0|1] whether to use iteratively re-weighted least squares
%     {default = 0}
%  'reweight_factor': tunable reweighting factor or IRLS
%     {default = 2}
%  'reweight_leverage_flag': [0|1] for IRLS, use "leverage"
%     {default = 1}
%  'reweight_leverage_max_flag': [0|1] for IRLS, take max leverage across sensors
%     otherwise take mean
%     {default = 1}
%
% Optional Parameters to specify which stimulus conditions to use:
%  'fstem_conds': stem of csv file containing condition information
%    {default = 'cond_info'}
%  'cond_offsets_flag': [0|1|2] whether to use cond_info containing
%    condition specific offsets (i.e. best fit offsets)
%    If 1, 'condoffsets' suffix is attached to 'fstem_conds'
%    If polarity_penalty does not equal 0, a string like 'polpenalty10.0'
%      is also attached
%    If 2, 'EEG' suffix is also attached to 'fstem_conds'
%    {default = 0}
%  'conditions': vector of condition numbers to analyze
%     if empty, use all conditions found in avg_data.averages
%     if supplied, hemivec, uplowvec, etc. still apply to further limit
%     {default = []}
%  'hemivec':  vector of hemifield indices (e.g. [1,2] or [2])
%    1 = right hemifield stimuli
%    2 = left hemifield stimuli
%    {default = [1 2]}
%  'uplowvec': vector of upper or lower field indices (e.g. [1,2] or [2])
%    1 = upper field stimuli
%    2 = lower field stimuli
%    {default = [1 2]}
%  'eccvec': vector of "ecc" (eccentricity) level indices (e.g. [1,2] or [2])
%    1 = eccentricity level nearest to fovea
%    N = eccentricity level farthest from fovea
%     (with N = the number of unique eccentricity levels in fname_conds)
%    If empty, use all available eccentricity levels
%    {default = []}
%  'thetavec': vector of "theta" (polar angle) level indices (e.g. [1,6,7,12])
%    1 = polar angle nearest to but greater than 0
%    N = polar angle closest to 360
%     (with N = the number of unique polar angles in fname_conds)
%    If empty, use all available polar angles
%    {default = []}
%  'contvec': vector of contrast level indices (e.g. [1,2,3])
%    If empty, use all available contrast levels
%    {default = []}
%  'sfvec': vector of spatial frequency level indices (e.g. [1,2,3])
%    If empty, use all available spatial frequency levels
%    {default = []}
%  'tfvec': vector of temporal frequency level indices (e.g. [1,2,3])
%    If empty, use all available temporal frequency levels
%    {default = []}
%  'colvec': vector of color type indices (e.g.[1,2,3])
%    If empty, use all available color types
%    {default = []}
%
% Created:  02/06/09 by Don Hagler
% Last Mod: 03/09/15 by Don Hagler
%

%% todo: change plot_areas and plot_areas_over to check for relative path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix','RCSE',[],...
  'infix',[],[],...
  'infix_flag',false,[false true],...
  'rootoutdir',pwd,[],...
  'datamatfile','proc_avg_data_subnull.mat',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'plotflag',0,[0 1 2 3],...
  'transfile','mri2head.trans',[],...
  'badchanfile',[],[],...
  'forward_matfile',[],[],...
  'indy_locs_flag',true,[false true],...
  'graymid_flag',false,[false true],...
  'calc_dipinfo_flag',false,[false true],...
  'inverse_type',1,[0:2],...
  'ncov_type',1,[0,1,2],...
  'SNR',10,[eps Inf],...
  'verbose',true,[false true],...
  'forceflag',0,[0:2],...
... % retfit
  'retfit_dir','retfit',[],...
  'retfit_stem','retfit',[],...
  'vf2ctx_flag',true,[false true],...
  'stim_type',2,[0:2],...
  'use_areas',[],[],...
  'r_offset',0,[-100,100],...
  'th_offset',0,[-180,180],...
  'r_max',12.5,[0,Inf],...
  'ecc_width',1,[0,100],...
  'theta_width',10,[0,360],...
  'w_thresh',0.01,[0,100],...
  'vfnorm_flag',true,[false true],...
  'w_thresh_patch',0,[0,1],...
  'single_vertex_flag',0,[0 1 2],...
  'rf_sizes',[1,1,1],[0.01,10],...
  'rf_slopes',[0.1,0.1,0.1],[0,10],...
  'rf_r_inter',6,[0,Inf],...
  'restrict_hemi_flag',false,[false true],...
  'restrict_uplow_flag',false,[false true],...
  'restrict_flag',false,[false true],...
  'retfit_data_flag',false,[false true],...
... % conditions
  'fstem_conds','cond_info',[],...
  'cond_offsets_flag',0,[0 1 2],...
  'conditions',[],[],...
  'hemivec',[1,2],[1,1,2],...
  'uplowvec',[1,2],[1,1,2],...
  'eccvec',[],[],...
  'thetavec',[],[],...
  'contvec',[],[],...
  'sfvec',[],[],...
  'tfvec',[],[],...
  'colvec',[],[],...
... % fitting
  'polarity_penalty',0,[-Inf,Inf],...
  'fit_range_flag',false,[false true],...
  'fit_time0',70,[-Inf,Inf],...
  'fit_time1',100,[-Inf,Inf],...
  'offset_niters',0,[0,Inf],...
  'max_offset_niters',Inf,[0,Inf],...
  'offset_mstarts',0,[0,Inf],...
  'mstart_rotsearch_flag',false,[false true],...
  'mstart_rfsearch_flag',false,[false true],...
  'mstart_randstart_flag',true,[false true],...
  'mstart_average_flag',true,[false true],...
  'reweight_flag',false,[false true],...
  'reweight_factor',2,[],...
  'reweight_maxiter',100,[],...
  'reweight_tol',1e-7,[],...
  'reweight_leverage_flag',true,[false true],...
  'reweight_leverage_max_flag',true,[false true],...
  'hybrid_flag',false,[false true],...
  'hybrid_pre_niters',0,[0,Inf],...
  'offset_niters_last',0,[0,Inf],...
  'offset_group_flag',1,[0 1 2 3 4],...
  'offset_group_patches_flag',false,[false true],...
  'offset_group_areas_flag',false,[false true],...
  'offset_const_areas',[],[],...
  'r_step',0.1,[1e-4,10],...
  'th_step',2,[1e-4,90],...
  'r_offset_range',0,[-100,100],...
  'th_offset_range',0,[-180,180],...
  'grid_offset_flag',false,[false true],...
  'grid_offset_smooth',0,[0,1000],...
  'rot_niters',0,[0,Inf],...
  'rot_step',1,[0,180],...
  'rot_max',10,[0,180],...
  'best_rot',[0 0 0],[-180 180],...
  'rscale_niters',0,[0,Inf],...
  'rscale_step',0.01,[0,1],...
  'rscale_range',[0.8,1.2],[0.1,10],...
  'best_rscale',1,[0.1,10],...
  'nbrhd_niters',0,[0,Inf],...
  'nbrhd_full_flag',false,[false true],...
  'prior_prefix',[],[],...
  'prior_weight',0.5,[0,1],...
  'prior_avg_flag',0,[0:3],...
  'prior_mindiff_flag',0,[0 1 2],...
  'prior_zscore_flag',false,[false true],...
  'prior_sem_min',0.5,[eps,100],...
  'prior_indy_wform_flag',false,[false true],...
  'corr_time0',0,[-Inf,Inf],...
  'corr_time1',170,[-Inf,Inf],...
  'best_retmap_prefix',[],[],...
  'err_prefix',[],[],...
  'rf_niters',0,[0,Inf],...
  'rf_min_sizes',[0.1,0.1,0.1],[0.01,10],...
  'rf_max_sizes',[3,3,3],[0.01,10],...
  'rf_min_slopes',[0,0,0],[0,10],...
  'rf_max_slopes',[0.2,0.2,0.2],[0,10],...
  'rf_size_step',0.1,[0.01,10],...
  'rf_slope_step',0.01,[0.001,1],...
... % undocumented
  'calc_scalefacts_flag',false,[false true],...
  'bem_flag',true,[false true],...
  'openmeeg_flag',false,[false true],...
  'conductivities',[0.3 0.012 0.3],[],...
  'EEG_gain_scalefact',1,[-Inf Inf],...
  'nlayers',3,[1:4],...
  'bem_surf_files',...
    {'bem/inner_skull.tri','bem/outer_skull.tri','bem/outer_scalp.tri'},[],...
  'grad_scalefact',10^13,[-Inf Inf],...
  'mag_scalefact',10^15,[-Inf Inf],...
  'EEG_scalefact',10^6,[-Inf Inf],...
  'usegrad_flag',true,[false true],...
  'usemag_flag',false,[false true],...
  'useEEG_flag',false,[false true],...
  'indy_smfact',0.999,[0,1],...
  'ecc_smfact',0,[0,1],...
  'upperlower_smfact',0,[0,1],...
  'hemi_smfact',0,[0,1],...
  'loose_flag',false,[false true],...
  'loose_tang_weight',0.5,[0,1],...
  'condF_thresh',0,[],...
  'dip_matfile',[],[],...
  'ret_dips_lh_dec_file',[],[],...
  'ret_dips_rh_dec_file',[],[],...
  'ret_dip_qfield_flag',true,[false true],...
  'nonret_dips_lh_dec_file',[],[],...
  'nonret_dips_rh_dec_file',[],[],...
  'ret_dips_weight',0,[0,1],...
  'nonret_dips_weight',0,[0,1],...
  'norm_weights_flag',2,[0 1 2 3],...
  'area_names',{'v1','v2','v3'},[],...
  'area_colors',{'b','g','r'},[],...
  'plot_outdir','RCSE_plots',[],...
  'plot_normflag',false,[false true],...
  'plot_ylim',[-30,30],[-Inf,Inf],...
  'save_avg_flag',true,[false true],...
  'write_err_flag',false,[false true],...
  'write_fit_flag',false,[false true],...
  'write_areafit_flag',false,[false true],...
  'write_fif_flag',false,[false true],...
  'template_fif',[],[],...
  'nnbrs',10,[1,1000],...
  'fmincon_flag',false,[false true],...
});

% local parameters not to be passed to rc_RCSE_func
local_tags = {'best_retmap_prefix' ...
  'infix' 'infix_flag' 'cond_offsets_flag' 'datamatfile'...
  'uplowvec' 'hemivec' 'eccvec' 'thetavec' 'contvec' 'sfvec' 'tfvec' 'colvec' ...
  'err_prefix' 'fstem_conds' 'graymid_flag' 'max_offset_niters' ...
  'plot_normflag' 'plot_outdir' 'prefix' 'r_offset' 'restrict_flag' ...
  'retfit_dir' 'retfit_stem' 'th_offset'};

% parameters to be set in parms (and passed to rc_RCSE_func)
extra_tags = {'best_retmap' 'best_rf_sizes' 'best_rf_slopes' 'cond_info' ...
  'event_codes' 'lh_dip_file' 'lh_dip_info' 'lh_verts'...
  'min_error' 'rh_dip_file' 'rh_dip_info' 'rh_verts' 'surfname'...
  'total_nbrhd_niters' 'total_offset_mstarts' 'total_offset_niters'...
  'total_rf_niters' 'total_rot_niters' 'total_rscale_niters'};

% combine extra tags and what is left after excluding local tags
tags = union(setdiff(fieldnames(parms),local_tags),extra_tags);

if ~isempty(parms.badchanfile)
  if mmil_isrelative(parms.badchanfile)
    parms.badchanfile = [parms.rootoutdir '/' parms.badchanfile];
  end;
  if ~exist(parms.badchanfile,'file')
    error('file %s not found',parms.badchanfile);
  end;
end;

if mmil_isrelative(parms.transfile)
  parms.transfile = [parms.rootoutdir '/' parms.transfile];
end;
if ~exist(parms.transfile,'file')
  error('file %s not found',parms.transfile);
end;

if mmil_isrelative(parms.retfit_stem)
  if mmil_isrelative(parms.retfit_dir)
    parms.retfit_dir = [parms.rootoutdir '/' parms.retfit_dir];
  end;
  parms.retfit_stem = [parms.retfit_dir '/matfiles/' parms.retfit_stem];
end;

if parms.restrict_flag
  parms.restrict_hemi_flag = 1;
  parms.restrict_uplow_flag = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parms

if ~isempty(parms.err_prefix)
  if mmil_isrelative(parms.err_prefix)
    parms.datamatfile = sprintf('%s/matfiles/%s_err.mat',...
      parms.rootoutdir,parms.err_prefix);
  else
    parms.datamatfile = sprintf('%s_err.mat',parms.err_prefix);
  end;
else
  if mmil_isrelative(parms.datamatfile)
    parms.datamatfile = sprintf('%s/matfiles/%s',...
      parms.rootoutdir,parms.datamatfile);
  end;
end;

if ~exist(parms.datamatfile)
  error('data mat file %s not found',parms.datamatfile);
end;

if any(parms.hemivec~=[1 2])
  parms.hemivec = unique(parms.hemivec);
end;
if any(parms.uplowvec~=[1 2])
  parms.uplowvec = unique(parms.uplowvec);
end;

if mmil_isrelative(parms.fstem_conds)
  parms.fstem_conds = [parms.rootoutdir '/' parms.fstem_conds];
end;
if ~parms.cond_offsets_flag
  parms.fname_conds = [parms.fstem_conds '.csv'];
else
  if parms.polarity_penalty~=0
    parms.fname_conds = sprintf('%s_polpenalty%0.0f',...
      parms.fname_conds,parms.polarity_penalty);
  end;
  if parms.cond_offsets_flag==2
    parms.fname_conds = [parms.fname_conds '_EEG.csv'];
  else
    parms.fname_conds = [parms.fname_conds '.csv'];
  end;
end;
if ~exist(parms.fname_conds,'file')
  error('condition info file %s not found',parms.fname_conds);
end;

parms.nareas = length(parms.area_names);
parms.cond_info = rc_read_cond_info(parms.fname_conds);
nconds = length(parms.cond_info);
if nconds==0
  error('no conditions found in condition info file %s',parms.fname_conds);
end;

if parms.graymid_flag
  parms.lh_dip_file = 'bem/lh_graymid.dip';
  parms.rh_dip_file = 'bem/rh_graymid.dip';
  parms.surfname = 'graymid';
else
  parms.lh_dip_file = 'bem/lh_white.dip';
  parms.rh_dip_file = 'bem/rh_white.dip';
  parms.surfname = 'white';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find unique locations from cond_info
same_location_conds = [1:nconds]';

if isfield(parms.cond_info,'contrast')
  contrasts = [parms.cond_info.contrast];
else
  contrasts = ones(1,nconds);
end;
ind_stim = find(contrasts>0); % all except null condition (or excluded)
uniq_contrasts = unique(contrasts(ind_stim));
ncontrasts = length(uniq_contrasts);

if isfield(parms.cond_info,'spatfreq')
  spatfreqs = [parms.cond_info.spatfreq];
elseif isfield(parms.cond_info,'sfreq')
  spatfreqs = [parms.cond_info.sfreq];
else
  spatfreqs = ones(1,nconds);
end;
uniq_spatfreqs = unique(spatfreqs(ind_stim));
nspatfreqs = length(uniq_spatfreqs);

if isfield(parms.cond_info,'tempfreq')
  tempfreqs = [parms.cond_info.tempfreq];
elseif isfield(parms.cond_info,'tfreq')
  tempfreqs = [parms.cond_info.tfreq];
else
  tempfreqs = ones(1,nconds);
end;
uniq_tempfreqs = unique(tempfreqs(ind_stim));
ntempfreqs = length(uniq_tempfreqs);

if isfield(parms.cond_info,'coltype')
  coltypes = [parms.cond_info.coltype];
else
  coltypes = ones(1,nconds);
end;
uniq_coltypes = unique(coltypes(ind_stim));
ncoltypes = length(uniq_coltypes);

if ncontrasts>1 || nspatfreqs>1 || ntempfreqs>1 || ncoltypes>1
  for i=nconds:-1:2
    r = parms.cond_info(i).ecc;
    th = parms.cond_info(i).theta;
    for j=1:i-1
      tmp_r = parms.cond_info(j).ecc;
      tmp_th = parms.cond_info(j).theta;
      if tmp_r==r & tmp_th==th
        same_location_conds(i) = j;
        break;
      end;
    end;
  end;
end;
uniq_location_conds = unique(same_location_conds);

% set location number for each cond, save in cond_info
for i=1:nconds
  parms.cond_info(i).loc_num = ...
    find(uniq_location_conds==...
         same_location_conds(i));
end;

% decide which conditions to use
if isempty(parms.conditions), parms.conditions = [1:nconds]; end;

hemi_conds = [];
for i=1:nconds
  for j=parms.hemivec
    th = parms.cond_info(i).theta;
    if j==1 && (th<90 || th>270) % right hemifield
      hemi_conds = [hemi_conds i];
    elseif j==2 && (th>90 & th<270) % left hemifield
      hemi_conds = [hemi_conds i];
    end;
  end;
end;

uplow_conds = [];
for i=1:nconds
  for j=parms.uplowvec
    th = parms.cond_info(i).theta;
    if j==1 && (th<180) % upper field
      uplow_conds = [uplow_conds i];
    elseif j==2 && (th>=180) % lower field
      uplow_conds = [uplow_conds i];
    end;
  end;
end;

eccs = cell2mat({parms.cond_info.ecc});
uniq_eccs = unique(eccs(ind_stim));
if isempty(parms.eccvec)
  ecc_conds = [1:nconds];
else
  ecc_conds = find(ismember(eccs,uniq_eccs(parms.eccvec)));
  if isempty(ecc_conds)
    error('value(s) in eccvec are invalid');
  end;
end;

thetas = cell2mat({parms.cond_info.theta});
uniq_thetas = unique(thetas(ind_stim));
if isempty(parms.thetavec)
  theta_conds = [1:nconds];
else
  theta_conds = find(ismember(thetas,uniq_thetas(parms.thetavec)));
  if isempty(theta_conds)
    error('value(s) in thetavec are invalid');
  end;
end;

% find conditions with allowed contrasts
if isempty(parms.contvec)
  cont_conds = ind_stim;
else
  cont_conds = ...
    find(ismember(contrasts,uniq_contrasts(parms.contvec)));
  if isempty(cont_conds)
    error('value(s) in contvec are invalid');
  end;
end;

% find conditions with allowed spatfreqs
if isempty(parms.sfvec)
  sf_conds = [1:nconds];
else
  sf_conds = ...
    find(ismember(spatfreqs,uniq_spatfreqs(parms.sfvec)));
  if isempty(sf_conds)
    error('value(s) in sfvec are invalid');
  end;
end;

% find conditions with allowed tempfreqs
if isempty(parms.tfvec)
  tf_conds = [1:nconds];
else
  tf_conds = ...
    find(ismember(tempfreqs,uniq_tempfreqs(parms.tfvec)));
  if isempty(tf_conds)
    error('value(s) in tfvec are invalid');
  end;
end;

% find conditions with allowed coltypes
if isempty(parms.colvec)
  col_conds = [1:nconds];
else
  col_conds = ...
    find(ismember(coltypes,uniq_coltypes(parms.colvec)));
  if isempty(col_conds)
    error('value(s) in colvec are invalid');
  end;
end;

parms.conditions = intersect(hemi_conds,parms.conditions);
parms.conditions = intersect(uplow_conds,parms.conditions);
parms.conditions = intersect(ecc_conds,parms.conditions);
parms.conditions = intersect(theta_conds,parms.conditions);
parms.conditions = intersect(cont_conds,parms.conditions);
parms.conditions = intersect(sf_conds,parms.conditions);
parms.conditions = intersect(tf_conds,parms.conditions);
parms.conditions = intersect(col_conds,parms.conditions);

parms.event_codes = [];

if isempty(parms.conditions)
  error('no usable conditions');
end;

if ~isempty(parms.best_retmap_prefix)
  parms.prior_prefix = [];
  parms.offset_niters = 0;
  parms.r_offset = 0;
  parms.th_offset = 0;
  parms.nbrhd_niters = 0;
  parms.rot_niters = 0;
  parms.rf_niters = 0;
end;

% restrict flags determine number of offsets
if parms.offset_niters>0
  parms.offset_infix_list = {''};
  if ~parms.restrict_hemi_flag
    parms.offset_infix_list{end+1} = '_ipsi';
  end;
  if ~parms.restrict_uplow_flag
    parms.offset_infix_list{end+1} = '_cross';
  end;
  if ~parms.restrict_hemi_flag && ~parms.restrict_uplow_flag
    parms.offset_infix_list{end+1} = '_ipsi_cross';
  end;
  parms.npatches = length(parms.offset_infix_list);
  parms.offset_types = {'r','th'};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set full output prefix
if parms.infix_flag && isempty(parms.infix)
  parms.infix = rc_RCSE_set_infix(parms);
end;
if ~isempty(parms.infix)
  parms.prefix = [parms.prefix '_' parms.infix];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.plotflag~=3
  % set options for rc_RCSE_func
  tmp_parms = parms;
  tmp_parms.plotflag = (parms.plotflag==2);
  tmp_parms.total_offset_niters = 0;
  tmp_parms.total_offset_mstarts = 0;
  tmp_parms.total_rot_niters = 0;
  tmp_parms.total_rscale_niters = 0;
  tmp_parms.total_nbrhd_niters = 0;
  tmp_parms.total_rf_niters = 0;
  tmp_parms.best_retmap = [];
  tmp_parms.min_error = [];

  % load best_retmap if best_retmap_prefix not empty
  %   or offset_niters>0  or nbrhd_niters>0
  if ~isempty(parms.best_retmap_prefix)
    if mmil_isrelative(parms.best_retmap_prefix)
      matfile = sprintf('%s/matfiles/%s_ret_forward.mat',...
        parms.rootoutdir,parms.best_retmap_prefix);
    else
      matfile = sprintf('%s_ret_forward.mat',parms.best_retmap_prefix);
    end;
    if ~exist(matfile,'file')
      error('best_retmap file %s not found',mfilename,matfile);
    end;
    fprintf('%s: loading previous best fit from %s...\n',mfilename,matfile);
    load(matfile)
    tmp_parms.best_retmap = retforward.retmap;
    tmp_parms.min_error = retforward.min_error;
    if isfield(retforward,'total_offset_niters')
      tmp_parms.total_offset_niters = retforward.total_offset_niters;
    end;
    if isfield(retforward,'total_offset_mstarts')
      tmp_parms.total_offset_mstarts = retforward.total_offset_mstarts;
    end;
    if isfield(retforward,'total_nbrhd_niters')
      tmp_parms.total_nbrhd_niters = retforward.total_nbrhd_niters;
    end;
    fprintf('%s: current minimum error = %0.4f\n',...
      mfilename,tmp_parms.min_error);
  elseif parms.offset_niters>0 || parms.nbrhd_niters>0
    matfile = sprintf('%s/matfiles/%s_ret_forward.mat',...
      parms.rootoutdir,parms.prefix);
    if exist(matfile,'file')
      fprintf('%s: loading previous best fit from %s...\n',mfilename,matfile);
      load(matfile)
      tmp_parms.best_retmap = retforward.retmap;
      tmp_parms.min_error = retforward.min_error;
      if isfield(retforward,'total_offset_niters')
        tmp_parms.total_offset_niters = retforward.total_offset_niters;
      end;
      if isfield(retforward,'total_offset_mstarts')
        tmp_parms.total_offset_mstarts = retforward.total_offset_mstarts;
      end;
      if isfield(retforward,'total_nbrhd_niters')
        tmp_parms.total_nbrhd_niters = retforward.total_nbrhd_niters;
      end;
      if parms.offset_niters>0
        if ~isfield(tmp_parms.best_retmap,'r_offset') &...
            isfield(tmp_parms.best_retmap,'r_offsets')
          tmp_parms.best_r_offsets = tmp_parms.best_retmap.r_offsets;
        else
          tmp_parms.best_r_offsets = ...
            zeros(parms.nareas,nconds,parms.npatches);
        end;
        if ~isfield(tmp_parms.best_retmap,'th_offset') &...
            isfield(tmp_parms.best_retmap,'th_offsets')
          tmp_parms.best_th_offsets = tmp_parms.best_retmap.th_offsets;
        else
          tmp_parms.best_th_offsets = ...
            zeros(parms.nareas,nconds,parms.npatches);
        end;
        event_codes = cell2mat({tmp_parms.best_retmap.orig_cond_info.event_code});
        % set offsets in cond_info
        for a=1:parms.nareas
          for i=1:nconds
            [usable,j] = ismember(tmp_parms.cond_info(i).event_code,event_codes);
            if ~usable, continue; end;
            for t=1:length(parms.offset_types)
              for p=1:parms.npatches
                fieldname = sprintf('%s_offset_%s%s',...
                  parms.offset_types{t},parms.area_names{a},...
                  parms.offset_infix_list{p});
                source_fieldname = ['best_' parms.offset_types{t} '_offsets'];
                tmp_parms.cond_info = setfield(tmp_parms.cond_info,{i},fieldname,{1},...
                  tmp_parms.(source_fieldname)(a,j,p));
              end;
            end;
          end;
        end;
      end;
      fprintf('%s: current minimum error = %0.4f\n',...
        mfilename,tmp_parms.min_error);
    else
      fprintf('%s: initializing best fit (%s not found)...\n',mfilename,matfile);
    end;
  end;

  if parms.rot_niters>0
    matfile = sprintf('%s/matfiles/%s_ret_forward.mat',...
      parms.rootoutdir,parms.prefix);
    if exist(matfile,'file')
      load(matfile);
      tmp_parms.best_rot = retforward.best_rot;
      tmp_parms.min_error = retforward.min_error;
      if isfield(retforward,'total_rot_niters')
        tmp_parms.total_rot_niters = retforward.total_rot_niters;
      end;
    end;
  end;

  if parms.rscale_niters>0
    matfile = sprintf('%s/matfiles/%s_ret_forward.mat',...
      parms.rootoutdir,parms.prefix);
    if exist(matfile,'file')
      load(matfile);
      tmp_parms.best_rscale = retforward.best_rscale;
      tmp_parms.min_error = retforward.min_error;
      if isfield(retforward,'total_rscale_niters')
        tmp_parms.total_rscale_niters = retforward.total_rscale_niters;
      end;
    end;
  end;

  best_rf_sizes = [];
  best_rf_slopes = [];
  if parms.rf_niters>0
    matfile = sprintf('%s/matfiles/%s_ret_forward.mat',...
      parms.rootoutdir,parms.prefix);
    if exist(matfile,'file')
      load(matfile);
      tmp_parms.best_rf_sizes = retforward.best_rf_sizes;
      tmp_parms.best_rf_slopes = retforward.best_rf_slopes;
      tmp_parms.min_error = retforward.min_error;
      if isfield(retforward,'total_rf_niters')
        tmp_parms.total_rf_niters = retforward.total_rf_niters;
      end;
    end;
  end;

  if tmp_parms.total_offset_niters >= parms.max_offset_niters
    fprintf('%s: number of offset iterations (%d) exceeds max_offset_niters (%d)\n',...
      mfilename,tmp_parms.total_offset_niters,parms.max_offset_niters);
    return;
  elseif tmp_parms.total_offset_niters + parms.offset_niters >...
         parms.max_offset_niters + 1
    tmp_parms.offset_niters = ...
      parms.max_offset_niters - tmp_parms.total_offset_niters + 1;
  end;

  % load retfit files
  retfit_results = rc_load_retfit_results(parms.retfit_stem);
  tmp_parms.lh_verts = find(retfit_results.lh_fit_data.th~=0 |...
                        retfit_results.lh_fit_data.r~=0);
  tmp_parms.rh_verts = find(retfit_results.rh_fit_data.th~=0 |...
                        retfit_results.rh_fit_data.r~=0);

  % load data
  fprintf('%s: loading data from %s...\n',mfilename,parms.datamatfile);
  load(parms.datamatfile);
  if ~isempty(parms.err_prefix)
    avg_data = err_data;
  end;
  avg_data.source_fname = parms.datamatfile;

  % run rc_RCSE_func
  start_time = cputime;
  args = mmil_parms2args(tmp_parms,tags);
  if ~isempty(parms.best_retmap_prefix)
    rc_RCSE_func(avg_data,subj,...
      'retmap',tmp_parms.best_retmap,...
      'prefix',parms.prefix,args{:});
  else
    rc_RCSE_func(avg_data,subj,...
      'retfit_results',retfit_results,...
      'r_offset',parms.r_offset,'th_offset',parms.th_offset,...
      'prefix',parms.prefix,args{:});
  end;
  elapsed_time = cputime - start_time;
  fprintf('%s: total elapsed time for RCSE = %0.1f seconds\n',...
    mfilename,elapsed_time);

  if parms.offset_niters>0
    % save best offsets in csv file
    [tmp,tmp_fstem] = fileparts(parms.fstem_conds);
    % NOTE: fstem_conds may be from other location, but fname_out
    %   should be in rootoutdir
    fname_out = [parms.rootoutdir '/' parms.prefix '_' tmp_fstem];
    if parms.grid_offset_flag
      fname_out = [fname_out '_offsets'];
    else
      fname_out = [fname_out '_rthoffsets'];
    end;      
    if ~isempty(parms.prior_prefix)
      fname_out = [fname_out '_prior'];
    end;
    if parms.polarity_penalty~=0
      fname_out = sprintf('%s_polpenalty%0.0f',...
        fname_out,parms.polarity_penalty);
    end;
    fname_out = [fname_out '.csv'];
    matfile = sprintf('%s/matfiles/%s_results.mat',...
      parms.rootoutdir,parms.prefix);
    load(matfile);
    cond_info = parms.cond_info;
    tmp_cond_info = results.retmap.cond_info;
    tmp_loc_nums = cell2mat({tmp_cond_info.loc_num});
    if isempty(parms.use_areas)
      use_areas = [1:parms.nareas];
    else
      use_areas = parms.use_areas;
    end;
    for i=1:length(cond_info)
      % find condition in tmp_cond_info with same stim location
      loc_num = cond_info(i).loc_num;
      j = find(loc_num==tmp_loc_nums);
      if isempty(j), continue; end;
      for a=use_areas
        for t=1:length(parms.offset_types)
          for p=1:parms.npatches
            fieldname = sprintf('%s_offset_%s%s',...
              parms.offset_types{t},parms.area_names{a},...
              parms.offset_infix_list{p});
            if isfield(tmp_cond_info,fieldname)
              val = tmp_cond_info(j).(fieldname);
            else
              val = 0;
            end;
            cond_info(i).(fieldname) = val;
          end;
        end;
      end;
    end;    
    rc_write_cond_info(cond_info,fname_out,1);
  end;
end;

if ismember(parms.plotflag,[1,3])
  cwd = pwd;
  cd(parms.rootoutdir);
  if mmil_isrelative(parms.plot_outdir)
    parms.plot_outdir = [parms.rootoutdir '/' parms.plot_outdir];
  end;
  mmil_mkdir(parms.plot_outdir);
  figure;
  rc_plot_sources('prefix',parms.prefix,...
    'area_names',parms.area_names,'area_colors',parms.area_colors,...
    'ylim',parms.plot_ylim,...
    'rootdir',parms.rootoutdir,'normflag',parms.plot_normflag);
  set(gcf,'Visible','off');
  fname_areas = [parms.plot_outdir '/' parms.prefix '_areas.tif'];
  print('-dtiff',fname_areas);
  close(gcf);
  figure;
  fname_err = rc_plot_fitvar(parms.prefix);
  set(gcf,'Visible','off');
  print('-dtiff',fname_err);
  close(gcf);

  % create montage
  fname_montage = sprintf('%s/%s_montage.tif',parms.plot_outdir,parms.prefix);
  img_areas = imread(fname_areas);
  img_err = imread(fname_err);
  img_mont = cat(2,img_areas,img_err);
  imwrite(img_mont,fname_montage);

  % remove originals
  cmd = sprintf('rm %s %s',fname_areas,fname_err);
  [status,result] = unix(cmd);
  if status
    error('failed to cleanup tif files:\n%s',result);
  end;

  cd(cwd);
end;

prefix = parms.prefix;

return;

