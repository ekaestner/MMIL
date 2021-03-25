function ts_dSPM(avg_data,subjname,varargin)
%function ts_dSPM(avg_data,subjname,[options])
%
% Usage:
%  ts_dSPM(avg_data,subjname, 'key1', value1,...)
%
% Required Input:
%  avg_data: average data structure
%    (see ts_process_fif_data and ts_avg_fif_data)
%  subjname: freesurfer subject name (must be found in $SUBJECTS_DIR)
%
% Optional parameters:
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'prefix': prefix of all output files
%    {default = 'dSPM'}
%  'rootoutdir': root output directory (several subdirectories will be created)
%    {default = pwd} (current working directory)
%  'conditions': vector of condition numbers to analyze
%    {default = []} (if empty, use all conditions found in avg_data.averages)
%  'badchans'   : vector of bad channel indices
%    {default = []}
%  'badchanfile': name of text file containing bad channel labels
%    {default = []}
%  'usegrad_flag': [0|1] use gradiometer data if available
%    {default = 1}
%  'usemag_flag': [0|1] use magnetometer data if available
%    {default = 1}
%  'useEEG_flag': [0|1] use EEG data if available
%    {default = 1}
%  'baseline_flag': [0|1] baseline correction of input data
%    {default = 1}
%  'baseline_start': start time of baseline period (msec)
%    if ncov_type = 1, use this time window for baseline correction
%      and noise covariance and scaling factor calculation from average data
%    relative to trigger onset; negative times occur before trigger
%    {default = -Inf} (start at beginning of prestimulus period)
%  'baseline_end'  : end time of baseline period (msec)
%    {default = 0} (end at trigger onset)
%  'ssp_projmat': Signal Space Projection matrix (nchan x nchan)
%    will be applied to both the input data and the forward gain matrix
%    {default = []}
%  'verbose': [0|1] display status messages
%     {default = 1}
%  'forceflag': [0|1|2] whether to overwrite files
%     0: do nothing if output exists
%     1: overwrite results and other output
%     2: overwrite all output, including forward solution
%     {default = 0}
%
% Optional parameters for forward model:
%  'forward_matfile': mat file containing gain matrix
%     This allows user to supply precalculated gain matrix from an earlier
%       run of this program -- for example, calculate gain matrix with
%       all channels (grad, mag, and EEG) and then rerun excluding certain
%       channels or channel types
%     Do not use this option if you change the number of dipoles
%     {default = []}
%  'calc_dipinfo_flag': [0|1] whether to calculate dipole info from surfaces
%     Otherwise use lh_dip_info and rh_dip_info or load from .dip files
%     {default = 1}
%  'lh_dip_info': matrix of left hemi dipole information (6 x ndips)
%    6 columns for x,y,z coordinates and nx,ny,nz normal vector
%    if specified, will ignore lh_dip_file
%    {default = []}
%  'rh_dip_info': matrix of right hemi dipole information (6 x ndips)
%    6 columns for x,y,z coordinates and nx,ny,nz normal vector
%    if specified, will ignore rh_dip_file
%    {default = []}
%  'lh_dip_file': name of left hemi freesurfer dipole file (from tksurfer)
%    {default = 'bem/lh_white.dip'} (relative to subject dir)
%  'rh_dip_file': name of right hemi freesurfer dipole file (from tksurfer)
%    {default = 'bem/rh_white.dip'} (relative to subject dir)
%  'nodec_flag': [0|1] if 1, use all dipoles; otherwise use only those
%     specified by lh_dec_dips and rh_dec_dips (or lh_dec_file and rh_dec_file)
%     {default = 0}
%  'lh_dec_dips': vector of 0's and 1's indicating subset of left hemi dipoles
%     length must match ndips (nvertices)
%     if specified, will ignore lh_dec_file
%     {default = []}
%  'rh_dec_dips': vector of 0's and 1's indicating subset of right hemi dipoles
%     if specified, will ignore rh_dec_file
%     {default = []}
%  'lh_dec_file': name of left hemi freesurfer dec file (from tksurfer)
%     {default = 'bem/lh_white_7.dec'} (relative to subject dir)
%     dec file contains 0's and 1's indicating subset of dipoles to use
%  'rh_dec_file': name of right hemi freesurfer dec file (from tksurfer)
%     {default = 'bem/rh_white_7.dec'} (relative to subject dir)
%  'bem_flag': [0|1] use boundary element method (BEM)
%     if 0, use spherical shell model instead
%     {default = 1}
%  'openmeeg_flag': [0|1] use OpenMEEG to calculate BEM solution
%    otherwise use M.X. Huang's BEM solution function
%    {default = 0}
%  'nlayers': number of layers for BEM or spherical shell model
%     For MEG only, can be 1 or 3; with EEG (alone or with MEG), must be 3
%     {default = 3}
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
%  'cen_sph': center of sphere (mm): vector of 3 numbers (x,y,z)
%     ignored if bem_flag=1
%     {default = []}
%  'radii': vector of three values specifying radii (mm) of spherical shell layers:
%     (1) inner skull (2) outer skull (3) outer scalp
%     Required if useEEG_flag=1 and bem_flag=0 (spherical shell model)
%     {default = []}
%  'conductivities': vector of three values specifying conductivity of:
%     (1) brain (2) skull  (3) scalp
%     {default = [0.3 0.012 0.3]}
%  'EEG_gain_scalefact': scaling factor applied to EEG gain matrix
%     an optimized value may be necessary for integration of MEG and EEG data
%     {default = 1}
%  'refEEG_coords': reference EEG electrode coordinates (e.g. [0.004,0.020,0.001])
%     if specified, forward solution for these coordinates is subtracted 
%       from forward solution for all other EEG electrodes
%     for data originally in fif format, these coordinates are already
%       stored in the "loc" matrix for each EEG electrode
%     {default = []}
%  'forward_only_flag': [0|1] calculate forward solution and quit
%     {default = 0}
%
% Optional parameters for coordinate transformation:
%  'trans': 4x4 matrix specifying mri2head transformation
%     {default = []}
%  'transfile': text file containing 4x4 mri2head transformation
%     {default = []}
%  'alignment_fif' fif file containing 4x4 mri2head transformation
%     note on reading from fif's:
%      This requires Uutela's fifaccess toolbox which only runs on 32 bit Matlab
%      and core dump will result if mri2head transform is missing from
%      specified file
%     {default = []}
%
% Note on transform files:
%   trans, trans_file, and alignment_fiff are different, mutually exclusive
%     ways to specify the mri2head transformation matrix.
%   This transformation matrix is a 4x4 afine transformation that
%     maps locations in the MRI (from dip files) to locations in the
%     MEG/EEG head space defined by the digitized cardinal points
%     (and HPI coils).
%   This transformation matrix can also be found in the avg_data structure.
%     The avg_data structure contains a field called coor_trans.
%     If the subfield "mri2head" exists and is not empty, this will be used.
%   If more than one transformation marix is supplied, the order of dominance is:
%     trans, trans_file, alignment_fiff, avg_data.coor_trans.mri2head
%     i.e., if trans is supplied directly, all other methods of
%     supplying the trans will be ignored.
%
% Optional parameters for noise covariance:
%  'ncov_type': specify what type of noise covariance matrix to use
%     if 0, use identity matrix
%       (assume uniform white noise, independently scaled for each sensor type)
%     if 1, calculate and use noise covariance matrix from average prestim
%     if 2, use noise covariance matrix calculated from single trials
%       (stored in avg_data.noise.covar)
%    {default = 2}
%  'ncov_conditions': vector of condition numbers to use for noise covariance
%    matrix calculation (if empty, will use all in 'conditions')
%    {default = []}
%  'raw_ncov_lambda': regularization factor used if ncov_type = 2
%    {default = 0.1}
%
% Optional parameters for fMRI bias:
%  'lh_sourcecov_file': full pathname of left hemi freesurfer mgh file
%     file can contain fMRI activations on surface for a priori weighting
%     if not supplied, will used identity matrix as source covariance
%     {default = []}
%  'rh_sourcecov_file': full pathname of right hemi freesurfer mgh file
%     {default = []}
%  'sourcecov_thresh_flag': [0|1] binarize input values with a threshold
%       and set values for subthreshold and suprathreshold elements
%     otherwise, use input values as they are (scaled so max is maxvar)
%    {default = 1}
%  'sourcecov_thresh': threshold applied to values in sourcecov mgh files
%     suprathreshold vertices (dipoles) will get weighting of sourcecov_maxvar
%     subthreshold vertices (dipoles) will get weighting of sourcecov_minvar
%     {default = 0}
%  'sourcecov_thresh_abs_flag': [0|1] use absolute value threshold
%     {default = 1}
%  'sourcecov_maxvar': value assigned to diagonal element for suprathreshold
%      dipoles in source covariance matrix
%     {default = 0.9}
%  'sourcecov_minvar': value assigned to diagonal element for suprathreshold
%      dipoles in source covariance matrix
%     {default = 0.09}
%
% Optional parameters for smoothness constraint:
%  'smooth_constr_flag': [0|1] use smoothness constraint
%     If 1, forces orient_constr_flag = 1
%     If 1, forces nodec_flag = 1
%       To reduce computation, use subject resampled to ico<7
%     {default = 0}
%  'smooth_constr_type': string describing type of smoothness constraint
%    'smooth': weights from iterative smoothing
%    'nbrhd': weights from neighborhood distances
%    'roi_nbrhd': weights from neighborhood distances within ROIs
%     {default = 'smooth'}
%  'smooth_constr_nsmooth': number of smoothing steps for smoothness constraint
%     only applies if smooth_constr_type = 'smooth'
%     {default = 10}
%  'smooth_constr_nbrhd_a0': variance weight for each vertex (in ROI)
%     only applies if smooth_constr_type = 'nbrhd' or 'roi_nbrhd'
%     {default = 1.0}
%  'smooth_constr_nbrhd_a1': covariance weight for first degree neighbor
%     only applies if smooth_constr_type = 'nbrhd' or 'roi_nbrhd'
%     {default = 0.5}
%  'smooth_constr_nbrhd_a2': covariance weight for second degree neighbor
%     only applies if smooth_constr_type = 'nbrhd' or 'roi_nbrhd'
%     {default = 0.25}
%  'smooth_constr_nbrhd_ax': covariance weight for all other vertices in ROIs
%     only applies if smooth_constr_type = 'roi_nbrhd'
%     {default = 0}
%  'smooth_constr_nbrhd_az': variance weight for all other vertices in surface
%     only applies if smooth_constr_type = 'roi_nbrhd'
%     {default = 0.1}
%  'smooth_constr_roi_files': cell array of label file names
%     with size = [nroi,nhemi], with nhemi = 2
%     only applies and must be supplied if smooth_constr_type = 'roi_nbrhd'
%     {default = []}
%  'smooth_constr_nbrhd_outfix': string included in smooth covar file name
%     only applies if smooth_constr_type = 'nbrhd'
%     {default = 'nbrhd'}
%  'smooth_constr_roi_nbrhd_outfix': string included in smooth covar file name
%     only applies if smooth_constr_type = 'roi_nbrhd'
%     {default = 'roi_nbrhd'}
%
% Optional parameters for inverse:
%  'inverse_matfile': mat file containing inverse operator
%     Do not use this option if you change the number of dipoles or SNR
%     {default = []}
%  'orient_constr_flag': [0|1] fix dipole orientations perpendicular to cortex
%     {default = 0}
%  'orient_tang': source covariance weight of tangential components
%       when using orientation constraint
%     Value of 0 gives a fixed orientation constraint
%     Value of 1 allows free orientations (equal weighting with normal vector)
%     Intermediate value gives a "loose" orientation constraint
%     Only applicable when orient_constr_flag=1
%     {default = 0}
%  'SNR': estimated signal-to-noise-ratio (for regularization parameter)
%     {default = 10}
%  'prewhiten_flag': [0|1] use Matti Hamalainen's pre-whitened inverse operator
%     otherwise use Anders Dale's original method
%     {default = 0}
%  'calc_scalefacts_flag': [0|1] calculate scaling factors
%     for different channel types (standard deviation of avg baseline)
%     if 1, calculate these scaling factors and apply to data and forward matrix
%     if 0, use default or user specified scaling factors
%     NOTE: ignored if prewhiten_flag=1
%     {default = 0}
%  'grad_scalefact': scaling factor applied to gradiometer data
%     purpose of scaling factors is to get data from different channel types
%     into roughly the same scale
%     NOTE: ignored if prewhiten_flag=1
%     {default = 10^13}
%  'mag_scalefact': scaling factor applied to magnetometer data
%     NOTE: ignored if prewhiten_flag=1
%     {default = 10^15}
%  'EEG_scalefact': scaling factor applied to EEG data
%     NOTE: ignored if prewhiten_flag=1
%     {default = 10^6}
%  'noisenorm_flag': [0|1] Noise-normalization (Dale et al., Neuron 2000)
%     If 1, source waveforms are noise sensitivity normalized sqrt(F) or z-stats
%     If 0, source waveforms are in nA*m (but are depth-biased)
%     {default = 1}
%  'noisenorm_identity_flag': [0|1] use identify matrix for noise-normalization
%     otherwise use matrix determined by ncov_type
%     NOTE: ignored if prewhiten_flag=1
%     {default = 1}
%  'depthweight_flag': [0|1] depth-weighting (Lin et al., NeuroImage 2006)
%     If 1, source covariance is depth-weighted
%     If 0, source covariance is not depth-weighted
%     {default = 0}
%  'depthweight_p': depth weighting parameter
%     {default = 0.5}
%  'signed_sources_flag': [0|1] output of source waveforms
%       as signed amplitudes along normal vectors
%     If 1, and noisenorm_flag=1, source waveforms are z-stats
%     If 0, and noisenorm_flag=1, source waveforms are sqrt(F-stats)
%       The F-stats are unsigned because they are the sum of the
%        squared components of dipole vector
%     Only available if orient_constr_flag = 1
%     {default = 1}
%
% Optional parameters for output:
%  'write_stc_flag': [0|1] output source waveforms as stc files
%    forced to 1 if write_mgh_flag = 1
%    {default = 1}
%  'stc_scalefact': scale factor multiplier for output stc file waveforms
%    {default = 1}
%  'write_mgh_flag': [0|1] additional output of source waveforms
%    as mgh files (displayable with FreeSurfer's tksurfer)
%    forced to 1 if resamp2ico_flag = 1
%    {default = 0}
%  'sparsesmooth': # of sparse smoothing steps applied before saving mgh file
%    see ts_stc2mgh for details on smoothing
%    {default = 10}
%  'postsmooth: # of normal smoothing steps applied when saving as mgh file
%    {default = 10}
%  'mbmask_flag': [0|1] whether to mask midbrain & corpus callosum
%    {default = 0}
%  'resamp2ico_flag': [0|1] resample mgh files to icosahedral sphere
%    necessary preparation for group analysis
%    (displayable on fsaverage subject with FreeSurfer's tksurfer)
%    {default = 0}
%  'icolevel': icosahedron order number:
%               Order  Number of Vertices
%                 0              12
%                 1              42
%                 2             162
%                 3             642
%                 4            2562
%                 5           10242
%                 6           40962
%                 7          163842
%    {default = 7}
%  'icosmooth': smoothing steps on surface after sampling to ico
%    {default = 3}
%  'write_fif_flag': [0|1] output fit and residual error
%    sensor waveforms as fif files (displayable with Neuromag's xplotter)
%    {default = 0}
%  'template_fif': full path name of online average fif file
%    must specify if write_fif_flag=1
%    {default = []}
%
% Created:  08/02/06 by Don Hagler
% Last Mod: 06/16/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
if ~mmil_check_nargs(nargin,2), return; end;
parms = check_input(avg_data,subjname,varargin);

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

% save parms to mat file for future reference
save_parms(parms);

% calculate forward solution
[parms,G_xyz] = calc_forward(parms,avg_data);
if parms.forward_only_flag, return; end;

% construct source covariance matrix
parms = construct_srccovar(parms,G_xyz);

% calculate inverse operator
[M,nnf] = calc_inverse(parms,G_xyz);

% calculate source estimates
calc_source_estimates(parms,avg_data,M,nnf,G_xyz);

if parms.verbose
  fprintf('%s: finished.\n',mfilename);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(avg_data,subjname,options)
  parms = mmil_args2parms(options,{...
    'subjname',subjname,[],...
  ...
    'subjdir',[],[],...
    'rootoutdir',pwd,[],...
    'prefix','dSPM',[],...
    'conditions',[],[],...
    'calc_dipinfo_flag',true,[false true],...
    'lh_dip_info',[],[],...
    'rh_dip_info',[],[],...
    'lh_dip_file','bem/lh_white.dip',[],...
    'rh_dip_file','bem/rh_white.dip',[],...
    'nodec_flag',false,[false true],...
    'lh_dec_dips',[],[],...
    'rh_dec_dips',[],[],...
    'lh_dec_file','bem/lh_white_7.dec',[],...
    'rh_dec_file','bem/rh_white_7.dec',[],...
    'ncov_type',2,[0 2],...
    'ncov_conditions',[],[],...
    'raw_ncov_lambda',0.1,[0,1],...
    'calc_scalefacts_flag',false,[false true],...
    'baseline_start',-Inf,[-Inf Inf],...
    'baseline_end',0,[-Inf Inf],...
    'baseline_flag',true,[false true],...
    'ssp_projmat',[],[],...
    'verbose',true,[false true],...
    'forceflag',0,[0:2],...
  ...
    'SNR',10,[eps Inf],...
    'noisenorm_flag',true,[false true],...
    'noisenorm_identity_flag',true,[false true],...
    'depthweight_flag',false,[false true],...
    'depthweight_p',0.5,[0,1],...
    'bem_flag',true,[false true],...
    'openmeeg_flag',false,[false true],...
    'radii',[],[],...
    'conductivities',[0.3 0.012 0.3],[],...
    'EEG_gain_scalefact',1,[-Inf Inf],...
    'nlayers',3,[1:4],...
    'badchans',[],[],...
    'badchanfile',[],[],...
    'usegrad_flag',true,[false true],...
    'usemag_flag',true,[false true],...
    'useEEG_flag',true,[false true],...
    'grad_scalefact',10^13,[-Inf Inf],...
    'mag_scalefact',10^15,[-Inf Inf],...
    'EEG_scalefact',10^6,[-Inf Inf],...
    'write_stc_flag',true,[false true],...
    'stc_scalefact',1,[eps Inf],...
    'write_mgh_flag',false,[false true],...
    'sparsesmooth',10,[0 1000],...
    'postsmooth',10,[0 1000],...
    'mbmask_flag',false,[false,true],...
    'resamp2ico_flag',false,[false true],...
    'icolevel',7,[1 7],...
    'icosmooth',3,[0 1000],...
    'write_fif_flag',false,[false true],...
    'template_fif',[],[],...
    'bem_surf_files',...
      {'bem/inner_skull.tri','bem/outer_skull.tri','bem/outer_scalp.tri'},[],...
    'cen_sph',[],[],...
    'trans',[],[],...
    'transfile',[],[],...
    'alignment_fif',[],[],...
    'lh_sourcecov_file',[],[],...
    'rh_sourcecov_file',[],[],...
    'sourcecov_thresh_flag',true,[false true],...
    'sourcecov_thresh',0,[-Inf Inf],...
    'sourcecov_thresh_abs_flag',true,[false true],...
    'sourcecov_maxvar',0.9,[0 1],...
    'sourcecov_minvar',0.09,[0 1],...
    'forward_matfile',[],[],...
    'inverse_matfile',[],[],...
    'refEEG_coords',[],[],...
    'forward_only_flag',false,[false true],...
    'prewhiten_flag',false,[false true],...
    'orient_constr_flag',false,[false true],...
    'orient_tang',0,[0 1],...
    'smooth_constr_flag',false,[false true],...
    'smooth_constr_type','smooth',{'smooth','nbrhd','roi_nbrhd'},...
    'smooth_constr_nsmooth',10,[1,100],...
    'smooth_constr_nbrhd_a0',1.0,[0.01,100.0],...
    'smooth_constr_nbrhd_a1',0.5,[0.01,100.0],...
    'smooth_constr_nbrhd_a2',0.25,[0.01,100.0],...
    'smooth_constr_nbrhd_ax',0,[0,100.0],...
    'smooth_constr_nbrhd_az',0.1,[0,100.0],...
    'smooth_constr_roi_files',[],[],...
    'smooth_constr_nbrhd_outfix','nbrhd',[],...
    'smooth_constr_roi_nbrhd_outfix','roi_nbrhd',[],...
    'signed_sources_flag',true,[false true],...
  ... % hidden
    'datatype','single',{'single','double'},...
    'save_results_flag',true,[false true],...
    'save_avg_flag',false,[false true],...
    'save_fiterr_flag',false,[false true],...
    'save_fiterr_struct_flag',false,[false true],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
  ...  % scaling factors to get forward matrix into same scale as fif data
    'grad_unitsfact',10^-13,[],... % fT/cm to T/m
    'mag_unitsfact',10^-15,[],... % fT to T
    'EEG_unitsfact',10^-6,[],... % uV to V
  ...
    'chan_tags',{'badchans','badchanfile','usegrad_flag','usemag_flag',...
                 'useEEG_flag','channames','verbose',},[],...
    'gainmat_tags',{'fname_lh_surf','fname_rh_surf',...
                    'lh_dip_info','rh_dip_info','trans','bem_surf_files',...
                    'bem_matfile','lh_dec_dips','rh_dec_dips',...
                    'bem_flag','openmeeg_flag','useMEG_flag','useEEG_flag',...
                    'conductivities','cen_sph','radii',...
                    'rootoutdir','verbose','forceflag','batchsize'},[],...
  });

  % check input data
  avg_data  = ts_checkdata_header(avg_data,'conditions',parms.conditions);

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

  % conditions to use
  if isempty(parms.conditions)
    parms.conditions = 1:length(avg_data.averages);
  else
    conditions = 1:length(avg_data.averages);
    parms.conditions = intersect(parms.conditions,conditions);
  end;
  if isempty(parms.conditions)
    error('no valid conditions specified');
  end;
  if isempty(parms.ncov_conditions)
    parms.ncov_conditions = parms.conditions;
  else
    conditions = 1:length(avg_data.averages);
    parms.ncov_conditions = intersect(parms.ncov_conditions,conditions);
  end;
  if isempty(parms.ncov_conditions)
    error('no valid ncov conditions specified');
  end;

  % mri2head transformation matrix
  T_mri2head = [];
  if ~isempty(parms.trans)
    if parms.verbose
      fprintf('%s: using transformation matrix supplied...\n',mfilename);
    end;
    T_mri2head = parms.trans;
  elseif ~isempty(parms.transfile)
    if ~exist(parms.transfile,'file')
      error('transfile %s not found',parms.transfile);
    end;
    if parms.verbose
      fprintf('%s: reading transformation matrix from %s...\n',...
        mfilename,parms.transfile);
    end;
    T_mri2head = ts_read_transfile(parms.transfile);
  elseif ~isempty(parms.alignment_fif)
    if ~exist(parms.alignment_fif,'file')
      error('alignment_fif %s not found',parms.alignment_fif);
    end
    if parms.verbose
      fprintf('%s: reading transformation matrix from %s...\n',...
        mfilename,parms.alignment_fif);
    end;
    T_mri2head=loadtrans(parms.alignment_fif,'MRI','HEAD');
  elseif isfield(avg_data.coor_trans,'mri2head')
    if parms.verbose
      fprintf('%s: using transformation matrix in avg_data structure...\n',...
        mfilename);
    end;
    if ~isempty(avg_data.coor_trans.mri2head)
      T_mri2head=avg_data.coor_trans.mri2head;
    end;
  end;
  if isempty(T_mri2head)
    error('no mri2head transformation was supplied');
  end;
  parms.trans = T_mri2head;
  if parms.verbose
    fprintf('%s: mri2head transformation matrix:\n',mfilename);
    disp(parms.trans);
  end;

  % determine which channels (sensors) to use
  parms = set_chans(parms,avg_data);

  % check signal-space-projection matrix
  if ~isempty(parms.ssp_projmat)
    pm_sz = size(parms.ssp_projmat);
    if length(pm_sz)~=2 | any(pm_sz~=[parms.nchans,parms.nchans])
      if parms.verbose
        fprintf('%s: ssp_projmat size:\n',mfilename);
        disp(pm_sz);
      end;
      error('ssp_projmat has wrong number of dimensions');
    end;
    parms.ssp_projmat = double(full(parms.ssp_projmat));
  end;

  % check forward solution parameters
  % check that nlayers is 1 or 3 if openmeeg ~= 1
  if ~parms.openmeeg_flag && ~ismember(parms.nlayers,[1,3])
    error('openmeeg_flag=1 required if number of layers not 1 or 3');    
  end;
  % check number of layers are consistent with MEG/EEG
  % EEG must have at least 3 layers
  if parms.nlayers<3 & parms.useEEG_flag
    error('must have at least 3 layers (BEM surfaces or spherical shells) for EEG');
  end;
  % bem vs. sphere
  if parms.bem_flag
    if iscell(parms.bem_surf_files)
      num_bem_surfs=length(parms.bem_surf_files);
    else
      parms.bem_surf_files={parms.bem_surf_files};
      num_bem_surfs=1;
    end;
    if num_bem_surfs > parms.nlayers % i.e. nlayers = 1
      parms.bem_surf_files = {parms.bem_surf_files{1:parms.nlayers}};
      num_bem_surfs=length(parms.bem_surf_files);
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

  % check dipole info files
  if ~parms.calc_dipinfo_flag
    if isempty(lh_dip_info)
      if mmil_isrelative(parms.lh_dip_file)
        parms.lh_dip_file = [parms.subjpath '/' parms.lh_dip_file];
      end;
      if ~exist(parms.lh_dip_file,'file')
        error('lh_dip_file %s not found',parms.lh_dip_file);
      end
    end;
    if isempty(rh_dip_info)
      if mmil_isrelative(parms.rh_dip_file)
        parms.rh_dip_file = [parms.subjpath '/' parms.rh_dip_file];
      end;
      if ~exist(parms.rh_dip_file,'file')
        error('rh_dip_file %s not found',parms.rh_dip_file);
      end
    end;
  end;  

  % check decimation files
  if ~parms.nodec_flag
    if isempty(parms.lh_dec_dips)
      if mmil_isrelative(parms.lh_dec_file)
        parms.lh_dec_file = [parms.subjpath '/' parms.lh_dec_file];
      end;
      if ~exist(parms.lh_dec_file,'file')
        error('lh_dec_file %s not found',parms.lh_dec_file);
      end;
    end
    if isempty(parms.rh_dec_dips)
      if mmil_isrelative(parms.rh_dec_file)
        parms.rh_dec_file = [parms.subjpath '/' parms.rh_dec_file];
      end;
      if ~exist(parms.rh_dec_file,'file')
        error('rh_dec_file %s not found',parms.rh_dec_file);
      end;
    end
  end;

  % check for smoothness constraint
  if parms.smooth_constr_flag
    if ~parms.nodec_flag
      if parms.verbose
        fprintf('%s: smooth_constr_flag=1, so forcing nodec_flag=1',...
          mfilename);
      end;
      parms.nodec_flag = 1;
    end;
    if ~parms.orient_constr_flag
      if parms.verbose
        fprintf('%s: smooth_constr_flag=1, so forcing orient_constr_flag=1',...
          mfilename);
      end;
      parms.orient_constr_flag = 1;
    end;
    % check ROI files
    if strcmp(parms.smooth_constr_type,'roi_nbrhd')
      [nroi,nhemi] = size(parms.smooth_constr_roi_files);
      if nhemi~=2
        error('smooth_constr_roi_files must be 2 (number of hemispheres)');
      end;
      parms.smooth_constr_nroi = nroi;
      % check that files exist
      for h=1:nhemi
        for r=1:nroi
          fname = parms.smooth_constr_roi_files{r,h};
          if ~exist(fname,'file')
            error('file %s not found',fname);
          end;
        end;
      end;
    end;
  end;

  % check template for writing fif files
  if parms.write_fif_flag && ~exist(parms.template_fif,'file')
    error('template_fif %s not found',parms.template_fif);
  end;

  % set output flags
  if parms.resamp2ico_flag
    parms.write_mgh_flag = 1;
  end;
  if parms.write_mgh_flag
    parms.write_stc_flag = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function create_output_dirs(parms)
  mmil_mkdir([parms.rootoutdir '/matfiles']);
  if parms.write_fif_flag
    mmil_mkdir([parms.rootoutdir '/fifs']);
  end;
  if parms.write_mgh_flag
    mmil_mkdir([parms.rootoutdir '/mghfiles']);
  end;
  if parms.write_stc_flag
    mmil_mkdir([parms.rootoutdir '/stcfiles']);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_dipinfo(parms)
  fname = sprintf('%s/matfiles/%s_dipinfo.mat',parms.rootoutdir,parms.prefix);
  if ~exist(fname,'file') || parms.forceflag
    % load dipole information
    if parms.calc_dipinfo_flag
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
      else
        lh_dip_info = parms.lh_dip_info;
      end;
      if isempty(parms.rh_dip_info)
        if parms.verbose
          fprintf('%s: loading right hemisphere dipole information from %s...\n',...
            mfilename,parms.rh_dip_file);
        end;
        rh_dip_info = ts_read_dip_file(parms.rh_dip_file);
      else
        rh_dip_info = parms.rh_dip_info;
      end;
    end;

    % load decimation information
    if parms.nodec_flag
      lh_dec_dips = ones(size(lh_dip_info,2),1);
      rh_dec_dips = ones(size(rh_dip_info,2),1);
    else
      if isempty(parms.lh_dec_dips)
        if parms.verbose
          fprintf('%s: loading lh dec file %s...\n',...
            mfilename,parms.lh_dec_file);
        end;
        lh_dec_dips = ts_read_dec_file(parms.lh_dec_file);
      else
        lh_dec_dips = parms.lh_dec_dips;
      end
      if isempty(parms.rh_dec_dips)
        if parms.verbose
          fprintf('%s: loading rh dec file %s...\n',...
            mfilename,parms.rh_dec_file);
        end;
        rh_dec_dips = ts_read_dec_file(parms.rh_dec_file);
      else
        rh_dec_dips = parms.rh_dec_dips;
      end
    end;
    % save dipole and decimation info
    if parms.verbose
      fprintf('%s: saving dipinfo to mat file %s...\n',mfilename,fname);
    end;
    save(fname,'lh_dip_info','rh_dip_info','lh_dec_dips','rh_dec_dips');
  else
    if parms.verbose
      fprintf('%s: loading dipinfo from mat file %s...\n',mfilename,fname);
    end;
    load(fname);
  end;
  parms.lh_dip_info = lh_dip_info;
  parms.rh_dip_info = rh_dip_info;
  parms.lh_dec_dips = lh_dec_dips;
  parms.rh_dec_dips = rh_dec_dips;

  % check for consistentcy
  parms.num_dips_lh=length(parms.lh_dip_info);
  parms.num_dips_rh=length(parms.rh_dip_info);
  parms.num_dec_dips_lh=length(find(parms.lh_dec_dips));
  parms.num_dec_dips_rh=length(find(parms.rh_dec_dips));
  if(parms.num_dips_lh==0 & parms.num_dips_rh==0)
    error('no dipole info found');
  end
  if(parms.num_dec_dips_lh==0 & parms.num_dec_dips_rh==0)
    error('no decimated dipoles specified');
  end
  if length(parms.lh_dec_dips) ~= parms.num_dips_lh
    error('length of lh_dec_dips must match number of lh dips');
  end
  if length(parms.rh_dec_dips) ~= parms.num_dips_rh
    error('length of rh_dec_dips must match number of rh dips');
  end
  if parms.verbose
    fprintf('%s: number of selected dipoles: lh=%d rh=%d\n',...
      mfilename,parms.num_dec_dips_lh,parms.num_dec_dips_rh);
  end;
  parms.ndips_lh = 3*parms.num_dec_dips_lh;
  parms.ndips_rh = 3*parms.num_dec_dips_rh;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % apply projection matrix to data
    if ~isempty(parms.ssp_projmat)
      if parms.verbose
        fprintf('%s: applying SSP projection matrix to data...\n',...
          mfilename);
      end;
      nconds = length(parms.conditions);
      for c=1:nconds
        cond = parms.conditions(c);
        avg_data.averages(cond).data = parms.ssp_projmat*avg_data.averages(cond).data;
      end;
    end;

    % baseline correction
    if parms.baseline_flag
      if parms.verbose
        fprintf('%s: subtracting baseline from average data...\n',...
          mfilename);
      end;
      avg_data = baseline_corr(parms,avg_data);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avg_data = baseline_corr(parms,avg_data)
  nconds = length(parms.conditions);
  for c=1:nconds
    cond = parms.conditions(c);
    time = avg_data.averages(cond).time;
    [tmp,baseline_start_samp] = min(abs(time - parms.baseline_start/1000));
    [tmp,baseline_end_samp] = min(abs(time - parms.baseline_end/1000));
    data = avg_data.averages(cond).data;
    baseline = data(:,baseline_start_samp:baseline_end_samp);
    mean_base=mean(baseline')';
    data = data - mean_base*ones(1,size(data,2));  % correct baseline
    avg_data.averages(cond).data = data;    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_parms(parms)
  fname = sprintf('%s/matfiles/%s_parms.mat',parms.rootoutdir,parms.prefix);
  if parms.verbose
    fprintf('%s: saving parameters to mat file %s...\n',mfilename,fname);
  end;
  save(fname,'parms');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_scalefacts(parms,avg_data)
  % calculate scaling factors for each channel (if necessary)
  if parms.calc_scalefacts_flag || parms.ncov_type==1
    parms = calc_channel_scalefacts(parms,avg_data);
  end;
  % calculate scaling factors for each channel type
  %   otherwise use default or user supplied scaling factors
  if parms.calc_scalefacts_flag && ~parms.prewhiten_flag
    parms = calc_chantype_scalefacts(parms);
  elseif parms.prewhiten_flag
    parms.grad_scalefact = 1;
    parms.mag_scalefact = 1;
    parms.EEG_scalefact = 1;
  end;
  % identify additional bad chans based on flat baseline
  parms = detect_flat_chans(parms,avg_data);
  if parms.verbose
    fprintf('%s: scaling factors: grad=%0.2g, mag=%0.2g, EEG=%0.2g\n',...
      mfilename,parms.grad_scalefact,parms.mag_scalefact,parms.EEG_scalefact);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_channel_scalefacts(parms,avg_data)
  if parms.verbose
    fprintf('%s: calculating scaling factors from data...\n',...
      mfilename);
  end;
  parms.scalefacts = zeros(parms.nchans,1);
  total_num_trials = 0;
  nconds = length(parms.ncov_conditions);
  for c=1:nconds
    cond = parms.ncov_conditions(c);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  % get scaling factors for each channel
  %   (either from raw ncov or from avg scale factors)
  C = create_noisecovar(parms,avg_data);
  scalefacts = diag(C);
  % multiply channel scale facts by channel type scale facts
  scale_vector = calc_scale_vector(parms);
  scalefacts = scale_vector .* scalefacts;
  scalefacts = scalefacts(parms.goodchans);
  thresh = mean(scalefacts)*eps;
  i_bad = find(scalefacts<thresh);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_noise(parms,avg_data)
  % get noise covariance matrix (either from raw or from avg scale factors)
  parms.noisecovar = create_noisecovar(parms,avg_data);

  % create noise covariance matrix for noise normalization
  if parms.noisenorm_flag && ~parms.prewhiten_flag
    if parms.noisenorm_identity_flag
      parms.noisecovar_noisenorm = ...
        create_identity_noisecovar(parms,parms.noisecovar);
    else
      parms.noisecovar_noisenorm = parms.noisecovar;
    end;
  else
    parms.noisecovar_noisenorm = [];
  end;

  % create matrix for scaling noise covariance matrix
  scale_matrix = calc_noise_scale_matrix(parms);

  % scale noise covariance matrix
  parms.noisecovar = parms.noisecovar.*scale_matrix;
  if ~isempty(parms.noisecovar_noisenorm)
    parms.noisecovar_noisenorm = parms.noisecovar_noisenorm.*scale_matrix;
  end;

  % exclude bad channels from noise covariance matrix
  parms.noisecovar = ...
    parms.noisecovar(parms.goodchans,parms.goodchans);
  if ~isempty(parms.noisecovar_noisenorm)
    parms.noisecovar_noisenorm = ...
      parms.noisecovar_noisenorm(parms.goodchans,parms.goodchans);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = create_noisecovar(parms,avg_data)
  % get noisecovar from avg_data.noise.covar  or create one from scalefacts
  if parms.ncov_type==0
    if parms.verbose
      fprintf('%s: using scaled identity matrix for noise covariance\n',...
        mfilename);
    end;
    C = create_identity_noisecovar(parms);
  elseif parms.ncov_type==1
    if parms.verbose
      fprintf('%s: using noise covariance calculated from average baseline\n',...
        mfilename);
    end;
    C = diag(parms.scalefacts.^2);
  elseif isempty(avg_data.noise.covar)
    if parms.verbose
      fprintf('%s: using scaled identity matrix for noise covariance\n',...
        mfilename);
    end;
    C = create_identity_noisecovar(parms);
  else
    if parms.verbose
      fprintf('%s: using noise covariance matrix from avg_data.noise.covar\n',...
        mfilename);
    end;
    C = avg_data.noise.covar;
    L = parms.raw_ncov_lambda; % regularization
    C = (1-L)*C + L*diag(diag(C));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = create_identity_noisecovar(parms,C)
  if ~exist('C','var'), C = []; end;
  if isempty(C)
    diagvals = zeros(parms.nchans,1);
    diagvals(parms.grad_chans) = 1.0/(parms.grad_scalefact^2);
    diagvals(parms.mag_chans) = 1.0/(parms.mag_scalefact^2);
    diagvals(parms.EEG_chans) = 1.0/(parms.EEG_scalefact^2);
  else
    diagvals = diag(C);
    diagvals(parms.grad_chans) = mean(diagvals(parms.grad_chans));
    diagvals(parms.mag_chans) = mean(diagvals(parms.mag_chans));
    diagvals(parms.EEG_chans) = mean(diagvals(parms.EEG_chans));
  end;
  C = diag(diagvals);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,G_xyz] = calc_forward(parms,avg_data)
  G_xyz = [];
  fname = sprintf('%s/matfiles/%s_forward_prep.mat',...
    parms.rootoutdir,parms.prefix);
  if ~exist(fname,'file') || parms.forceflag
    if ~isempty(parms.forward_matfile)
      if mmil_isrelative(parms.forward_matfile)
        forward_matfile = sprintf('%s/matfiles/%s',...
          parms.rootoutdir,parms.forward_matfile);
      else
        forward_matfile = parms.forward_matfile;
      end;
    else
      forward_matfile = sprintf('%s/matfiles/%s_forward.mat',...
        parms.rootoutdir,parms.prefix);
      parms.forward_matfile = forward_matfile;
    end;
    if ~exist(forward_matfile,'file') || parms.forceflag>1
      if parms.verbose
        fprintf('%s: calculating forward solution gain matrix...\n',mfilename);
      end;
      if parms.bem_flag
        parms.bem_matfile = sprintf('%s/matfiles/%s_bem.mat',...
          parms.rootoutdir,parms.prefix);
      end;
      args = mmil_parms2args(parms,parms.gainmat_tags);
      G_xyz = ts_calc_gainmat(avg_data,args{:});
      if isempty(G_xyz)
        error('failed to calculate gain matrix');
      end;
      save(forward_matfile,'-v7.3','G_xyz');
    else
      if parms.verbose
        fprintf('%s: loading pre-calculated forward solution from %s...\n',...
          mfilename,forward_matfile);
      end;
      load(forward_matfile);
    end;
    % prepare forward matrix: apply ssp, scale gain matrix,
    %                         apply orientation constraint, remove bad channels
    [G_xyz,num_sensors,num_sources] = prep_forward(parms,G_xyz);
    save(fname,'-v7.3','G_xyz','num_sensors','num_sources');
  else
    if parms.verbose
      fprintf('%s: loading pre-calculated prepared forward solution from %s...\n',...
        mfilename,fname);
    end;
    load(fname);
  end;
  parms.num_sensors = num_sensors;
  parms.num_sources = num_sources;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G_xyz,num_sensors,num_sources] = prep_forward(parms,G_xyz)
  % apply ssp projection matrix
  if ~isempty(parms.ssp_projmat)
    if parms.verbose
      fprintf('%s: applying SSP projection matrix to gain matrix...\n',...
        mfilename);
    end;
    G_xyz = parms.ssp_projmat*G_xyz;
  end;

  % scale gain matrix
  if parms.verbose
    fprintf('%s: scaling gain matrix...\n',mfilename);
  end;
  scale_matrix = calc_forward_scale_matrix(parms,G_xyz);
  G_xyz = G_xyz.*scale_matrix;

  if parms.orient_constr_flag
    % convert xyz to normal and two orthogonal tangential components
    if parms.verbose
      fprintf('%s: converting xyz components to normal and tangentials\n',...
        mfilename);
    end;
    [G_norm,G_tang1,G_tang2]=ts_gain_xyz2norm(G_xyz,parms.lh_dip_info,...
      parms.rh_dip_info,parms.lh_dec_dips,parms.rh_dec_dips,parms.trans);
    tmp = zeros(size(G_xyz));
    tmp(:,1:3:end) = G_norm;
    tmp(:,2:3:end) = G_tang1;
    tmp(:,3:3:end) = G_tang2;
    G_xyz = tmp;
  end;

  % resize to remove unwanted and bad chans
  if parms.verbose
    fprintf('%s: resizing gain matrix to remove bad channels...\n',mfilename);
  end;
  G_xyz = G_xyz(parms.goodchans,:);
  [num_sensors,num_sources]=size(G_xyz);
  if parms.verbose
    fprintf('%s: size(G_xyz)=[%d,%d]\n',mfilename,num_sensors,num_sources);
  end;
  num_sources = num_sources/3; % for 3 orientations
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = construct_srccovar(parms,G_xyz)
  % construct source covariance matrix
  fname = sprintf('%s/matfiles/%s_sourcecov.mat',...
    parms.rootoutdir,parms.prefix);
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: constructing source covariance matrix...\n',...
        mfilename);
    end;
    % source covariance for smoothness constraint
    if parms.smooth_constr_flag
      [Rsmooth_lh,Rsmooth_rh] = calc_sourcecov_smooth(parms);
    end;
    if ~isempty(parms.lh_sourcecov_file) || ~isempty(parms.rh_sourcecov_file)
      % create source covariance from mgh files (e.g. fMRI bias)
      [R_lh,R_rh] = load_sourcecov_files(parms);
      if parms.smooth_constr_flag
        [R_lh,R_rh] = combine_sourcecov_smooth(parms,...
                        R_lh,R_rh,Rsmooth_lh,Rsmooth_rh);
      end;
    else
      % source covariance as identity matrix
      if parms.smooth_constr_flag
        R_lh = Rsmooth_lh;
      else
        R_lh = speye(parms.ndips_lh,parms.ndips_lh);
      end;
      if parms.smooth_constr_flag
        R_rh = Rsmooth_rh;
      else
        R_rh = speye(parms.ndips_rh,parms.ndips_rh);
      end;
    end;
    R = combine_sourcecov_hemis(parms,R_lh,R_rh);
    % orient_tang for 2nd and 3rd components
    if parms.orient_constr_flag
      R = calc_sourcecov_orient_constr(parms,R);
    end;
    % depth weighting
    if parms.depthweight_flag
      R = calc_sourcecov_depthweight(parms,R,G_xyz);
    end;
    save(fname,'R');
  else
    if parms.verbose
      fprintf('%s: loading source covariance matrix from %s...\n',...
        mfilename,fname);
    end;
    load(fname);
  end;
  parms.srccovar = R;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rsmooth_lh,Rsmooth_rh] = calc_sourcecov_smooth(parms)
  surfname = 'white';
  switch parms.smooth_constr_type
    case 'smooth'
      fname = sprintf('%s/matfiles/%s_Rsmooth%d.mat',...
        parms.rootoutdir,parms.subjname,parms.smooth_constr_nsmooth);
    case 'nbrhd'
      fname = sprintf('%s/matfiles/%s_R%s.mat',...
        parms.rootoutdir,parms.subjname,parms.smooth_constr_nbrhd_outfix);
    case 'roi_nbrhd'
      fname = sprintf('%s/matfiles/%s_R%s.mat',...
        parms.rootoutdir,parms.subjname,parms.smooth_constr_roi_nbrhd_outfix);
  end;
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: creating %s source covariance matrix for smoothness constraint...\n',...
        mfilename,parms.smooth_constr_type);
    end;
    for h=1:length(parms.hemilist)
      hemi = parms.hemilist{h};
      surf = fs_load_subj(parms.subjname,hemi,[],[],parms.subjdir);
      switch parms.smooth_constr_type
        case 'smooth'
          R = ts_smoothconstr_covar(surf,parms.smooth_constr_nsmooth);
        case 'nbrhd'
          R = ts_nbrhd_covar(surf,...
            parms.smooth_constr_nbrhd_a0,...
            parms.smooth_constr_nbrhd_a1,...
            parms.smooth_constr_nbrhd_a2);
        case 'roi_nbrhd'
          nroi = parms.smooth_constr_nroi;
          roi_verts = cell(1,nroi);
          for r=1:nroi
            fname_roi = parms.smooth_constr_roi_files{r,h};
            roi_verts{r} = fs_read_label(fname_roi);
          end;
          R = ts_roi_nbrhd_covar(surf,roi_verts,...
            parms.smooth_constr_nbrhd_a0,...
            parms.smooth_constr_nbrhd_a1,...
            parms.smooth_constr_nbrhd_a2,...
            parms.smooth_constr_nbrhd_ax,...
            parms.smooth_constr_nbrhd_az);
      end;
      switch hemi
        case 'lh'
          Rsmooth_lh = R;
        case 'rh'
          Rsmooth_rh = R;
      end;
    end;
    save(fname,'Rsmooth_lh','Rsmooth_rh');
  else
    load(fname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_lh,R_rh] = load_sourcecov_files(parms)
  parms.sourcecov_thresh = max(parms.sourcecov_thresh,0);
  if ~isempty(parms.lh_sourcecov_file)
    if parms.verbose
      fprintf('%s: loading lh sourcecov mgh %s...\n',...
        mfilename,parms.lh_sourcecov_file);
    end;
    R_lh = ts_mgh2sourcecov(parms.lh_sourcecov_file,...
      'decdips',parms.lh_dec_dips,...
      'thresh_flag',parms.sourcecov_thresh_flag,...
      'thresh',parms.sourcecov_thresh,...
      'absflag',parms.sourcecov_thresh_abs_flag,...
      'maxvar',parms.sourcecov_maxvar,...
      'minvar',parms.sourcecov_minvar);
  else
    R_lh = speye(parms.ndips_lh,parms.ndips_lh)*parms.sourcecov_minvar;
  end
  if ~isempty(parms.rh_sourcecov_file)
    if parms.verbose
      fprintf('%s: loading rh sourcecov mgh %s...\n',...
        mfilename,parms.rh_sourcecov_file);
    end;
    R_rh = ts_mgh2sourcecov(parms.rh_sourcecov_file,...
      'decdips',parms.rh_dec_dips,...
      'thresh_flag',parms.sourcecov_thresh_flag,...
      'thresh',parms.sourcecov_thresh,...
      'absflag',parms.sourcecov_thresh_abs_flag,...
      'maxvar',parms.sourcecov_maxvar,...
      'minvar',parms.sourcecov_minvar);
  else
    R_rh = speye(parms.ndips_rh,parms.ndips_rh)*parms.sourcecov_minvar;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_lh,R_rh] = combine_sourcecov_smooth(parms,R_lh,R_rh,Rsmooth_lh,Rsmooth_rh)
  if parms.verbose
    fprintf('%s: combining source covariance with smoothness constraint...\n',...
      mfilename);
  end;
  %% todo: requires testing
  tmp_R_lh = full(diag(R_lh));
  for i=1:length(tmp_R_lh)
    Rsmooth_lh(i,:) = Rsmooth_lh(i,:)*tmp_R_lh(i);
  end;
  tmp_R_rh = full(diag(R_rh));
  for i=1:length(tmp_R_rh)
    Rsmooth_rh(i,:) = Rsmooth_rh(i,:)*tmp_R_rh(i);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = combine_sourcecov_hemis(parms,R_lh,R_rh)
  if parms.verbose
    fprintf('%s: combining left and right source covariance matrices...\n',...
      mfilename);
  end;
  % combine R_lh and R_rh
  if parms.smooth_constr_flag
    ndips = parms.ndips_lh + parms.ndips_rh;
    R=sparse(ndips,ndips); %% todo: may be faster to use non-sparse matrix
    R(1:parms.ndips_lh,1:parms.ndips_lh) = R_lh;
    R(parms.ndips_lh+1:end,parms.ndips_lh+1:end) = R_rh;
    % takes about 17 seconds
  else  % diagonal only
    ndips = parms.ndips_lh + parms.ndips_rh;
    R=sparse(ndips,ndips);
    for i=1:parms.ndips_lh
      R(i,i)=R_lh(i,i);
    end
    j=1;
    for i=parms.ndips_lh+1:ndips
      R(i,i)=R_rh(j,j);
      j=j+1;
    end
    % takes < 1 second
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = calc_sourcecov_orient_constr(parms,R)
  if parms.verbose
    fprintf('%s: setting diagonal tangential components of source covariance matrix...\n',...
      mfilename);
  end;
  ndips = parms.ndips_lh + parms.ndips_rh;
  n = ndips/3;
  j = 1;
  for i=1:n
    for k=j+1:j+2
      R(k,k) = R(j,j)*parms.orient_tang;
    end;
    j = j+3;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = calc_sourcecov_depthweight(parms,R,G_xyz)
  if parms.verbose
    fprintf('%s: calculating depth weighted source covariance...\n',...
      mfilename);
  end;
  R_dw = ts_depthweight_covar(G_xyz,parms.depthweight_p);
  if parms.smooth_constr_flag
    if parms.verbose
      fprintf('%s: combining smoothness constraint and depthweighting...\n',...
        mfilename);
    end;
    %% todo: requires testing
    tmp_R_dw = full(diag(R_dw));
    for i=1:length(tmp_R_dw)
      R(i,:) = R(i,:)*tmp_R_dw(i);
    end;
  else
    R = R.*R_dw;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,nnf] = calc_inverse(parms,G_xyz)
  if ~isempty(parms.inverse_matfile)
    fname = parms.inverse_matfile;
  else
    fname = sprintf('%s/matfiles/%s_inverse.mat',...
      parms.rootoutdir,parms.prefix);
  end;
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating inverse operator...\n',mfilename);
    end;
    [M,nnf] = ts_calc_inverse_operator(G_xyz,...
      'prewhiten_flag',parms.prewhiten_flag,...
      'SNR',parms.SNR,'C',parms.noisecovar,...
      'R',parms.srccovar,'C_nn',parms.noisecovar_noisenorm);
    if isempty(M)
      error('failed to calculate inverse operator');
    end;
    save(fname, 'M', 'nnf');
  else
    if parms.verbose
      fprintf('%s: loading pre-calculated inverse operator from %s...\n',...
        mfilename,fname);
    end;
    load(fname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_source_estimates(parms,avg_data,M,nnf,G_xyz)
  % initialize avg_data structure for fit and err
  if parms.save_fiterr_struct_flag
    % initialize avg_data structure for least-squares fit
    if parms.verbose
      fprintf('%s: initializing best fit structure...\n',mfilename);
    end;
    fit = init_fit_struct(parms,avg_data);
    % initialize avg_data structure for residual error
    if parms.verbose
      fprintf('%s: initializing residual error structure...\n',mfilename);
    end;
    err = init_fit_struct(parms,avg_data);
  else
    fit = []; err = [];
  end;

  % calculate scale matrix
  scale_matrix = calc_data_scale_matrix(parms,avg_data);

  % apply inverse to data for each condition
  for c=1:length(parms.conditions)
    cond = parms.conditions(c);
    % check for potential output files
    if ~parms.forceflag
      allexist = check_output(parms,cond);
      if allexist, continue; end;        
    end;
    [S_xyz,S,F,Y] = apply_inverse(parms,cond,avg_data,scale_matrix,M,nnf);
    if parms.save_results_flag
      [fit,err] = save_results(parms,cond,S_xyz,S,F,Y,G_xyz,fit,err);
    end;
    if parms.write_stc_flag
      write_stc_files(parms,cond,S,F,avg_data);
      if parms.write_mgh_flag
        write_mgh_files(parms,cond);
      end;
    end;
  end;

  % save fit and err sensor waveforms
  if parms.save_fiterr_struct_flag
    fname = sprintf('%s/matfiles/%s_fiterr.mat',parms.rootoutdir,parms.prefix);
    if ~exist(fname) || parms.forceflag
      save(fname,'fit', 'err');
    end;
    if parms.write_fif_flag
      % save fit
      outstem = sprintf('%s/fifs/%s_fit',parms.rootoutdir,parms.prefix);
      if parms.verbose
        fprintf('%s: writing fitted average to fif...\n',mfilename);
      end;
      ts_avg2fif(fit,parms.template_fif,outstem,[],parms.forceflag);
      % save err
      outstem = sprintf('%s/fifs/%s_err',parms.rootoutdir,parms.prefix);
      if parms.verbose
        fprintf('%s: writing residual error to fif...\n',mfilename);
      end;
      ts_avg2fif(err,parms.template_fif,outstem,[],parms.forceflag);
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit = init_fit_struct(parms,avg_data)
  fit = [];
  fit.num_sensors = avg_data.num_sensors;
  fit.sfreq = avg_data.sfreq;
  fit.sensor_info = avg_data.sensor_info;
  fit.coor_trans = avg_data.coor_trans;
  fit.noise = avg_data.noise;
  for c=1:length(avg_data.averages)
    fit.averages(c).event_code=avg_data.averages(c).event_code;
    fit.averages(c).num_trials=avg_data.averages(c).num_trials;
    fit.averages(c).num_rejects=avg_data.averages(c).num_rejects;
    fit.averages(c).time=avg_data.averages(c).time;
    fit.averages(c).data=zeros(size(avg_data.averages(c).data));
    fit.averages(c).stdev=zeros(size(avg_data.averages(c).data));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S_xyz,S,F,Y] = apply_inverse(parms,cond,avg_data,scale_matrix,M,nnf)
  % initialize variables
  num_tpoints = size(avg_data.averages(1).data,2);
  S = zeros(num_tpoints,parms.num_sources,parms.datatype);
  F = zeros(num_tpoints,parms.num_sources,parms.datatype);
  % apply inverse to data
  if parms.verbose
    fprintf('%s: applying inverse for cond %d...\n',mfilename,cond);
  end;
  % first apply scale factors to data
  data = avg_data.averages(cond).data.*scale_matrix;
  Y = data(parms.goodchans,:)';
  S_xyz = (M*Y')';
  % calculate amplitudes and do noise normalization
  if parms.verbose
    fprintf('%s: calculating source amplitudes...\n',...
      mfilename);
  end;
  % L should be number of trials in average, unless noise covariance
  %   calculated from avg data, in which case it should be 1
  if ismember(parms.ncov_type,[0,1]) || isempty(parms.noisecovar)
    L = 1;
  else
    L = avg_data.averages(cond).num_trials;
  end;
  nnf2L = (nnf.^2)/L;
  nnfL = nnf/L;
  for s=1:parms.num_sources
    j = 3*(s-1) + 1;
    k = j + 2;
    if parms.signed_sources_flag && parms.orient_constr_flag
      S(:,s) = S_xyz(:,j);
      if parms.noisenorm_flag
        % noise normalization
        F(:,s) = S(:,s)/nnfL(j);
      end;
    else
      source = S_xyz(:,j:k);
      S(:,s) = sqrt(sum(source.^2,2));
      if parms.noisenorm_flag
        % noise normalization
        F(:,s) = (S(:,s).^2)/sum(nnf2L(j:k));
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function allexist = check_output(parms,cond)
  allexist = 1;
  output_files = []; j=1;
  output_files{j} = sprintf('%s/matfiles/%s_results_cond%02d.mat',...
    parms.rootoutdir,parms.prefix,cond);
  j=j+1;
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    if parms.write_stc_flag
      output_files{j} = sprintf('%s/stcfiles/%s_cond%02d-%s.stc',...
        parms.rootoutdir,parms.prefix,cond,hemi);
      j=j+1;
    end;
    if parms.write_mgh_flag
      output_files{j} = sprintf('%s/mghfiles/%s_cond%02d-spsm%d-sm%d-%s.mgh',...
        parms.rootoutdir,parms.prefix,cond,...
        parms.sparsesmooth,parms.postsmooth,...
        hemi);
    end;
  end;
  for j=1:length(output_files)
    if ~exist(output_files{j},'file')
      allexist = 0;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fit,err] = save_results(parms,cond,S_xyz,S,F,Y,G_xyz,fit,err)
  if parms.save_fiterr_struct_flag
    % calculate fit and error
    if parms.verbose
      fprintf('%s: calculating fit and residual error...\n',mfilename);
    end;
    % calculate fit and residual error
    Yfit = (G_xyz*S_xyz')';
    E = Y-Yfit;
    % add to fitted average and error structures
    fit.averages(cond).data(parms.goodchans,:)=Yfit';
    err.averages(cond).data(parms.goodchans,:)=E';
  end;

  fname = sprintf('%s/matfiles/%s_results_cond%02d.mat',...
    parms.rootoutdir,parms.prefix,cond);
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: saving results to mat file %s...\n',mfilename,fname);
    end;
    if parms.save_fiterr_flag
      if ~parms.save_fiterr_struct_flag
        % calculate fit and residual error
        Yfit = (G_xyz*S_xyz')';
        E = Y-Yfit;
      end;
      % calculate separate error for each channel type
      norm_var_E = 0;
      ntypes = 0;
      if parms.useEEG_flag
        [good_chans,ind,ind_good_eeg] = ...
          intersect(parms.EEG_chans,parms.goodchans);
        Yeeg = Y(:,ind_good_eeg);
        Yfiteeg = Yfit(:,ind_good_eeg);
        Eeeg = Yeeg - Yfiteeg;
        var_Eeeg = var(Eeeg,0,2);
        var_Yeeg = var(Yeeg,0,2);
        var_Yfiteeg = var(Yfiteeg,0,2);
        tmp_var_Yeeg = var_Yeeg;
        tmp_var_Yeeg(find(~var_Yeeg))=1;
        norm_var_Eeeg = var_Eeeg./var_Yeeg;
        norm_var_E = norm_var_E + norm_var_Eeeg;
        ntypes = ntypes + 1;
      else
        Yeeg = [];
        Yfiteeg = [];
        Eeeg = [];
        var_Eeeg = [];
        var_Yeeg = [];
        norm_var_Eeeg = [];
      end;
      if parms.usegrad_flag
        [good_chans,ind,ind_good_grad] =...
          intersect(parms.grad_chans,parms.goodchans);
        Ygrad = Y(:,ind_good_grad);
        Yfitgrad = Yfit(:,ind_good_grad);
        Egrad = Ygrad - Yfitgrad;
        var_Egrad = var(Egrad,0,2);
        var_Ygrad = var(Ygrad,0,2);
        tmp_var_Ygrad = var_Ygrad;
        tmp_var_Ygrad(find(~var_Ygrad))=1;
        norm_var_Egrad = var_Egrad./var_Ygrad;
        norm_var_E = norm_var_E + norm_var_Egrad;
        ntypes = ntypes + 1;
      else
        Ygrad = [];
        Yfitgrad = [];
        Egrad = [];
        var_Egrad = [];
        var_Ygrad = [];
        norm_var_Egrad = [];
      end;
      if parms.usemag_flag
        [good_chans,ind,ind_good_mag] = ...
          intersect(parms.mag_chans,parms.goodchans);
        Ymag = Y(:,ind_good_mag);
        Yfitmag = Yfit(:,ind_good_mag);
        Emag = Ymag - Yfitmag;
        var_Emag = var(Emag,0,2);
        var_Ymag = var(Ymag,0,2);
        tmp_var_Ymag = var_Ymag;
        tmp_var_Ymag(find(~var_Ymag))=1;
        norm_var_Emag = var_Emag./var_Ymag;
        norm_var_E = norm_var_E + norm_var_Emag;
        ntypes = ntypes + 1;
      else
        Ymag = [];
        Yfitmag = [];
        Emag = [];
        var_Emag = [];
        var_Ymag = [];
        norm_var_Emag = [];
      end;
      norm_var_E = norm_var_E/ntypes;
      save(fname,'F','S','S_xyz',...
        'Y', 'Yfit','E','norm_var_E',...
        'Yeeg','Yfiteeg','Eeeg','var_Eeeg','var_Yeeg','norm_var_Eeeg',...
        'Ygrad','Yfitgrad','Egrad','var_Egrad','var_Ygrad','norm_var_Egrad',...
        'Ymag','Yfitmag','Emag','var_Emag','var_Ymag','norm_var_Emag');
    else
      save(fname,'F','S','S_xyz');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stc_files(parms,cond,S,F,avg_data)
  % write to lh stc file
  fname = sprintf('%s/stcfiles/%s_cond%02d-lh.stc',parms.rootoutdir,parms.prefix,cond);
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: saving to stc file %s...\n',mfilename,fname);
    end;
    lh_vnums=find(parms.lh_dec_dips==1)-1;
    num_lh_verts = length(lh_vnums);
    time = avg_data.averages(cond).time(1)*1000;
    % write_stc wants sfreq in Hz, t0 in ms
    if parms.signed_sources_flag && parms.orient_constr_flag &&...
       parms.noisenorm_flag
      tmp = F(:,1:num_lh_verts)'*parms.stc_scalefact;
    elseif parms.noisenorm_flag
      tmp = sqrt(F(:,1:num_lh_verts))'*parms.stc_scalefact;
    else
      tmp = S(:,1:num_lh_verts)'*parms.stc_scalefact;
    end;
    ts_write_stc(fname,lh_vnums,avg_data.sfreq,...
      time,tmp);
  end;

  % write to rh stc file
  fname = sprintf('%s/stcfiles/%s_cond%02d-rh.stc',parms.rootoutdir,parms.prefix,cond);
  if ~exist(fname,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: saving to stc file %s...\n',mfilename,fname);
    end;
    lh_vnums=find(parms.lh_dec_dips==1)-1;
    num_lh_verts = length(lh_vnums);
    rh_vnums=find(parms.rh_dec_dips==1)-1;
    time = avg_data.averages(cond).time(1)*1000;
    % write_stc wants sfreq in Hz, t0 in ms
    if parms.signed_sources_flag && parms.orient_constr_flag &&...
       parms.noisenorm_flag
      tmp = F(:,num_lh_verts+1:end)'*parms.stc_scalefact;
    elseif parms.noisenorm_flag
      tmp = sqrt(F(:,num_lh_verts+1:end))'*parms.stc_scalefact;
    else
      tmp = S(:,num_lh_verts+1:end)'*parms.stc_scalefact;
    end;
    ts_write_stc(fname,rh_vnums,avg_data.sfreq,...
      time,tmp);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_mgh_files(parms,cond)
  hemilist = {'lh' 'rh'};
  for h=1:length(hemilist)
    hemi = hemilist{h};
    stcfile = sprintf('%s/stcfiles/%s_cond%02d-%s.stc',...
      parms.rootoutdir,parms.prefix,cond,hemi);
    mghfile = sprintf('%s/mghfiles/%s_cond%02d-spsm%d-sm%d-%s.mgh',...
      parms.rootoutdir,parms.prefix,cond,...
      parms.sparsesmooth,parms.postsmooth,...
      hemi);
    if ~exist(mghfile,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: converting %s stc file to mgh file...\n',mfilename,hemi);
      end;
      ts_stc2mgh(stcfile,mghfile,parms.subjname,hemi,...
        parms.sparsesmooth,parms.postsmooth,...
        [],[],parms.subjdir,parms.mbmask_flag,parms.forceflag);
    end;
    if parms.resamp2ico_flag
      infile = mghfile;
      outfile = fs_surf2ico(infile,parms.subjname,...
        'icolevel',parms.icolevel,'smooth_out',parms.icosmooth,...
        'subjdir',parms.subjdir,'forceflag',parms.forceflag);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: flag to specify whether avg_data.noise.covar
%%       is from average or from raw
%% todo: option to pass in noise covar as parameter
%% todo: option to read data from fif file (need new function ts_fif2avg)
%% todo: extra dips (e.g. subcortical)
%% todo: optionally resample stc's with spline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

