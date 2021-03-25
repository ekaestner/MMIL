function dti_transform_fibers(fiber_dir,varargin)
%function dti_transform_fibers(fiber_dir,[options])
%
% Purpose: modify existing fiber ROIs by resampling to T1,
%   excluding CSF and gray matter, excluding voxels with multiple fibers,
%   or creating a fiber segmentation volume
%
% Required Input:
%   fiber_dir: full path of directory containing fiber ROI files
%
% Optional Parameters:
%   'outdir': output directory
%     {default = pwd}
%   'fiber_outdir': output directory for resampled or masked fibers
%     full path or relative to outdir
%     {default = 'fibers'}
%   'fibers': fiber numbers to resample or mask
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so, what type
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'resT1flag': [0|1] resample fibers to T1 resolution
%     {default = 0}
%   'M_T1_to_DTI': registration matrix between DTI data and high-res T1 data
%     {default = eye(4)}
%   'volsz_T1': dimensions of high-res T1 volume
%     {default = [256 256 256]}
%   'M_T1': vox2ras matrix for high-res T1 volume
%     {default = [-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1]}
%   'save_mgz_flag': [0|1] save fibers in mgz format in addition to sparse
%     {default = 0}
%   'verbose': display status messages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Optional Parameters for aseg masking:
%   'xcg_flag': [0|1] generate fiber ROIs excluding CSF and gray matter
%     as defined by FreeSurfer aseg
%     {default = 0}
%   'fname_aseg': full path name of aseg file
%     expected to have T1 resolution
%     required only if xcg_flag=1
%     {default = []}
%   'xcg_suffix':  suffix attached to output fiber file names
%       after excluding CSF and gray matter
%     {default = 'xcg'}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     Implicitly sets fseg_flag=1
%     {default = 0}
%   'masksf_suffix': suffix attached to output fiber file names
%       after excluding voxels with multiple fibers
%     {default = 'masksf'}
%
% Optional Parameters for dispersion weighting:
%   'disp_flag': [0|1] generate fiber ROIs excluding CSF and gray matter
%     as defined by FreeSurfer aseg
%     {default = 0}
%   'fname_vals': full path of value file (with dimensions matching fibers)
%     required for dispersion weighting
%     {default = []}
%   'disp_suffix':  suffix attached to output fiber file names
%       after calculating dispersion weighting
%     {default = 'dwtd'}
%   'disp_scalefact': scaling factor applied to values in fname_vals
%     {default = 1}
%   'dispvec': vector of dispersion values (MAD estimates)
%     may be a single value or vector with size matching number of fibers
%     {default = 0.1}
%   'dispfact': multiplicative factor applied to dispersion values
%     {default = 4}
%   'disp_xcg_flag': [0|1] interaction between disp weighting and xcg
%     0: calculate dispersion weighting everywhere, mask out exclude xcg voxels
%     1: calculate dispersion weighting only in xcg voxels
%     {default = 0}
%
% Optional Parameters for thresholding:
%   'thresh_FA': fractional anisotropy threshold applied to fiber ROIs
%     {default = 0}
%   'fname_FA': full path name of FA (fractional anisotropy) image file
%     Required if thresh_FA>0
%     {default = []}
%   'thresh_prob': fiber probability threshold applied to fiber ROIs
%     {default = 0}
%
% Optional Parameters for fiber segmentation volume
%   'fseg_flag': [0|1] generate fiber segmentation volume
%     {default = 1}
%   'fname_fseg': file name of output segmentation volume
%     unless absolute path is given, is relative to outdir
%     {default = 'fseg.mgz'}
%   'fseg_fibers': fiber numbers to include in fiber segmentation
%     by default, excludes Fmaj (121), Fmin (122), SLF (133,134), SCS (141,142)
%     {default = [103:110,115:120,123,135:138,1011,1021]}
%   'fseg_resT1flag': [0|1] gemerate fiber segmentation in T1 resolution
%     but do not save individual resampled fibers (if resT1flag = 0)
%     ignored if resT1flag = 1
%     {default = 0}
%   'fseg_xcg_flag': [0|1] exclude CSF and gray matter from fseg
%     but do not save individual masked fibers (if xcg_flag = 0)
%     ignored if resT1flag = 1
%     {default = 0}
%   'fseg_thresh_prob': fiber probability threshold applied to fseg
%     ignored if thresh_prob > 0
%     {default = 0.08}
%
% Optional Parameters for streamline path generation
%   'create_paths_flag': [0|1] generate DTI Studio format fiber paths
%     ignored if resT1flag = 1
%     {default = 1}
%   'paths_outdir': output directory for fiber paths
%     full path or relative to outdir
%     {default = 'paths'}
%   'fname_V0': full path name of V0 (primary Eigenvector) image file
%     Required if create_paths_flag = 1
%     {default = []}
%   'min_fiberlen': minimum fiber length for fiber path generation
%     {default = 12}
%   'thresh_angle': maximum turning angle for fiber path generation
%     {default = 70}
%   'path_suffix': string attached to fiber paths
%     if empty, will attach a long string including information about parameters
%       e.g. pthresh0.02_threshFA0.00_minflen10_loc_countatlas_path
%     {default = []}
%
% Created:  10/11/12 by Don Hagler
% Last Mod: 02/25/13 Don Hagler
%

%% todo: create paths even if resT1flag = 1
%%       need to resample FA and V0 to T1 resolution
%%        but also problematic because fact_suabe assumes axial slices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
parms = check_input(fiber_dir,varargin);

if parms.resT1flag
  % resample fibers to T1 resolution
  resamp_fibers(parms);
  parms.fiber_dir = parms.fiber_outdir;
  %% todo: if create_paths_flag, resample FA and V0 to T1 resolution
end;

% use aseg to exclude CSF and gray matter or
%   calculate dispersion weighting or
%   apply probability or FA thresholds
if parms.xcg_flag || parms.disp_flag ||...
   parms.thresh_prob>0 || parms.thresh_FA>0
  mask_fibers(parms);
end;

% create fiber segmentation volume
if parms.fseg_flag
  create_fseg(parms);
end;

% use fiber segmentation to exclude voxels with multiple fibers
if parms.masksf_flag
  mask_fibers(parms,1);
end;

% generate DTI Studio format fiber paths
if parms.create_paths_flag
  create_paths(parms);
end;
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fiber_dir,args)
  parms_filter = {...
    'fiber_dir',fiber_dir,[],...
...
    'outdir',pwd,[],...
    'fiber_outdir','fiber_maps',[],...
    'atlas_flag',2,[0:4],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'resT1flag',false,[false true],...
    'M_T1_to_DTI',eye(4),[],...
    'volsz_T1',[256,256,256],[],...
    'M_T1',[-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1],[],...
    'save_mgz_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
...
    'xcg_flag',false,[false true],...
    'fname_aseg',[],[],...
    'xcg_suffix','xcg',[],...
    'masksf_flag',false,[false true],...
    'masksf_suffix','masksf',[],...
...
    'disp_flag',false,[false true],...
    'fname_vals',[],[],...
    'disp_suffix','dwtd',[],...
    'dispvec',0.1,[],...
    'disp_scalefact',1,[1e-10,1e10],...
    'dispfact',4,[1e-6,1e6],...
    'disp_xcg_flag',false,[false true],...
...
    'thresh_FA',0,[0 1],...
    'fname_FA',[],[],...
    'thresh_prob',0,[0 1],...
...
    'fseg_flag',true,[false true],...
    'fseg_fibers',[101:110,115:120,123,135:138,143:150],[],...
    'fseg_resT1flag',false,[false true],...
    'fseg_xcg_flag',false,[false true],...
    'fseg_thresh_prob',0.08,[0 1],...
    'fname_fseg','fseg.mgz',[],...
...
    'create_paths_flag',true,[false true],...
    'fname_V0',[],[],...
    'min_fiberlen',12,[1 1000],...
    'thresh_angle',70,[30 90],...
    'path_suffix',[],[],...
    'paths_outdir','fiber_paths',[],...
... % hidden parameters
    'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
    'roicode_base',10000,[0,Inf],... % value added to fiber codes for fseg
    'count_flag',true,[false true],... % only if atlas_flag = 0
    'prob_exponent',1,[],...
    'orient_ref','LPS',[],...
... % parameters to pass to various functions
    'resT1_tags',{'outdir','fibers','M_T1_to_DTI','volsz_T1',...
                  'M_T1','atlas_flag','save_mgz_flag','verbose','forceflag',...
                  'fiber_outext','count_flag','interpm','bclamp'},[],...
    'mask_tags',{'outdir','fibers','resT1flag','atlas_flag','mask_suffix',...
                 'save_mgz_flag','verbose','forceflag',...
                 'fname_aseg','M_T1_to_DTI','aseg_codes','exclude_flag',...
                 'fname_vals','dispvec','disp_scalefact',...
                 'dispfact','disp_mask_flag',...
                 'thresh_FA','fname_FA','thresh_prob',...
                 'count_flag'},[],...
    'fseg_tags',{'fname_fseg','outdir','fibers','resT1flag',...
                 'M_T1_to_DTI','volsz_T1','M_T1',...
                 'atlas_flag','thresh_prob','thresh_FA',...
                 'fname_FA','xcg_flag','fname_aseg','verbose','forceflag',...
                 'count_flag','roicode_multiple','roicode_base','interpm',...
                 'bclamp','xcg_codes'},[],...
    'infix_tags',{'atlas_flag','count_flag','resT1flag',...
                  'xcg_flag','xcg_suffix',...
                  'masksf_flag','masksf_suffix',...
                  'disp_flag','disp_suffix','dispfact',...
                  'thresh_prob','thresh_FA'},[],...
    'suffix_tags',{'xcg_flag','xcg_suffix',...
                   'masksf_flag','masksf_suffix',...
                   'disp_flag','disp_suffix','dispfact',...
                   'thresh_prob','thresh_FA'},[],...
    'track_tags',{'rois','fname_prob','prob_exponent','outstem','roidir',...
                  'thresh_prob','thresh_FA','thresh_angle','min_fiberlen',...
                  'flipx_flag','flipy_flag','flipz_flag','verbose',...
                  'imgorient','imgseq','volsz','voxsz'},[],...
  };
  parms = mmil_args2parms(args,parms_filter);
  % check if input exists
  if ~exist(parms.fiber_dir,'dir')
    error('fiber_dir %s not found',parms.fiber_dir);
  end;
  if parms.xcg_flag
    if isempty(parms.fname_aseg)
      error('must specify fname_aseg if xcg_flag=1');
    end;
    if ~exist(parms.fname_aseg,'file')
      error('file %s not found',parms.fname_aseg);
    end;
  end;
  if parms.disp_flag
    if isempty(parms.fname_vals)
      error('must specify fname_vals if disp_flag=1');
    end;
    if ~exist(parms.fname_vals,'file')
      error('file %s not found',parms.fname_vals);
    end;
  end;
  if parms.thresh_FA>0
    if isempty(parms.fname_FA)
      error('must specify fname_FA if thresh_FA>0');
    end;
    if ~exist(parms.fname_FA,'file')
      error('file %s not found',parms.fname_FA);
    end;
  end;
  if parms.create_paths_flag
    if parms.resT1flag
      fprintf('%s: WARNING: creating paths not currently supported for resT1flag = 1\n',...
        mfilename);
      parms.create_paths_flag = 0;
    end;
    if isempty(parms.fname_FA) || isempty(parms.fname_V0)
      error('must specify fname_FA and fname_V0 if create_paths_flag = 1');
    end;
    if ~exist(parms.fname_FA,'file')
      error('file %s not found',parms.fname_FA);
    end;
    if ~exist(parms.fname_V0,'file')
      error('file %s not found',parms.fname_V0);
    end;
  end;
  if mmil_isrelative(parms.fiber_outdir)
    parms.fiber_outdir = [parms.outdir '/' parms.fiber_outdir];
  end;
  if mmil_isrelative(parms.paths_outdir)
    parms.paths_outdir = [parms.outdir '/' parms.paths_outdir];
  end;
  if mmil_isrelative(parms.fname_fseg)
    parms.fname_fseg = [parms.outdir '/' parms.fname_fseg];
  end;
  switch parms.atlas_flag
    case 0
      parms.locflag = 0;
      parms.countflag = 0;
    case 1
      parms.locflag = 1;
      parms.countflag = 1;
    case 2
      parms.locflag = 0;
      parms.countflag = 1;
    case 3
      parms.locflag = 1;
      parms.countflag = 0;
    case 4
      parms.locflag = 0;
      parms.countflag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_fibers(parms)
  orig_parms = parms;
  parms.outdir = parms.fiber_outdir;
  args = mmil_parms2args(parms,parms.resT1_tags);
  dti_resample_fibers(parms.fiber_dir,args{:})
  parms = orig_parms;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask_fibers(parms,masksf_flag)
  if ~exist('masksf_flag','var') || isempty(masksf_flag)
    masksf_flag = 0;
  end;
  parms.outdir = parms.fiber_outdir;
  parms.masksf_flag = masksf_flag;
  if parms.masksf_flag
    parms.fname_aseg = parms.fname_fseg;
    parms.aseg_codes = parms.fseg_fibers + parms.roicode_base;
    parms.exclude_flag = 0;
    parms.disp_flag = 0;
  elseif parms.xcg_flag
    parms.aseg_codes = parms.xcg_codes;
    parms.exclude_flag = 1;
  else
    parms.fname_aseg = [];
  end;
  if ~parms.disp_flag, parms.fname_vals = []; end;
  if parms.disp_flag && parms.xcg_flag
    parms.disp_mask_flag = parms.disp_xcg_flag;
  end;
  args = mmil_parms2args(parms,parms.suffix_tags);
  parms.mask_suffix = dti_set_fiber_suffix(args{:});
  args = mmil_parms2args(parms,parms.mask_tags);
  dti_mask_fibers(parms.fiber_dir,args{:})
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_fseg(parms)
  parms.fibers = parms.fseg_fibers;
  if parms.thresh_prob <= 0
    parms.thresh_prob = parms.fseg_thresh_prob;
  end;
  if parms.fseg_resT1flag, parms.resT1flag = 1; end;
  if parms.fseg_xcg_flag, parms.xcg_flag = 1; end;
  args = mmil_parms2args(parms,parms.fseg_tags);
  dti_fiber_segmentation(parms.fiber_dir,args{:})
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_paths(parms)
  mmil_mkdir(parms.paths_outdir);
  % set fiber_infix and fiber_ext
  args = mmil_parms2args(parms,parms.infix_tags);
  [fiber_infix,fiber_ext] = dti_set_fiber_infix(args{:});
  for f=parms.fibers
    parms.fname_prob = sprintf('%s/fiber_%02d_%s%s',...
      parms.fiber_dir,f,fiber_infix,fiber_ext);
    if isempty(parms.path_suffix)
      parms.outstem = sprintf('%s/fiber_%02d_%s_minflen%d',...
        parms.paths_outdir,f,fiber_infix,parms.min_fiberlen);
      parms.outstem = [parms.outstem '_path'];
    else
      parms.outstem = sprintf('%s/fiber_%02d_%s',...
        parms.paths_outdir,f,parms.path_suffix);
    end;
    fname_out = sprintf('%s.grp',parms.outstem);
    if ~exist(fname_out,'file') || parms.forceflag
      if ~exist(parms.fname_prob,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,parms.fname_prob);
      else
        orient = fs_read_orient(parms.fname_FA);
        [permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
        flipflags = flipvec<0;
        parms.flipx_flag = flipflags(1);
        parms.flipy_flag = flipflags(2);
        parms.flipz_flag = flipflags(3);
        args = mmil_parms2args(parms,parms.track_tags);
        dti_track_fibers(parms.fname_FA,parms.fname_V0,args{:});
      end;
    end;
  end;
return;

%% todo: make fact_suabe work with different orientations than axial
%%       or reorient volumes to axial,
%%         but need to reorient streamlines too
