function dti_fiber_segmentation(fiber_dir,varargin)
%function dti_fiber_segmentation(fiber_dir,varargin)
%
% Required Input:
%   fiber_dir: full path of directory containing fiber ROI files
%
% Optional Parameters:
%   'fname_fseg': file name of output segmentation volume
%     unless absolute path is given, is relative to outdir
%     {default = 'fseg.mgz'}
%   'outdir': output directory
%     {default = pwd}
%   'fibers': fiber numbers to convert
%     by default, excludes Fmaj (121), Fmin (122), SLF (133,134), SCS (141,142)
%     {default = [103:110,115:120,123,135:138,1011,1021]}
%   'atlas_flag': whether to use atlas fibers and if so, what type
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'resT1flag': [0|1] whether to use fibers resampled to T1 resolution
%     if resampled fibers not found, will resample fibers using M_T1_to_DTI
%     {default = 0}
%   'xcg_flag': [0|1] use fibers that have CSF and gray-mattter excluded
%     if xcg fibers not found, will exclude those voxels using fname_aseg
%     {default = 0}
%   'M_T1_to_DTI': registration matrix between DTI data and high-res T1 data
%     used to resample fibers to T1 resolution (if resT1flag=1)
%     or to resample aseg to DTI resolution if (resT1flag=0 and xcg_flag=1)
%     {default = eye(4)}
%   'volsz_T1': dimensions of high-res T1 volume
%     {default = [256 256 256]}
%   'M_T1': vox2ras matrix for high-res T1 volume
%     {default = [-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1]}
%   'thresh_prob': fiber probability threshold
%     {default = 0.08}
%   'thresh_FA': fractional anisotropy threshold
%     {default = 0}
%   'fname_FA': full path name of FA image file
%     Required if thresh_FA>0
%     {default = []}
%   'fname_aseg': full path name of aseg file
%     required only if xcg_flag=1 and xcg fibers do not already exist
%     {default = []}
%   'verbose': display status messages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  10/10/12 by Don Hagler
% Last Mod: 04/15/13 Don Hagler
%

% based on DTI_MMIL_AsegMask_Fibers, created 06/03/09 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: allow fiber_infix and fiber_ext as parameters, with [] as default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters and check for problems
parms = check_input(fiber_dir,varargin);

% resample aseg to DTI resolution if necessary
if parms.fseg_xcg_flag
  parms = resamp_aseg(parms);
end;

% check whether output file exists
if exist(parms.fname_fseg,'file') && ~parms.forceflag, return; end;

% load fibers, combine into segmentation file
create_segmentation(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fiber_dir,options)
  parms_filter = {...
    'fiber_dir',fiber_dir,[],...
  ...
    'fname_fseg','fseg.mgz',[],...
    'outdir',pwd,[],...
    'fibers',[101:110,115:120,123,135:138,143:150],[],... % excluded Fmaj (121), Fmin (122), SLF (133,134), SCS (141,142)
    'resT1flag',false,[false true],...
    'M_T1_to_DTI',eye(4),[],...
    'volsz_T1',[256,256,256],[],...
    'M_T1',[-1,0,0,129;0,0,1,-129;0,-1,0,129;0,0,0,1],[],...
    'atlas_flag',2,[0:4],...
    'thresh_prob',0.08,[0 1],...
    'thresh_FA',0,[0 1],...
    'fname_FA',[],[],...
    'xcg_flag',false,[false true],...
    'fname_aseg',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden parameters
    'count_flag',true,[false true],...
    'roicode_multiple',200,[0,Inf],... % base roi code for voxels with multiple fibers
    'roicode_base',10000,[0,Inf],... % value added to all roi codes
    'interpm',1,[0:5],...
    'bclamp',true,[false,true],...
    'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % check if input exists
  if ~exist(fiber_dir,'dir')
    error('fiber_dir %s not found',fiber_dir);
  end;
  % check fname_FA if thresh_FA > 0
  if parms.thresh_FA>0
    if isempty(parms.fname_FA)
      error('must specify fname_FA if thresh_FA>0');
    end;
    if ~exist(parms.fname_FA,'file')
      error('file %s not found',parms.fname_FA);
    end;
  end;
  % set full output file name
  if mmil_isrelative(parms.fname_fseg)
    parms.fname_fseg = [parms.outdir '/' parms.fname_fseg];
  else
    parms.outdir = fileparts(parms.fname_fseg);
  end;

  % check for existing resT1 fibers
  parms.fseg_resT1flag = parms.resT1flag;
  if parms.resT1flag
    [fiber_infix,fiber_ext] = dti_set_fiber_infix(...
      'resT1flag',parms.resT1flag,...
      'count_flag',parms.count_flag,'atlas_flag',parms.atlas_flag);
    flist = dir(sprintf('%s/fiber_*_%s%s',...
      parms.fiber_dir,fiber_infix,fiber_ext));
    if ~isempty(flist)
      parms.fseg_resT1flag = 0; % use existing resampled fibers
      parms.fiber_infix = fiber_infix;
      parms.fiber_ext = fiber_ext;
    end;
  end;
  
  % check for original fibers
  if ~parms.resT1flag || parms.fseg_resT1flag
    % set fiber_infix and fiber_ext
    [parms.fiber_infix,parms.fiber_ext] = dti_set_fiber_infix(...
      'count_flag',parms.count_flag,'atlas_flag',parms.atlas_flag);
    % check fiber files
    flist = dir(sprintf('%s/fiber_*_%s%s',...
      parms.fiber_dir,parms.fiber_infix,parms.fiber_ext));
  end;
  if isempty(flist)
    error('no %s fiber files with fiber infix %s found in %s\n',...
      parms.fiber_ext,parms.fiber_infix,parms.fiber_dir);
  end;
  
  % check for existing xcg fibers
  parms.fseg_xcg_flag = parms.xcg_flag;
  if parms.xcg_flag
    [fiber_infix,fiber_ext] = dti_set_fiber_infix(...
      'resT1flag',parms.resT1flag,'xcg_flag',parms.xcg_flag,...
      'count_flag',parms.count_flag,'atlas_flag',parms.atlas_flag);
    flist = dir(sprintf('%s/fiber_*_%s%s',...
      parms.fiber_dir,fiber_infix,fiber_ext));
    if ~isempty(flist)
      parms.fseg_xcg_flag = 0; % use existing xcg fibers
      parms.fiber_infix = fiber_infix;
      parms.fiber_ext = fiber_ext;
    end;
  end;

  % check fname_aseg if needed
  if parms.fseg_xcg_flag
    if isempty(parms.fname_aseg)
      error('must specify fname_aseg if xcg_flag=1 and xcg fibers do not already exist');
    end;
    if ~exist(parms.fname_aseg,'file')
      error('file %s not found',parms.fname_aseg);
    end;
  end;
  % load header info from a fiber file to get dimensions
  parms.volsz = [];
  parms.M = [];
  parms.nfibers = length(parms.fibers);
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_fiber = sprintf('%s/fiber_%d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    if ~exist(fname_fiber,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_fiber);
      continue;
    end;  
    [parms.volsz,parms.M] = read_volsz(fname_fiber,parms);
    break;
  end;
  if isempty(parms.volsz)
    error('no %s fiber files with fiber infix %s found in %s\n',...
      parms.fiber_ext,fiber_infix,fiber_dir);
  end;
  if length(parms.volsz_T1)>3
    parms.volsz_T1 = parms.volsz_T1(1:3);
  end;
  if parms.fseg_resT1flag
    parms.volsz = parms.volsz_T1;
    parms.M = parms.M_T1;
  end;
  parms.nvox = prod(parms.volsz);
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_aseg(parms)
  % compare size and vox2ras matrix of fibers and aseg
  [volsz_aseg,M_aseg] = read_volsz(parms.fname_aseg,parms);
  if all(volsz_aseg==parms.volsz) && all(M_aseg(:)==parms.M(:))
    % no resampling needed
    return;
  end;
  if parms.resT1flag
    error('resT1flag=1 but fname_aseg not in T1 resolution');
  end;
  [tmp,fstem,fext] = fileparts(parms.fname_aseg);
  fname_aseg_resDTI = sprintf('%s/%s_resDTI%s',...
    parms.outdir,fstem,fext);
  if ~exist(fname_aseg_resDTI,'file') || parms.forceflag
    [vol,M] = fs_load_mgh(parms.fname_aseg);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',parms.volsz,'M_ref',parms.M,...
      'interpm',0,...
      'M_reg',inv(parms.M_T1_to_DTI));
    fs_save_mgh(vol,fname_aseg_resDTI,M);
  end;
  parms.fname_aseg = fname_aseg_resDTI;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_segmentation(parms)
  % load FA file if needed and resample to T1 resolution
  if parms.thresh_FA>0
    if parms.verbose
      fprintf('%s: loading FA file %s...\n',mfilename,parms.fname_FA);
    end;
    vec_FA = load_DTI_vec(parms.fname_FA,parms);
  end;

  % load aseg file if needed
  if parms.fseg_xcg_flag
    if parms.verbose
      fprintf('%s: loading aseg file %s...\n',mfilename,parms.fname_aseg);
    end;
    vec_aseg = load_vec(parms.fname_aseg,parms);
  end;

  % load fiber files
  fibercount = zeros(parms.nvox,1);
  fibernums = zeros(parms.nvox,1);
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_fiber = sprintf('%s/fiber_%02d_%s%s',...
      parms.fiber_dir,fnum,parms.fiber_infix,parms.fiber_ext);
    if ~exist(fname_fiber,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_fiber);
      continue;
    end;
    if parms.verbose
      fprintf('%s: loading fiber file %s...\n',mfilename,fname_fiber);
    end;
    vec_fiber = load_DTI_vec(fname_fiber,parms);

    % apply probability threshold
    if parms.verbose
      fprintf('%s: applying probability threshold %0.3f...\n',...
        mfilename,parms.thresh_prob);
    end;
    vec_fiber = 1.0*(vec_fiber>=parms.thresh_prob);

    % apply FA threshold
    if parms.thresh_FA>0
      if parms.verbose
        fprintf('%s: applying FA threshold %0.3f...\n',...
          mfilename,parms.thresh_FA);
      end;
      vec_fiber(vec_FA<parms.thresh_FA) = 0;
    end;

    fibercount = fibercount + vec_fiber;
    fibernums = fibernums + vec_fiber*fnum;
  end;

  clear vec_fiber vec_FA;

  % find voxels with 1 or more fibers
  if parms.verbose
    fprintf('%s: finding voxels with 1 or more fibers...\n',mfilename);
  end;
  ind_one = find(fibercount==1);
  ind_multiple = find(fibercount>1);

  % create segmentation volume
  if parms.verbose
    fprintf('%s: creating segmentation volume...\n',mfilename);
  end;
  vec_fseg = zeros(parms.nvox,1);
  vec_fseg(ind_one) = fibernums(ind_one) + parms.roicode_base;
  vec_fseg(ind_multiple) = fibercount(ind_multiple) +...
    parms.roicode_multiple + parms.roicode_base;

  % exclude CSF and gray matter
  if parms.fseg_xcg_flag
    ind = find(ismember(vec_aseg,parms.xcg_codes));
    vec_fseg(ind) = 0;
  end;

  % save segmentation volume
  vol_fseg = reshape(vec_fseg,parms.volsz);
  fs_save_mgh(vol_fseg,parms.fname_fseg,parms.M);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M] = resample_vol(vol,M,parms)
  % resample vol using M_T1_to_DTI
  if parms.verbose
    fprintf('%s: resampling volume to T1 resolution...\n',mfilename);
  end;
  [vol,M] = mmil_resample_vol(vol,M,...
    'nvox_ref',parms.volsz_T1,'M_ref',parms.M_T1,...
    'interpm',parms.interpm,'bclamp',parms.bclamp,...
    'M_reg',parms.M_T1_to_DTI);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vec = load_DTI_vec(fname,parms)
  vec = [];
  [vol,M] = load_vol(fname);
  if parms.fseg_resT1flag
    vol = resample_vol(vol,M,parms);
  end;
  vec = reshape(vol,[parms.nvox,1]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vec = load_vec(fname,parms)
  vec = [];
  [vol,M] = load_vol(fname);
  vec = reshape(vol,[parms.nvox,1]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz,M] = read_volsz(fname,parms)
  volsz = []; M = [];
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [tmp,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
  end;
  volsz = volsz(1:3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M,volsz] = load_vol(fname)
  vol = []; M = []; volsz = [];
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [vol,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [vol,M,tmp,volsz] = fs_load_mgh(fname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_vol(vol,fname,M)
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      mmil_save_sparse(vol,fname,M);
    case {'.mgh','.mgz'}
      fs_save_mgh(vol,fname,M);
  end;
return;

