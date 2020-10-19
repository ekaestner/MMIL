function [fname_out,fname_B0dx,fname_motion_B0uw,fname_regref_B0uw] =...
            bold_unwarpB0(fname,varargin)
%function [fname_out,fname_B0dx,fname_motion_B0uw,fname_regref_B0uw] =...
%            bold_unwarpB0(fname,[options])
%
% Purpose: correct BOLD images for distortions due to B0 inhomogeneity
%
% Required Parameters:
%   fname: full path name of mgh/mgz file containing 4D BOLD volume
%     if empty, will estimate B0 distortion only
%
% Optional Parameters:
%   'fname_out': output file fname for corrected data
%     If empty, will append '_B0uw.mgz' to file stem of fname
%     {default = []}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': overwrite existing output files
%     {default: 0}
%
% Optional Parameters for B0 unwarping (uses Dominic Holland's unwarpB0):
%   'optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%     {default = 0}
%   'fname_B0dx'; name of input/output B0dx file with estimated displacements
%     If not found, will create using fname_for and fname_rev
%     If empty, will not do B0 unwarping
%     {default = []}
%   'fname_for'; name of input b=0 image with "forward" phase-encode polarity
%      Required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'fname_rev'; name of input b=0 image with "reverse" phase-encode polarity
%      Required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'EchoSpacing': echo spacing for fname
%     If supplied and different from EchoSpacing_for, will scale dx field
%     {default = []}
%   'AcquisitionColumns': acquisition columns for fname
%     If supplied and different from EchoSpacing_for, will scale dx field
%     {default = []}
%   'EchoSpacing_for': echo spacing for fname_for
%     {default = []}
%   'AcquisitionColumns_for': acquisition columns for fname_for
%     {default = []}
%   'EchoSpacing_rev': echo spacing for fname_rev
%     If supplied and differenct from EchoSpacing_for, will give warning
%     {default = []}
%   'AcquisitionColumns_rev': acquisition columns for fname_rev
%     If supplied and differenct from EchoSpacing_for, will give warning
%     {default = []}
%   'PhaseDir': phase encode direction (obtained from dicom header)
%     'COL' or 'ROW'
%     {default = 'COL'}
%   'revflag': [0|1] whether to treat input data as "reverse" 
%     phase-encode polarity for applying B0 unwarping
%     {default = 0}
%
% Optional Parameters for estimating motion:
%  'regref_B0uw_flag': [0|1] whether to register reference image (fname_ref)
%      to fname_for (if revflag = 0) or fname_rev (if revflag = 1)
%      and apply transformation to B0dx
%    {default = 1}
%  'motion_B0uw_flag': [0|1] whether to estimate head motion and apply to B0dx
%    {default = 0}
%  'motion_B0uw_iters': number of iterations to estimate motion and B0 displacement
%    {default = 2}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    If empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%  'fname_reg': output file name for between-scan registration matrix
%     {default = []}
%  'ref_frame': reference frame number (1-based)
%    {default = 1}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'fname_motion_B0uw': output file name of mat file containing estimated motion
%    If empty, will append '_motion_B0uw.mat' to file stem of fname
%    Ignored if motion_B0uw_flag = 0
%    {default = []}
%  'fname_regref_B0uw': output file name of mat file containing registration
%     between fname_ref and fname_for or fname_ref
%    Ignored if regref_B0uw_flag = 0
%    If empty, will append '_regreg_B0uw.mat' to file stem of fname
%    {default = []}
%
% Output:
%   fname_out: mgz file containing 4D BOLD volume corrected for B0 distortion
%   fname_B0dx: mgz file containing volume of estimated displacement
%     along phase encode direction
%   fname_motion_B0uw: mat file containing estimated motion applied to B0dx:
%     'M_motion': cell array of 4x4 transformation matrices defining
%       registration from reference to each frame
%     'motion_tseries': matrix of motion estimates [nframes,6]
%       with this order: rotz, rotx, roty, dz, dx, dy
%                        to match AFNI 3dvolreg's motion.1D files
%           also called: roll, pitch, yaw, dS, dL, dP
%   fname_regref_B0uw: mat file containing registration matrix
%     'M_reg': cell array of 4x4 transformation matrices defining
%       registration from reference to fname_for or fname_rev
%
% Created:  08/28/12 by Don Hagler
% Last Mod: 07/24/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_B0dx = []; fname_motion_B0uw = []; fname_regref_B0uw = [];

if ~mmil_check_nargs(nargin, 1), return; end;

if ~isempty(fname) && ~exist(fname,'file'), error('file %s not found',fname); end;

parms = check_input(fname,varargin);

fname_out = parms.fname_out;
fname_B0dx = parms.fname_B0dx;
fname_motion_B0uw = parms.fname_motion_B0uw;
fname_regref_B0uw = parms.fname_regref_B0uw;

if output_exists(parms), return; end;

% estimate B0 distortion
parms = estimate_distortion(parms);

if isempty(parms.fname_in), return; end;

% estimate motion between reference scan and field map
if parms.regref_B0uw_flag
  parms = register_reference(parms);
end;

% scale displacements according to differences in echo spacing, etc.
parms = scale_displacements(parms);

if parms.motion_B0uw_flag
  % estimate motion
  parms = estimate_motion(parms,1);
  for i=2:parms.motion_B0uw_iters
    % apply B0 correction to reference scan
    if ~strcmp(parms.fname_in,parms.fname_ref)
      parms = apply_correction_to_ref(parms);
    end;
    % apply B0 correction, correcting for motion
    parms = apply_correction(parms,0); % temporary output
    % re-estimate motion using B0 corrected data
    parms = estimate_motion(parms,0); % use temporary output
  end;
  % save motion estimate to mat file
  save_motion(parms);
end;

parms = apply_correction(parms);

% remove temporary files
cleanup_tmpdir(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
    'fname_out',[],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % B0 unwarping
    'optimize_B0uw_flag',false,[false true],...
    'fname_B0dx',[],[],...
    'fname_for',[],[],...
    'fname_rev',[],[],...
    'revflag',false,[false,true],...
    'EchoSpacing',[],[],...
    'AcquisitionColumns',[],[],...
    'AcquisitionRows',[],[],...
    'PhaseDir','COL',{'ROW','COL'},...
    'EchoSpacing_for',[],[],...
    'AcquisitionColumns_for',[],[],...
    'AcquisitionRows_for',[],[],...
    'EchoSpacing_rev',[],[],...
    'AcquisitionColumns_rev',[],[],...
    'AcquisitionRows_rev',[],[],...
  ... % motion correction for B0uw
    'regref_B0uw_flag',true,[false true],...
    'mask_thresh',0.82,[0,1],...
    'motion_B0uw_flag',false,[false true],...
    'motion_B0uw_iters',2,[1:10],...
    'fname_ref',fname,[],...
    'fname_mask',[],[],...
    'ref_frame',1,[1 Inf],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'fname_motion_B0uw',[],[],...
    'fname_regref_B0uw',[],[],...
    'M_reg',[],[],...
  ... % B0uw optimization
    'kernelWidthMax',25,[1:100],...
    'lambda2',1100,[1:10000],...
    'kernelWidthMax_vec',[25,31,35],[1:100],...
    'lambda2_vec',[1100,1500,1900],[1:10000],...
    'multi_opt_flag',false,[false true],...
  ... % other
    'outfix','B0uw',[],...
    'outext','.mgz',{'.mgh','.mgz'},...
    'tmpext','.mgh',{'.mgh','.mgz'},...
    'cleanupflag',true,[false true],...
  ...
    'motion_tags',{'intra_flag','fname_ref','ref_frame','fname_mask',...
                   'min_trans','min_rot','verbose','mstep','scales'},[],...
  };

  parms = mmil_args2parms(options,parms_filter);
  % determine whether to swap rows and columns
  if strcmp(parms.PhaseDir,'ROW')
    parms.swapxy_flag = 1;
  else
    parms.swapxy_flag = 0;
  end;
  % set fname_out
  if isempty(parms.fname_in)
    parms.fname_out = [];
    parms.fname_motion_B0uw = [];
    parms.fname_regref_B0uw = [];
    if isempty(parms.fname_B0dx)
      parms.outdir = fileparts(parms.fname_for);
    else
      parms.outdir = fileparts(parms.fname_B0dx);
    end;
  else
    [tpath,tstem,text] = fileparts(parms.fname_in);
    if isempty(parms.fname_out)
      parms.fname_out = [tpath '/' tstem '_' parms.outfix parms.outext];
    end;
    parms.outdir = fileparts(parms.fname_out);
    parms.tmpdir = [parms.outdir '/tmp_B0uw'];
    % set fname_motion_B0uw
    if parms.motion_B0uw_flag && isempty(parms.fname_motion_B0uw)
      parms.fname_motion_B0uw = ...
        [parms.outdir '/' tstem '_motion_B0uw.mat'];
    elseif ~parms.motion_B0uw_flag
      parms.fname_motion_B0uw = [];
    end;
    % set fname_regref_B0uw
    if parms.regref_B0uw_flag && isempty(parms.fname_regref_B0uw)
      parms.fname_regref_B0uw = ...
        [parms.outdir '/' tstem '_regref_B0uw.mat'];
    elseif ~parms.regref_B0uw_flag
      parms.fname_regref_B0uw = [];
    end;
    % get info about input file
    [parms.M,parms.volsz] = ...
      mmil_load_mgh_info(parms.fname_in,parms.forceflag,parms.outdir);
  end;
  % check that fname_for and fname_rev exist
  if isempty(parms.fname_B0dx) || ~exist(parms.fname_B0dx,'file')
    if isempty(parms.fname_for) || isempty(parms.fname_rev)
      error('either fname_B0dx or fname_for and fname_rev must be supplied');
    end;
    if ~exist(parms.fname_for,'file')
      error('file %s not found',parms.fname_for);
    end;
    if ~exist(parms.fname_rev,'file')
      error('file %s not found',parms.fname_rev);
    end;
    % warn about differences in acquisition between for and rev
    check_acquisition_diffs(parms);
  end;
  % set fname_B0dx
  if isempty(parms.fname_B0dx)
    require_for_rev_flag = 1;
    [tpath,tstem_rev,text] = fileparts(parms.fname_rev);
    [tpath,tstem_for,text] = fileparts(parms.fname_for);
    parms.fname_B0dx = ...
      [parms.outdir '/' tstem_for 'VS' tstem_rev '_B0dx' parms.outext];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_acquisition_diffs(parms)
  if ~isempty(parms.EchoSpacing_for) &...
     ~isempty(parms.EchoSpacing_rev) &...
     parms.EchoSpacing_for ~= parms.EchoSpacing_rev
    fprintf('%s: WARNING: EchoSpacing does not match (%0.2f vs. %0.2f)\n',...
      mfilename,parms.EchoSpacing_for,parms.EchoSpacing_rev);
  end;
  if ~isempty(parms.AcquisitionColumns_for) &...
     ~isempty(parms.AcquisitionColumns_rev) &...
     parms.AcquisitionColumns_for ~= parms.AcquisitionColumns_rev
    fprintf('%s: WARNING: AcquisitionColumns do not match (%0.2f vs. %0.2f)\n',...
      mfilename,parms.AcquisitionColumns_for,parms.AcquisitionColumns_rev);
  end;
  if ~isempty(parms.AcquisitionRows_for) &...
     ~isempty(parms.AcquisitionRows_rev) &...
     parms.AcquisitionRows_for ~= parms.AcquisitionRows_rev
    fprintf('%s: WARNING: AcquisitionRows do not match (%0.2f vs. %0.2f)\n',...
      mfilename,parms.AcquisitionRows_for,parms.AcquisitionRows_rev);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~isempty(parms.fname_out) && ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if ~exist(parms.fname_B0dx,'file')
    exist_flag = 0;
  end;
  if parms.motion_B0uw_flag &&...
     ~isempty(parms.fname_motion_B0uw) && ~exist(parms.fname_motion_B0uw,'file')
    exist_flag = 0;
  end;
  if parms.regref_B0uw_flag &&...
     ~isempty(parms.fname_regref_B0uw) && ~exist(parms.fname_regref_B0uw,'file')
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = estimate_distortion(parms)
  % estimate B0 distortion
  if ~exist(parms.fname_B0dx,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: estimating B0 distortion from\n%s\nand\n%s...\n',...
        mfilename,parms.fname_for,parms.fname_rev);
    end;
    if parms.optimize_B0uw_flag
      if parms.verbose
        fprintf('%s: optimizing B0 unwarping parameters...\n',mfilename);
        tic
      end;
      [parms.kernelWidthMax,parms.lambda2] = epi_optimizeB0(...
        parms.fname_for,parms.fname_rev,...
        'swapxy_flag',parms.swapxy_flag,...
        'kernelWidthMax_vec',parms.kernelWidthMax_vec,...
        'lambda2_vec',parms.lambda2_vec,...
        'lambda2_init',parms.lambda2,...
        'multi_opt_flag',parms.multi_opt_flag);
      if parms.verbose, toc; end;
    end;
    [tpath,tstem] = fileparts(parms.fname_for);
    parms.fname_for_B0uw = [parms.outdir '/' tstem '_f0_B0uw.mgz'];
    [tpath,tstem] = fileparts(parms.fname_rev);
    parms.fname_rev_B0uw = [parms.outdir '/' tstem '_f0_B0uw.mgz'];
    if parms.verbose, tic; end;
    epi_estimateB0(parms.fname_for,parms.fname_rev,...
      'fname_dx',parms.fname_B0dx,...
      'fname_for_out',parms.fname_for_B0uw,...
      'fname_rev_out',parms.fname_rev_B0uw,...
      'swapxy_flag',parms.swapxy_flag,...
      'kernelWidthMax',parms.kernelWidthMax,...
      'lambda2',parms.lambda2,...
      'forceflag',parms.forceflag);
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = register_reference(parms)
  % check if fname_in is either fname_for or fname_rev
  if ismember(parms.fname_in,{parms.fname_for,parms.fname_rev})
    return;
  end;
  % set file names of volumes to be registered
  if parms.revflag
    fname_SE = parms.fname_rev;
  else
    fname_SE = parms.fname_for;
  end;
  fname_GE = parms.fname_ref;
  [fpath,fstem_SE] = fileparts(fname_SE);
  [fpath,fstem_GE] = fileparts(fname_GE);
  if strcmp(fname_SE,fname_GE)
    parms.M_reg = [];
    if parms.verbose
      fprintf('%s: skipping registration between %s and %s...\n',mfilename,fstem_SE,fstem_GE);
    end;
    return;
  end;
  if parms.verbose
    fprintf('%s: registering reference images %s and %s...\n',mfilename,fstem_SE,fstem_GE);
  end;
  if ~exist(parms.fname_regref_B0uw,'file') || parms.forceflag
    tparms = [];
    tparms.outdir = sprintf('%s/reg_GE_%s_SE_%s',parms.tmpdir,fstem_GE,fstem_SE);
    tparms.mask_thresh = parms.mask_thresh;
    tparms.cleanup_flag = 0;
    tparms.forceflag = 0;
    args = mmil_parms2args(tparms);
    M_GE_to_SE = mmil_jpdfreg_GESE(fname_GE,fname_SE,args{:});
    save(parms.fname_regref_B0uw,'M_GE_to_SE','fname_SE','fname_GE');
  else
    tmp = load(parms.fname_regref_B0uw);
    M_GE_to_SE = tmp.M_GE_to_SE;
  end;
  parms.M_reg = M_GE_to_SE;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = scale_displacements(parms)
  parms.fname_B0dx_orig = parms.fname_B0dx;
  dx_sf = 1;
  % longer echo spacing results in larger distortions
  if ~isempty(parms.EchoSpacing_for) &...
     ~isempty(parms.EchoSpacing) &...
     parms.EchoSpacing_for>0 &...
     parms.EchoSpacing>0 &...
     parms.EchoSpacing_for ~= parms.EchoSpacing
    dx_sf = dx_sf*parms.EchoSpacing/parms.EchoSpacing_for;
  end;
  % more phase encode steps results in larger distortions
  if strcmp(parms.PhaseDir,'ROW')
    if ~isempty(parms.AcquisitionRows_for) &...
       ~isempty(parms.AcquisitionRows) &...
       parms.AcquisitionRows_for>0 &...
       parms.AcquisitionRows>0 &...
       parms.AcquisitionRows_for ~= parms.AcquisitionRows
      dx_sf = dx_sf*parms.AcquisitionRows/parms.AcquisitionRows_for;
    end;
  else
    if ~isempty(parms.AcquisitionColumns_for) &...
       ~isempty(parms.AcquisitionColumns) &...
       parms.AcquisitionColumns_for>0 &...
       parms.AcquisitionColumns>0 &...
       parms.AcquisitionColumns_for ~= parms.AcquisitionColumns
      dx_sf = dx_sf*parms.AcquisitionColumns/parms.AcquisitionColumns_for;
    end;
  end;
  % displacements are in units of voxels, so larger voxel size requires
  %   smaller displacements
  M_dx = fs_read_header(parms.fname_B0dx);
  voxsize_dx = fs_voxsize(M_dx);
  voxsize = fs_voxsize(parms.M);
  if voxsize_dx(1)~=voxsize(1)
    voxratio = voxsize_dx(1)/voxsize(1);
    dx_sf = dx_sf*voxratio;
  else
    voxratio = 1;
  end;
  if dx_sf~=1 || voxratio~=1 || ~isempty(parms.M_reg)
    [tpath,tstem] = fileparts(parms.fname_in);
    [tpath_dx,tstem_dx] = fileparts(parms.fname_B0dx);
    fname_B0dx_scaled = sprintf('%s/%s_scaled_%s.mgz',...
      tpath_dx,tstem_dx,tstem);
    if ~exist(fname_B0dx_scaled,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: scaling B0 displacement field...\n',mfilename);
      end;
      if voxratio~=1 || ~isempty(parms.M_reg)
        [vol_dx,M_dx] = fs_load_mgh(parms.fname_B0dx);
        vol_dx = dx_sf*vol_dx;
        % resample displacement field
        vol_dx_res = mmil_resample_vol(vol_dx,M_dx,...
          'M_reg',parms.M_reg,...
          'M_ref',parms.M,'nvox_ref',parms.volsz(1:3),...
          'bclamp',0);
        % smooth resampled displacement field
        %   if displacement field is lower res
        if voxratio>1
          sigma = 1.5*voxratio; % 1.5 is arbitrary, but seems to work (and 1.0 is not enough)
          vol_dx_res = mmil_smooth3d(vol_dx_res,...
            sigma,sigma,sigma);
        end;
        fs_save_mgh(vol_dx_res,fname_B0dx_scaled,parms.M);
      else
        fs_mri_convert(parms.fname_B0dx,fname_B0dx_scaled,...
          'options',sprintf('-sc %f',dx_sf));
      end;
    end;
    parms.fname_B0dx = fname_B0dx_scaled;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = estimate_motion(parms,init_flag)
  % estimate head motion
  if ~exist('init_flag','var'), init_flag = 1; end;
  tparms = parms;
  if init_flag
    fname_in = parms.fname_in;
  else
    % re-estimate motion using B0 corrected images
    fname_in = parms.fname_out_tmp;
    if strcmp(parms.fname_in,parms.fname_ref)
      tparms.fname_ref = fname_in; % no between-scan motion correction
    else
      tparms.fname_ref = parms.fname_ref_B0uw;
    end;
  end;
  % motion estimation
  if parms.verbose
    fprintf('%s: estimating motion for %s with %s as reference...\n',...
      mfilename,fname_in,tparms.fname_ref);
  end;
  args = mmil_parms2args(tparms,parms.motion_tags);
  [parms.M_motion,parms.motion_tseries] = ...
    bold_estimate_motion(fname_in,args{:});
  parms.nframes = length(parms.M_motion);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_motion(parms)
  M_motion = parms.M_motion;
  motion_tseries = parms.motion_tseries;
  save(parms.fname_motion_B0uw,'M_motion','motion_tseries');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = apply_correction_to_ref(parms)
  if parms.verbose
    fprintf('%s: applying B0 distortion correction to reference image %s...\n',...
      mfilename,parms.fname_ref);
  end;
  [tpath,tstem] = fileparts(parms.fname_ref);
  [M,volsz] = mmil_load_mgh_info(parms.fname_ref,parms.forceflag,parms.outdir);
  if volsz(4)>1
    fname_in = [parms.outdir '/' tstem '_f0' parms.outext];
    fname_out = [parms.outdir '/' tstem '_f0_B0uw' parms.outext];
    if ~exist(fname_in,'file') || parms.forceflag
      [vol,M]=fs_load_mgh(parms.fname_ref,[],1);
      fs_save_mgh(vol,fname_in,M);
    end;
  else
    fname_in = parms.fname_ref;
    fname_out = [parms.outdir '/' tstem '_B0uw' parms.outext];
  end;
  epi_unwarpB0(fname_in,parms.fname_B0dx,...
    'fname_out',fname_out,...
    'revflag',parms.revflag,...
    'swapxy_flag',parms.swapxy_flag,...
    'forceflag',parms.forceflag);
  parms.fname_ref_B0uw = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = apply_correction(parms,final_flag)
  % apply B0 distortion correction displacements
  if ~exist('final_flag','var'), final_flag = 1; end;
  if final_flag
    fname_out = parms.fname_out;
    forceflag = parms.forceflag;
  else
    [tpath,tstem] = fileparts(parms.fname_out);
    fname_out = [parms.tmpdir '/' tstem parms.tmpext];
    parms.fname_out_tmp = fname_out;
    forceflag = 1; % iterative, so force overwrite
  end;
  if ~exist(fname_out,'file') || forceflag
    if parms.verbose
      fprintf('%s: applying B0 distortion correction to %s...\n',...
        mfilename,parms.fname_in);
    end;
    if ~parms.motion_B0uw_flag,
      if parms.verbose, tic; end;
      epi_unwarpB0(parms.fname_in,parms.fname_B0dx,...
        'fname_out',fname_out,...
        'revflag',parms.revflag,...
        'swapxy_flag',parms.swapxy_flag,...
        'forceflag',parms.forceflag);
      if parms.verbose, toc; end;
    else
      mmil_mkdir(parms.tmpdir);
      [tpath,tstem_in] = fileparts(parms.fname_in);
      [tpath,tstem_dx] = fileparts(parms.fname_B0dx);
      [tpath,tstem_out] = fileparts(parms.fname_out);
      % save each frame of fname_in to temporary files
      write_frames(parms);
      vol_out = zeros(parms.volsz);
      vol_dx = [];
      if parms.verbose, tic; end;
      for f=1:parms.nframes
        M_reg = parms.M_motion{f};
        fname_in_tmp = sprintf('%s/%s_f%d%s',...
          parms.tmpdir,tstem_in,f,parms.tmpext);
        fname_dx_tmp = sprintf('%s/%s_f%d%s',...
          parms.tmpdir,tstem_dx,f,parms.tmpext);
        fname_out_tmp = sprintf('%s/%s_f%d%s',...
          parms.tmpdir,tstem_out,f,parms.tmpext);
        if isempty(vol_dx)
          [vol_dx,M_dx] = fs_load_mgh(parms.fname_B0dx);
        end;
        % resample using inverse of registration matrix
        vol_dx_res = mmil_resample_vol(vol_dx,M_dx,...
          'M_reg',inv(M_reg),'bclamp',0);
        fs_save_mgh(vol_dx_res,fname_dx_tmp,M_dx);
        % apply distortion correction
        epi_unwarpB0(fname_in_tmp,fname_dx_tmp,...
          'fname_out',fname_out_tmp,...
          'revflag',parms.revflag,...
          'swapxy_flag',parms.swapxy_flag,...
          'forceflag',parms.forceflag);
        [vol_out_tmp,M] = fs_load_mgh(fname_out_tmp);
        vol_out(:,:,:,f) = vol_out_tmp;
      end;
      if parms.verbose, toc; end;
      fs_save_mgh(vol_out,fname_out,parms.M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = write_frames(parms)
  all_exist = 1;
  [tpath,tstem_in] = fileparts(parms.fname_in);
  if ~parms.forceflag
    for f=1:parms.nframes
      fname_in_tmp = sprintf('%s/%s_f%d%s',...
        parms.tmpdir,tstem_in,f,parms.tmpext);
      if ~exist(fname_in_tmp,'file')
        all_exist = 0;
        break;
      end;  
    end;
    if all_exist, return; end;
  end;
  vol = fs_load_mgh(parms.fname_in);
  for f=1:parms.nframes
    fname_in_tmp = sprintf('%s/%s_f%d%s',...
      parms.tmpdir,tstem_in,f,parms.tmpext);
    if ~exist(fname_in_tmp,'file')
      fs_save_mgh(vol(:,:,:,f),fname_in_tmp,parms.M);
    end;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanupflag && exist(parms.tmpdir,'dir')
    % delete tmp dir
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result] = unix(cmd);
    if status
      error('deletion of temporary directory %s failed:\n%s\n',...
        parms.tmpdir,result);
    end;
  end;
return;

