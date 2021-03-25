function [fname_out,fname_B0dx,fname_motion_B0uw] = dti_unwarpB0(fname,varargin)
%function [fname_out,fname_B0dx,fname_motion_B0uw] = dti_unwarpB0(fname,[options])
%
% Purpose: correct diffusion images for distortions due to B0 inhomogeneity
%
% Required Parameters:
%   fname: full path name of mgh/mgz file containing 4D DTI volume
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
%   'inorm_B0uw_flag': [0|1] intensity normalize 'forward' and 'reverse'
%     polarity images {default = 0}
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
%   'PhaseDir': phase encode direction (obtained from dicom header)
%     'COL' or 'ROW'
%     {default = 'COL'}
%   'revflag': [0|1] whether to treat input data as "reverse" 
%     phase-encode polarity for applying B0 unwarping
%     {default = 0}
%
% Optional Parameters for estimating motion:
%  'motion_B0uw_flag': [0|1] whether to estimate head motion and apply to B0dx
%    {default = 1}
%  'motion_B0uw_iters': number of iterations to estimate motion and B0 displacement
%    {default = 2}
%  'qmat': matrix of diffusion direction vectors
%    if empty, will do between-scan estimation only
%    {default = []}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    If empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%  'fname_reg': output file name for between-scan registration matrix
%     {default = []}
%  'censor_niter': number of iterations of censoring
%     {default = 1}
%  'censor_min_ndirs': minimum number of diffusion directions (not including
%    b=0 images) required for tensor fit after censoring
%    will not do censoring if it means reducing number of directions below min
%    {default = 6}
%  'censor_thresh': error threshold for censoring bad frames
%    normalized to median error for each slice
%    higher values mean less censoring
%    {default = 3.2}
%  'censor_mat': matrix of slices by frame to exclude from fit
%    ignored if censor_niter = 0
%    {default = []}
%  'bvals': vector of b values (one for each diffusion direction)
%    If single value supplied, will use same for all
%    {default = [1000]}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'fname_motion_B0uw': output file name of mat file containing estimated motion
%    If empty, will append '_motion_B0uw.mat' to file stem of fname
%    {default = []}
%
% Output:
%   fname_out: mgz file containing 4D DTI volume corrected for B0 distortion
%   fname_B0dx: mgz file containing volume of estimated displacement
%     along phase encode direction
%   fname_motion_B0uw: mat file containing estimated motion applied to B0dx:
%     'M_motion': cell array of 4x4 transformation matrices defining
%       registration from reference to each frame
%     'motion_tseries': matrix of motion estimates [nframes,6]
%       with this order: rotz, rotx, roty, dz, dx, dy
%                        to match AFNI 3dvolreg's motion.1D files
%           also called: roll, pitch, yaw, dS, dL, dP
%
% Created:  02/11/13 by Don Hagler
% Last Mod: 10/10/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_B0dx = []; fname_motion_B0uw = [];

if ~mmil_check_nargs(nargin, 1), return; end;

if ~isempty(fname) && ~exist(fname,'file'), error('file %s not found',fname); end;

parms = check_input(fname,varargin);

fname_out = parms.fname_out;
fname_B0dx = parms.fname_B0dx;
fname_motion_B0uw = parms.fname_motion_B0uw;

if output_exists(parms), return; end;

if parms.inorm_B0uw_flag
  parms = inorm_B0uw(parms);
end;

parms = estimate_distortion(parms);

if isempty(parms.fname_in), return; end;

if parms.motion_B0uw_flag
  parms = estimate_motion(parms,1);
  for i=2:parms.motion_B0uw_iters
    parms = apply_correction(parms,0); % temporary output
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
    'inorm_B0uw_flag',true,[false true],...
    'fname_B0dx',[],[],...
    'fname_for',[],[],...
    'fname_rev',[],[],...
    'revflag',false,[false,true],...
    'PhaseDir','COL',{'ROW','COL'},...
  ... % motion correction for B0uw
    'motion_B0uw_flag',true,[false true],...
    'motion_B0uw_iters',2,[1:10],...
    'qmat',[],[],...
    'fname_ref',fname,[],...
    'fname_mask',[],[],...
    'fname_reg',[],[],...
    'censor_niter',1,[0,100],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
    'censor_mat',[],[],...
    'bvals',1000,[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'fname_motion_B0uw',[],[],...
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
    'motion_tags',{'qmat','fname_ref','fname_mask','fname_reg',...
                   'censor_niter','censor_min_ndirs','censor_thresh',...
                   'censor_mat','bvals','min_trans','min_rot',...
                   'interpm','verbose','forceflag'},[],...
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
        [parms.outdir '/' tstem '_' parms.outfix '_motion_B0uw.mat'];
    elseif ~parms.motion_B0uw_flag
      parms.fname_motion_B0uw = [];
    end;
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
  end;
  % set fname_B0dx
  if isempty(parms.fname_B0dx)
    require_for_rev_flag = 1;
    [tpath,tstem_rev,text] = fileparts(parms.fname_rev);
    [tpath,tstem_for,text] = fileparts(parms.fname_for);
    parms.fname_B0dx = ...
      [parms.outdir '/' tstem_for 'VS' tstem_rev '_B0dx' parms.outext];
  end;
  % set fname_for_B0uw
  [tpath,tstem] = fileparts(parms.fname_for);
  parms.fname_for_B0uw = [tpath '/' tstem '_B0uw_f0.mgz'];
  if ~exist(parms.fname_for_B0uw,'file')
    parms.fname_for_B0uw = [parms.outdir '/' tstem '_B0uw_f0.mgz'];
  end;
  % set fname_rev_B0uw
  [tpath,tstem] = fileparts(parms.fname_rev);
  parms.fname_rev_B0uw = [tpath '/' tstem '_B0uw_f0.mgz'];
  if ~exist(parms.fname_rev_B0uw,'file')
    parms.fname_rev_B0uw = [parms.outdir '/' tstem '_B0uw_f0.mgz'];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = inorm_B0uw(parms)
  % set output name
  [fpath,fstem,fext] = fileparts(parms.fname_rev);
  fname_rev_inorm = sprintf('%s/%s_inorm%s',parms.outdir,fstem,fext);
  if ~exist(fname_rev_inorm,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: normalizing intensity of rev phase-encode polarity image...\n',...
        mfilename);
    end;
    parms.fname_rev = epi_inorm_B0uw(parms.fname_for,parms.fname_rev,...
      'fname_rev_inorm',fname_rev_inorm);
  else
    parms.fname_rev = fname_rev_inorm;
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
      tparms.fname_ref = fname_in;
    else
      if ~parms.revflag
        tparms.fname_ref = parms.fname_for_B0uw;
      else
        tparms.fname_ref = parms.fname_rev_B0uw;
      end;
    end;
  end;
  % motion estimation
  if parms.verbose
    fprintf('%s: estimating motion...\n',mfilename);
  end;
  args = mmil_parms2args(tparms,parms.motion_tags);
  [parms.M_motion,parms.motion_tseries] = dti_estimate_motion(fname_in,args{:});
  parms.nframes = length(parms.M_motion);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_motion(parms)
  M_motion = parms.M_motion;
  motion_tseries = parms.motion_tseries;
  save(parms.fname_motion_B0uw,'M_motion','motion_tseries');
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
      fprintf('%s: applying B0 distortion correction...\n',mfilename);
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

