function [M_motion,motion_tseries] = bold_estimate_motion(fname,varargin)
%function [M_motion,motion_tseries] = bold_estimate_motion(fname,[options])
%
% Purpose: Estimate head motion for BOLD scans using
%            rigid-body volume registration
%
% Usage:
%  [M_motion,motion_tseries] = bold_estimate_motion(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Input:
%  'intra_flag': estimate within-scan motion
%     {default = 1}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    If empty, will use fname as reference
%    {default = []}
%  'ref_frame': reference frame number (1-based)
%    {default = 1}
%  'fname_msk': name of mask volume
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'verbose': [0|1] display status messages
%    {default = 1}
%
% Output:
%   M_motion: cell array containing 4x4 rigid-body transformation matrices
%             for each frame
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%
% Created:  08/27/12 by Don Hagler
% Last Mod: 03/23/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_motion = []; motion_tseries = [];
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters, check files exist
parms = check_input(fname,varargin);

% load input volume, reference, and mask
parms = load_data(parms);

% estimate head motion
M_motion = estimate_motion(parms);

% derive motion parameter time series from registration matrices
motion_tseries = set_motion_tseries(parms,M_motion);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'intra_flag',true,[false true],...
    'fname_ref',fname,[],...
    'ref_frame',1,[1 Inf],...
    'fname_mask',[],[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'verbose',true,[false true],...
  ...
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;
  if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file'), 
    error('file %s not found',parms.fname_mask);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_data(parms)
  % load input file
  [parms.vol,parms.M] = fs_load_mgh(parms.fname_in);
  parms.nframes = size(parms.vol,4);

  % check that reference image has enough frames
  if parms.ref_frame>1
    [parms.M_ref,parms.volsz_ref] = fs_read_header(parms.fname_ref);
    if parms.ref_frame > parms.volsz_ref(4)
      error('ref_frame (%d) is greater than number of frames in fname_ref (%d)',...
        parms.ref_frame,parms.volsz_ref(4));
    end;
  end;

  % load reference volume
  if ~strcmp(parms.fname_in,parms.fname_ref)
    [vol_ref,M_ref] =...
      fs_load_mgh(parms.fname_ref,[],parms.ref_frame);
  else
    vol_ref = parms.vol(:,:,:,parms.ref_frame);
    M_ref = parms.M;
  end;
  parms.vol_ref_ctx = ctx_mgh2ctx(vol_ref,M_ref);

  % load mask volume or create from reference volume
  if ~isempty(parms.fname_mask)
    parms.vol_mask_ctx = ctx_load_mgh(parms.fname_mask);
  else
    parms.vol_mask_ctx = mmil_quick_brainmask(parms.vol_ref_ctx);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_motion = estimate_motion(parms)
  M_motion = cell(parms.nframes,1);
  for f=1:parms.nframes
    if f==1 || parms.intra_flag
      vol_ctx = ctx_mgh2ctx(squeeze(parms.vol(:,:,:,f)),parms.M);
      M_reg = rbreg_vol2vol_icorr_mask(parms.vol_ref_ctx,vol_ctx,...
        parms.vol_mask_ctx,0,parms.mstep,0,parms.scales,...
        parms.min_trans,parms.min_rot*pi/180);
    end;
    M_motion{f} = M_reg;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_tseries = set_motion_tseries(parms,M_motion)
  motion_tseries = zeros(parms.nframes,6);
  for f=1:parms.nframes
    tmp_motion = mmil_M_mat2vec(inv(M_motion{f}));
    motion_tseries(f,:) = tmp_motion;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

