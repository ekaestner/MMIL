function [M_motion,motion_tseries,qmat,DTfit] = dti_estimate_motion(fname,varargin)
%function [M_motion,motion_tseries,qmat,DTfit] = dti_estimate_motion(fname,[options])
%
% Purpose: Estimate head motion for DTI scans using
%            rigid-body volume registration between data and tensor fit
%
% Usage:
%  [M_motion,motion_tseries,qmat,DTfit] = dti_estimate_motion(fname,...
%                                                       'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Input:
%  'intra_flag': estimate within-scan motion
%     {default = 1}
%  'qmat': matrix of diffusion direction vectors
%    if empty, will do between-scan estimation only
%    {default = []}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    if empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%  'fname_reg': output file name for between-scan registration matrix
%     {default = []}
%  'censor_niter': number of iterations of censoring
%     {default = 1}
%  'censor_thresh': error threshold for censoring bad frames
%    normalized to median error for each slice
%    higher values mean less censoring
%    {default = 3.2}
%  'censor_min_ndirs': minimum number of diffusion directions (not including
%    b=0 images) required for tensor fit after censoring
%    will not do censoring if it means reducing number of directions below min
%    {default = 6}
%  'censor_mat': matrix of slices by frame to exclude from fit
%    ignored if censor_niter = 0
%    {default = []}
%  'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 0}
%  'b0_thresh': threshold used for considering a b-value to be 0
%    {default = 10}
%  'bvals': vector of b values (one for each diffusion direction)
%    If single value supplied, will use same for all
%    {default = [1000]}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': overwrite existing output files (e.g. fname_mask, fname_reg)
%     {default: 0}
%
% Output:
%   M_motion: cell array containing 4x4 rigid-body transformation matrices
%             M_ref_to_orig for each frame
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%   qmat: diffusion direction matrix adjusted for motion
%   DTfit: struct containing tensor fit output
%
% Created:  02/11/13 by Don Hagler
% Last Mod: 03/23/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_motion = []; motion_tseries = []; qmat = []; DTfit = [];
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters, check files exist
parms = check_input(fname,varargin);

% between-scan registration
parms = register_scans(parms);

% within-scan registration
if ~isempty(parms.qmat)
  % registration between data and tensor fit
  parms = load_vol(parms);
  if parms.intra_flag
    [parms,DTfit] = tensor_fit(parms);
  end;
  [M_motion,qmat] = estimate_motion(parms);
else
  % use between-scan motion estimate only
  M_motion = set_M_motion(parms);
end;

% derive motion parameter time series from registration matrices
motion_tseries = set_motion_tseries(parms,M_motion);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'intra_flag',true,[false true],...
    'qmat',[],[],...
    'fname_ref',fname,[],...
    'fname_mask',[],[],...
    'fname_reg',[],[],...
    'censor_niter',1,[0,100],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
    'censor_mat',[],[],...
    'nonlin_flag',false,[false true],...
    'b0_thresh',10,[0,500],...
    'bvals',1000,[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'smf',1e-5,[1e-100,1e-1],...
  ...
    'tensor_tags',{'bvals','M','censor_niter','censor_min_ndirs',...
                   'censor_thresh','censor_mat','nonlin_flag','b0_thresh'},[],...
    'censor_tags',{'volmask','censor_min_ndirs','nb0'},[],...
    'reg_tags',{'fname_mask','fname_reg','rigid_flag','affine_flag',...
                'forceflag','tmpdir','cleanupflag',},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;

  % get info about input volume
  [parms.M,parms.volsz] = ...
    mmil_load_mgh_info(parms.fname_in,parms.forceflag);
  parms.nframes = parms.volsz(4);

  % extend bvals to nf elements if only 1 supplied
  if numel(parms.bvals)==1
    bvals = parms.bvals*ones(parms.nframes,1);
  elseif numel(parms.bvals)==parms.nframes
    bvals = reshape(parms.bvals,[parms.nframes,1]);
  else
    error('number of bvals must match frames in vol');
  end;
  
  % check dimensions of qmat
  if isempty(parms.qmat)
    parms.intra_flag = 0;
  else
    if size(parms.qmat,1)~=parms.nframes
      error('number of directions in qmat must match frames in vol');
    end;
    if size(parms.qmat,2)~=3
      error('qmat must have 3 columns');
    end;
    qlength = sqrt(sum(parms.qmat.^2,2));
    parms.i_b0 = find(qlength<parms.smf | bvals<=parms.b0_thresh);
    parms.nb0 = length(parms.i_b0);
    parms.i_bx = setdiff([1:parms.nframes],parms.i_b0);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = register_scans(parms)
  if strcmp(parms.fname_in,parms.fname_ref)
    parms.M_reg = eye(4);
  else
    if parms.verbose
      fprintf('%s: registering %s to %s...\n',...
        mfilename,parms.fname_in,parms.fname_ref);
    end;
    % between scan registration
    args = mmil_parms2args(parms,parms.reg_tags);
    parms.M_reg = epi_register_scans(parms.fname_ref,parms.fname_in,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_vol(parms)
  [parms.vol,parms.M] = fs_load_mgh(parms.fname_in);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,DTfit] = tensor_fit(parms)
  if parms.verbose
    fprintf('%s: calculating tensor fit for %s...\n',mfilename,parms.fname_in);
  end;
  % fit tensor, censoring bad slice/frame combinations
  args = mmil_parms2args(parms,parms.tensor_tags);
  DTfit = dti_fit_tensor(parms.vol,parms.qmat,args{:});
  % first B0 and non-B0 image should not be censored for this motion correction method.
  DTfit.censor_mat(:,parms.i_b0(1)) = 0;
  DTfit.censor_mat(:,parms.i_bx(1)) = 0;
  % synthesize volume using tensor fit
  parms.vol_synth = dti_synth_vol(DTfit);
  % replace censored slices in vol
  if parms.censor_niter>0
    parms.volmask = DTfit.volmask;
    args = mmil_parms2args(parms,parms.censor_tags);
    parms.vol = dti_censor_vol(parms.vol,parms.vol_synth,...
      DTfit.censor_mat,args{:});
  end;
  % use dilated brain mask created by dti_fit_tensor
  parms.volmask = DTfit.volmask_dilated;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_motion,qmat] = estimate_motion(parms)
  M_motion = cell(parms.nframes,1);
  qmat = parms.qmat;
  if parms.intra_flag
    vol_mask = ctx_mgh2ctx(parms.volmask,parms.M);
  end;
  % register each frame to synthesized frame
  %   but keep registration for b=0 images locked to first b=0 frame
  %   and keep registration for b>0 images locked to first b>0 frame
  if parms.verbose && parms.intra_flag
    if parms.nframes>1
      plurstr = 's';
    else
      plurstr = [];
    end;      
    fprintf('%s: registering %d frame%s...\n',...
      mfilename,parms.nframes,plurstr);
  end;
  for f=1:parms.nframes
    if parms.intra_flag
      vol_data = ctx_mgh2ctx(squeeze(parms.vol(:,:,:,f)),parms.M);
      vol_synth = ctx_mgh2ctx(squeeze(parms.vol_synth(:,:,:,f)),parms.M);

      % register images
      M_synth_to_orig = rbreg_vol2vol_icorr_mask(vol_synth,vol_data,...
        vol_mask,0,parms.mstep,0,parms.scales,...
        parms.min_trans,parms.min_rot*pi/180);

      % make registration relative to first b=0 or b>0 image
      if ismember(f,parms.i_b0) % if b=0, register to first b=0 image
        if f==parms.i_b0(1)
          M_ref_to_synth_b0 = inv(M_synth_to_orig);
        end;
        M_ref_to_synth = M_ref_to_synth_b0;
      elseif ismember(f,parms.i_bx) % if b>0, register to first b>0 image
        if f==parms.i_bx(1)
          M_ref_to_synth_bx = inv(M_synth_to_orig);
        end;
        M_ref_to_synth = M_ref_to_synth_bx;
      end
      M_ref_to_orig = M_synth_to_orig*M_ref_to_synth;
    else
      M_ref_to_orig = eye(4);    
    end;
 
     % adjust for between-scan motion
    M_ref_to_orig = M_ref_to_orig * parms.M_reg;

    % store registration in M_motion
    M_motion{f} = M_ref_to_orig;

    % adjust qmat for motion
    qmat(f,:) = dti_rotate_qmat(qmat(f,:),inv(M_ref_to_orig));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_motion = set_M_motion(parms)
  M_motion = cell(parms.nframes,1);
  for f=1:parms.nframes
    M_motion{f} = parms.M_reg;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_tseries = set_motion_tseries(parms,M_motion)
  % derive motion parameter time series from registration matrices
  motion_tseries = zeros(parms.nframes,6);
  for f=1:parms.nframes
    tmp_motion = mmil_M_mat2vec(inv(M_motion{f}));
    motion_tseries(f,:) = tmp_motion;
  end;
return;
