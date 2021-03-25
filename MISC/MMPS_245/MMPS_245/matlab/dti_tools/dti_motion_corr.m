function [fname_out,fname_qmat] = dti_correct_motion(fname,varargin)
%function [fname_out,fname_qmat] = dti_correct_motion(fname,[options])
%
% Purpose: Correct head motion for DTI scans using
%            rigid-body volume registration between data and tensor fit
%
% Usage:
%  [fname_out,fname_qmat] = dti_correct_motion(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Parameters:
%  'fname_out': output file name
%    if empty, will add '_mc' to file stem of input file name
%    {default = []}
%  'fname_qmat': output file name containing adjusted qmat
%     and motion estimates
%    if empty, will add '_mc_qmat' to file stem of input file name
%    {default = []}
%  'interpm' - interpolation method
%     0 = nearest neighbor, 1 = linear, 2 = cubic
%     3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%    { default: 2 }
%  'maskoutput': apply brain mask to set background to zero
%    {default = 0}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': overwrite existing output files
%     {default: 0}
%
% Optional Parameters for motion estimation:
%  'intra_flag': perform within-scan motion correction
%     {default = 1}
%  'qmat': matrix of diffusion direction vectors
%    if empty, will do between-scan estimation only
%    {default = []}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    if empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%    if empty, will generate one from fname_ref using mmil_quick_brainmask
%    {default = []}
%  'fname_reg': output file name for between-scan registration matrix
%    {default = []}
%  'censor_niter': number of iterations of censoring
%    {default = 1}
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
%
% Output:
%   fname_out: motion correct mgh/mgz file containing 4D diffusion volume
%   fname_qmat: mat file containing:
%    'qmat': diffusion direction matrix (ndirs x 3) adjusted for rotation
%      from motion correction
%    'M_motion': cell array containing registration matrices M_ref_to_orig
%      for each frame
%
% Created:  02/11/13 by Don Hagler
% Last Mod: 03/23/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_qmat = [];
if ~mmil_check_nargs(nargin,2), return; end;

parms = check_input(fname,varargin);

if output_exists(parms)
  fname_out = parms.fname_out;
  fname_qmat = parms.fname_qmat;
  return;
end;

[parms,fname_qmat] = estimate_motion(parms);

fname_out = correct_motion(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'fname_out',[],[],...
    'fname_qmat',[],[],...
    'maskoutput',false,[false true],...
    'interpm',2,[0:5],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
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
  ... % hidden
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'smf',1e-5,[1e-100,1e-1],...
    'motion_radius',50,[],... % for calculating distance from angle
  ...
    'estimate_tags',{'intra_flag','qmat',...
                     'fname_ref','fname_mask','fname_reg','censor_niter',...
                     'censor_min_ndirs','censor_thresh','censor_mat',...
                     'nonlin_flag','b0_thresh',...
                     'bvals','min_trans','min_rot','verbose',...
                     'forceflag','mstep','scales','smf'},[],...
    'censor_tags',{'volmask','censor_min_ndirs','nb0'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;

  [fpath,fstem,fext] = fileparts(parms.fname_in);
  if ~isempty(parms.fname_out)
    parms.outdir = fileparts(parms.fname_out);
  else
    parms.outdir = fpath;
  end;
  if isempty(parms.fname_out)
    parms.fname_out = [parms.outdir '/' fstem '_mc' fext];
  end;
  if isempty(parms.fname_qmat)
    parms.fname_qmat = [parms.outdir '/' fstem '_mc_qmat.mat'];
  end;

  % get info about input volume
  [parms.M,parms.volsz] = ...
    mmil_load_mgh_info(parms.fname_in,parms.forceflag);
  parms.nframes = parms.volsz(4);

  % get info about input volume
  [parms.M_ref,parms.volsz_ref] = ...
    mmil_load_mgh_info(parms.fname_ref,parms.forceflag);

  % create brainmask if needed
  if parms.maskoutput
    if ~isempty(parms.fname_mask)
      [fpath,fstem_ref,fext] = fileparts(parms.fname_ref);
      parms.fname_mask = [parms.outdir '/' fstem_ref '_mask' fext];
    end;
    if ~exist(parms.fname_mask,'file') || parms.forceflag
      mmil_quick_brainmask([],...
        'fname_in',parms.fname_ref,'fname_out',parms.fname_mask);
    end;
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
  if ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if ~exist(parms.fname_qmat,'file')
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname_qmat] = estimate_motion(parms)
  if ~exist(parms.fname_qmat,'file') || parms.forceflag
    args = mmil_parms2args(parms,parms.estimate_tags);
    [M_motion,motion_tseries,qmat,DTfit] = ...
      dti_estimate_motion(parms.fname_in,args{:});
    if ~isempty(parms.qmat)
      bvals = parms.bvals;
      orig_qmat = parms.qmat;
      % calculate absolute difference in head position or rotation
      [mean_motion,mean_trans,mean_rot] = ...
        mmil_mean_motion(motion_tseries,parms.motion_radius);
      save(parms.fname_qmat,'M_motion','motion_tseries',...
        'mean_motion','mean_trans','mean_rot',...
        'qmat','orig_qmat','bvals','DTfit');
    end;
  else
    load(parms.fname_qmat);
  end;
  parms.M_motion = M_motion;
  parms.DTfit = DTfit;
  fname_qmat = parms.fname_qmat;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = correct_motion(parms)
  if ~exist(parms.fname_out,'file') || parms.forceflag
    [vol_orig,M] = fs_load_mgh(parms.fname_in);
    % replace censored slices in vol with those synthesized from tensor
    if parms.censor_niter>0 && ~isempty(parms.DTfit)
      vol_synth = dti_synth_vol(parms.DTfit);
      parms.volmask = parms.DTfit.volmask;
      parms.nb0 = length(parms.DTfit.i_b0);
      args = mmil_parms2args(parms,parms.censor_tags);
      vol_orig = dti_censor_vol(vol_orig,vol_synth,...
                                parms.DTfit.censor_mat,args{:});
    end;
    % initialize output
    vol = zeros([parms.volsz_ref(1:3),parms.nframes],'single');
    if parms.maskoutput
      volmask = fs_load_mgh(parms.fname_mask);
    end;
    if parms.verbose
      if parms.nframes>1
        plurstr = 's';
      else
        plurstr = [];
      end;      
      fprintf('%s: resampling %d frame%s...\n',...
        mfilename,parms.nframes,plurstr);
    end;
    for f = 1:parms.nframes
      vol_tmp = squeeze(vol_orig(:,:,:,f));
      vol_tmp = mmil_resample_vol(vol_tmp,M,...
        'M_ref',parms.M_ref,'nvox_ref',parms.volsz_ref,...
        'M_reg',parms.M_motion{f},...
        'interpm',parms.interpm);
      if parms.maskoutput
        vol_tmp(~volmask) = 0;
      end;
      vol(:,:,:,f) = vol_tmp;
    end;
    fs_save_mgh(vol,parms.fname_out,parms.M_ref);
  end;
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

