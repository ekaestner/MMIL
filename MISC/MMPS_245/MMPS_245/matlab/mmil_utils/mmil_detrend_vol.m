function vol = mmil_detrend_vol(vol,varargin)
%function vol = mmil_detrend_vol(vol,varargin)
%
% Purpose: normalize, detrend, and regress motion from timeseries volume
%
% Usage:
%  vol = mmil_detrend_vol(vol,'key1', value1,...); 
%
% Required Parameters:
%   vol: time series (4D) volume with size [nx,ny,nz,nt]
%
% Optional parameters:
%  'skipTRs': number of initial repetitions to remove
%    {default = 0}
%  'norm_flag': [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100)
%    {default = 1}
%  'thresh': when normalize to mean, set voxels with original values
%    less than this to 0
%    {default = 10}
%  'detrend': [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    {default = 2}
%  'regressors': matrix of time series to regress out of vol time series
%    size should be [nt,nr], where nr = number of regressors
%    {default = []}
%  'fname_motion': text file containing 6 columns of motion estimates
%     to be used as regressors (e.g. from AFNI's 3dvolreg)
%    {default = []}
%   ind_valid: vector of time point index numbers
%     to use in calculation of beta coefficients
%
% Created:  08/20/10 Don Hagler
% Prev Mod: 03/23/16 by Don Hagler
% Last Mod: 11/21/17 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'skipTRs',0,[0 Inf],...
  'norm_flag',true,[false true],...
  'normval',100,[1,1e10],...
  'thresh',10,[0 Inf],...
  'detrend',2,[0:2],...
  'regressors',[],[],...
  'fname_motion',[],[],...
  'ind_valid',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get size of volume
[nx,ny,nz,nframes] = size(vol);

%  load fname_motion
if ~isempty(parms.fname_motion)
  parms.motion_data = mmil_load_motion_1D(parms.fname_motion,...
    'skipTRs',parms.skipTRs,'nframes',nframes);
end;

% check regressors
if ~isempty(parms.regressors)
  [nt,nr] = size(parms.regressors);
  if nt ~= nframes
    error('first dim of regressors must match last dim of vol');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% strip dummy TRs
if parms.skipTRs>0 && parms.skipTRs<nframes
  vol = vol(:,:,:,parms.skipTRs+1:end);
  [nx,ny,nz,nframes] = size(vol);
  if ~isempty(parms.regressors)
    parms.regressors = parms.regressors(parms.skipTRs+1:end,:);
  end;
  if ~isempty(parms.ind_valid)
    parms.ind_valid = parms.ind_valid - parms.skipTRs;
    parms.ind_valid = parms.ind_valid(parms.ind_valid>0);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize volume time-series to mean
if parms.norm_flag
  vol = mmil_norm_vol(vol,parms.normval,parms.thresh);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% regress out trend, motion, and other regressors
if parms.detrend || ~isempty(parms.fname_motion) || ~isempty(parms.regressors)
  t = [0:nframes-1]';
  X = ones(size(t));
  if parms.detrend
    for d=1:parms.detrend
      X = [X t.^d];
    end;
  end;
  if ~isempty(parms.fname_motion)
    X = [X parms.motion_data];
  end;
  if ~isempty(parms.regressors)
    X = [X parms.regressors];
  end;
  vol = mmil_regress_vol(vol,X,parms.ind_valid);
end;

