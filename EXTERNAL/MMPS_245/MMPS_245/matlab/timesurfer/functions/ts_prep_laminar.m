function data = ts_prep_laminar(data,varargin)
%function data = ts_prep_laminar(data,[options])
%
% Purpose: prepare laminar data for analysis
%   by replacing bad channels and applying Gaussian smoothing
%
% Required Input:
%   data: 2D data matrix, with size = [nchans,ntpoints]
%
% Optional Parameters:
%   'badchans': vector of bad channels
%     will be replaced with weighted average of surrounding good channels
%     {default = []}
%   'chan_wf': exponential decay factor (0 to 1)
%     used to weight neighbors of bad channels based on relative distance
%     {default = 0.1}
%   'smoothing': smoothing sigma (# of channels)
%     {default = 0}
%
% Output:
%   data: matrix with size = [nchans,ntpoints]
%     with replacement of badchans and smoothing across channels
%     units = mV
%
% Created:  09/01/15 by Don Hagler
% Last Mod: 09/01/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters
[parms,data] = check_input(data,varargin);

% replace bad channels with average of surrounding channels
data = replace_chans(data,parms);

% apply smoothing 1D smoothing across channels
data = smooth_chans(data,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = check_input(data,options)
  parms = mmil_args2parms(options,{...
    'badchans',[],[],...
    'chan_wf',0.1,[0,1],...
    'smoothing',0,[],...
  });
  if ndims(data)~=2
    error('data matrix must be 2D (nchans x ntpoints)');
  end;
  [parms.nchans,parms.ntpoints] = size(data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = replace_chans(data,parms)
  if ~isempty(parms.badchans)
    goodchans = setdiff(1:parms.nchans,parms.badchans);
    T = eye(parms.nchans);
    for i=1:length(parms.badchans)
      badchan = parms.badchans(i);
      t = zeros(1,parms.nchans);
      ind_good = find(goodchans<badchan);
      if ~isempty(ind_good)
        goodchan = goodchans(ind_good(end));
        % calculate weight proportional to distance
        w = calc_chan_weight(goodchan,badchan,parms.nchans,parms.chan_wf);
        t(goodchan) = w;
      end;
      ind_good = find(goodchans>badchan);
      if ~isempty(ind_good)
        goodchan = goodchans(ind_good(1));
        % calculate weight proportional to distance
        w = calc_chan_weight(goodchan,badchan,parms.nchans,parms.chan_wf);
        t(goodchan) = w;
      end;
      % adjust weights to sum to 1
      t = t/sum(t);
      % insert vector into matrix
      T(badchan,:) = t;
    end;
    % apply transformation matrix to data
    data = T*data;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = smooth_chans(data,parms)
  if parms.smoothing>0
    % create Gaussian blurring kernel
    h = fspecial('gaussian',parms.nchans,parms.smoothing);
    % apply kernel to identity matrix
    T = imfilter(eye(parms.nchans),h,'replicate');
    % normalize rows to sum to 1
    T = bsxfun(@rdivide,T,sum(T,2));
    % apply transformation matrix to data
    data = T*data;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = calc_chan_weight(a,b,n,f)
  % distance relative to maximum distance
  d = abs(a-b)/(n-1);
  % exponential decay function
  w = exp(-d/f);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

