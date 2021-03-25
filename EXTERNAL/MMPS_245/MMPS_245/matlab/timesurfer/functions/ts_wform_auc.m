function auc = ts_wform_auc(wforms,varargin)
%function auc = ts_wform_auc(wforms,[options])
%
% Purpose: calculate area under curve
%
% Required Input:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nareas,nconditions]
%
% Optional Parameters:
%   'auc_range': time range (msec) within which to calculate area under curve
%     {default = [0,350]}
%   'auc_nbins': number of bins in which to subdivide auc_range
%     {default = 1}
%   'auc_baseline_flag': whether to subtract baseline area under curve
%     {default = 1}
%   'auc_baseline_range': time range (msec) to calculate baseline area under curve
%     {default = [-100,0]}
%   'powerflag': whether to use magnitude or power
%     0: raw values      (signed magnitude)
%     1: absolute values (unsigned magnitude)
%     2: squared values  (power)
%     {default = 2}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 300}
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%
% Output:
%   auc: matrix of size [nareas,nconditions] with area under curve values
%     if auc_nbins > 1, auc will have size [nareas,nconditions,nbins+1]
%     extra bin is for sum across all bins (for convenience)
%
% Created:  02/08/12 by Don Hagler
% Rcnt Mod: 03/22/12 by Don Hagler
% Last Mod: 09/13/12 by Don Hagler
%

%% todo:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nconditions] or
%                              [ntpoints,nareas,nconditions]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'auc_range',[0,350],[],...
  'auc_nbins',1,[1,100],...
  'auc_baseline_flag',true,[false true],...
  'auc_baseline_range',[-100,0],[],...
  'powerflag',2,[0 1 2],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
};
auc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);
auc = get_auc(wforms,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nareas = size(wforms,2);
  parms.nconditions = size(wforms,3);
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  if parms.auc_nbins>1
    parms.auc_bins = zeros(parms.auc_nbins+1,2); % extra bin for entire range
    range_dur = range(parms.auc_range);
    bin_dur = range_dur / parms.auc_nbins;
    for n=1:parms.auc_nbins
      parms.auc_bins(n,1) = parms.auc_range(1) + (n-1)*bin_dur;
      parms.auc_bins(n,2) = parms.auc_bins(n,1) + bin_dur;
    end;
    parms.auc_nbins = parms.auc_nbins + 1; % extra bin for entire range
    parms.auc_bins(parms.auc_nbins,:) = parms.auc_range;
  else
    parms.auc_bins = parms.auc_range;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auc = get_auc(wforms,parms)
  auc = nan(parms.nareas,parms.nconditions,parms.auc_nbins);
  switch parms.powerflag
    case 1
      wforms = abs(wforms);
    case 2
      wforms = wforms.^2;
  end;
  for n=1:parms.auc_nbins
    [tmp,ind0] = min(abs(parms.time - parms.auc_bins(n,1)));
    [tmp,ind1] = min(abs(parms.time - parms.auc_bins(n,2)));
    auc_range = [ind0:ind1];
    n_auc_range = length(auc_range);
    if parms.auc_baseline_flag
      [tmp,b_ind0] = min(abs(parms.time - parms.auc_baseline_range(1)));
      [tmp,b_ind1] = min(abs(parms.time - parms.auc_baseline_range(2)));
      auc_b_range = [b_ind0:b_ind1];
    end;
    j = 0;
    for c=1:parms.nconditions
      for a=1:parms.nareas
        % area under curve (integrate over time range)
        if parms.auc_baseline_flag
          % average baseline across conditions
          auc_b = mean(mmil_rowvec(wforms(auc_b_range,a,:)));
        else
          auc_b = 0;
        end;
        auc(a,c,n) = mean(squeeze(wforms(auc_range,a,c))) - auc_b;
      end;
    end;
  end;
return;


