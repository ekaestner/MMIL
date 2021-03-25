function wforms = ts_baseline_wforms(wforms,varargin)
%function wforms = ts_baseline_wforms(wforms,[options])
%
% Purpose: subtract baseline from waveform matrix
%
% Required Input:
%   wforms: 2D waveform matrix
%     size of wforms should be [ntpoints,nconds]
%
% Optional Parameters:
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 300}
%   'offset': subtract this value from waveforms
%     {default = 0}
%   'baseline_flag': [0|1|2] subtract mean of baseline from waveforms
%     0: no baseline subtraction
%     1: subtract baseline of each condition
%     2: subtract baseline averaged across all conditions
%     {default = 1}
%   'baseline_t0': start time of baseline period
%     {default = -Inf}
%   'baseline_t1': end time of baseline period
%     {default = 0}
%
% Created:  03/21/12 by Don Hagler
% Last Mod: 03/21/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'time',[],[],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'xlim',[],[],...
  'ylim',[],[],...
  'offset',0,[-Inf,Inf],...
  'baseline_flag',1,[0 1 2],...
  'baseline_t0',-Inf,[],...
  'baseline_t1',0,[],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);

if parms.baseline_flag || parms.offset
  wforms = baseline_waveforms(wforms,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nconds = size(wforms,2);
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  if parms.baseline_flag
    [tmp,parms.baseline_s0] = min(abs(parms.time - parms.baseline_t0));
    [tmp,parms.baseline_s1] = min(abs(parms.time - parms.baseline_t1));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wforms,parms] = baseline_waveforms(wforms,parms)
  if parms.baseline_flag == 2
    tmp_baseline = wforms(parms.baseline_s0:parms.baseline_s1,:);
    tmp_baseline = mean(tmp_baseline(:));
  end;
  for c=1:parms.nconds
    tmp_wform = wforms(:,c);
    if parms.baseline_flag == 1
      tmp_baseline = mean(tmp_wform(parms.baseline_s0:parms.baseline_s1));
    elseif parms.baseline_flag == 0;
      tmp_baseline = 0;
    end;
    if parms.offset
      tmp_baseline = tmp_baseline + parms.offset;
    end;
    % adjust wforms
    wforms(:,c) = tmp_wform - tmp_baseline;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
