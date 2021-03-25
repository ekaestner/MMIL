function results = ts_wform_peaks(wforms,varargin)
%function results = ts_wform_peaks(wforms,[options])
%
% Purpose: detect peaks and summarize measures for source estimate waveforms
%   including peak amplitude and latency
%
% Required Input:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nrois,nconditions]
%
% Optional Parameters:
%   'smooth_sigma': Gaussian blurring kernel applied to waveforms
%     {default = 0}
%   'peak_range': time range (msec) within which to select peaks
%     {default = [50,150]}
%   'peak_pol': peak polarity (-1 or 1)
%     {default = 1}
%   'peak_mindiff': minimum difference for peak detection
%     {default = 0.5}
%   'firstflag': use first peak in time range, otherwise use largest amplitude
%     {default = 0}
%   'contrast_latency_flag': [0|1] require longer latency for 
%     lower condition values (e.g. for low luminance contrast)
%     {default = 0}
%   'latency_allowance': number of milliseconds to allow lower condition values
%      to have earlier latency (if contrast_latency_flag = 1)
%     {default = 5}
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
%   'resfact': temporal ressampling factor
%     use values greater than 1 (e.g. 10) to avoid step-like
%       discretization of latency
%     {default = 1}
%
% Output:
%   results: struct containing fields:
%     amplitude: peak amplitude; [nrois, nconditions]
%     deflection: peak-to-avg-peak deflection; [nrois,nconditions]
%     latency: peak latency; [nrois, nconditions]
%     nrois
%     nconditions
%
% Created:  02/08/12 by Don Hagler
% Last Mod: 04/07/14 by Don Hagler
%

%% todo: for selected peak, go back in time and find onset

%% todo:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nconditions] or
%                              [ntpoints,nrois,nconditions]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'smooth_sigma',0,[0,1000],...
  'peak_range',[50,150],[],...
  'peak_pol',1,[-1,1,1],...
  'peak_mindiff',0.5,[0,Inf],...
  'firstflag',false,[false true],...
  'contrast_latency_flag',false,[false true],...
  'latency_allowance',5,[0,100],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
  'resfact',1,[1,1000],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);
results = init_results(parms);
if parms.smooth_sigma
  wforms = smooth_waveforms(wforms,parms);
end;
results = get_peaks(wforms,parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nrois = size(wforms,2);
  parms.nconditions = size(wforms,3);
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  if parms.resfact~=1
    parms.t0 = parms.time(1);
    parms.t1 = parms.time(end);
    sd = parms.time(2) - parms.time(1);
    sd_res = sd/parms.resfact;
    parms.time_res = [parms.t0:sd_res:parms.t1];
  else
    parms.time_res = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  fnames = {'nrois','nconditions','time','time_res'};
  for f=1:length(fnames)
    results.(fnames{f}) = parms.(fnames{f});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wforms = smooth_waveforms(wforms,parms)
  for r=1:parms.nrois
    for c=1:parms.nconditions
      wforms(:,r,c) = mmil_smooth(squeeze(wforms(:,r,c)),parms.smooth_sigma);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = get_peaks(wforms,parms,results)
  % find minima and maxima
  clear minima maxima;
  for c=1:results.nconditions
    for r=1:results.nrois
      % select time range of wform
      tmp_wform = wforms(:,r,c);
      if parms.resfact~=1
        tmp_wform = spline(parms.time,tmp_wform,parms.time_res);
        time = parms.time_res;
      else
        time = parms.time;
      end;
      [maxtab,mintab] = mmil_peakdet(tmp_wform,parms.peak_mindiff,time);
      tmp_minima = get_minmax(parms,mintab);
      tmp_maxima = get_minmax(parms,maxtab);
      tmp_minima = get_deflection(tmp_minima,tmp_maxima);
      tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
      minima(r,c) = tmp_minima;
      maxima(r,c) = tmp_maxima;
    end;
  end;
  results.minima = minima;
  results.maxima = maxima;

  % get peak amplitude, deflection, and latency
  amplitude = nan(results.nrois,results.nconditions);
  deflection = amplitude;
  latency = amplitude;
  % peaks
  if parms.peak_pol>0
    extrema = results.maxima;
  else
    extrema = results.minima;
  end;

  cond_order = [1:results.nconditions];
  if parms.contrast_latency_flag
    cond_order = fliplr(cond_order);
  end;

  for cond=1:results.nconditions
    c = cond_order(cond);
    for r=1:results.nrois
      % peak amplitude, deflection, and latency
      peaks = extrema(r,c); % for one waveform (roi x contrast)
      npeaks = length(peaks.latency);
      % select one peak
      if npeaks
        ind = 1:npeaks;
        % require that latency increases with decreasing condition value
        %   NOTE: this follows the parvo pathway
        if parms.contrast_latency_flag && cond>1
          min_latency = latency(r,cond_order(cond-1)) - parms.latency_allowance;
          if ~isnan(min_latency)
            ind = find(peaks.latency>=min_latency);
          end;
        end;
        npeaks = length(ind);
      end;
      if ~npeaks
%        fprintf('%s: missing peak for condition %d, area %d\n',...
%          mfilename,c,a);
        continue;
      end;
      % select first or largest peak
      if parms.firstflag || npeaks==1
        ind = ind(1);
      else
        [tmp,tmp_ind] = max(parms.peak_pol*peaks.amplitude(ind));
        ind = ind(tmp_ind);
      end;
      deflection(r,c) = peaks.deflection(ind);
      amplitude(r,c) = peaks.amplitude(ind);
      latency(r,c) = peaks.latency(ind);
    end;
  end;
  results.amplitude = amplitude;
  results.deflection = deflection;
  results.latency = latency;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax(parms,mtab)
  minmax = struct('latency',[],'amplitude',[],'deflection',[]);
  j = 1;
  for i=1:size(mtab,1)
    latency = mtab(i,1);
    amplitude = mtab(i,2);
    if latency<parms.peak_range(1) ||...
       latency>parms.peak_range(2)
      continue;
    end;
    minmax.latency(j) = latency;
    minmax.amplitude(j) = amplitude;
    j = j + 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax1 = get_deflection(minmax1,minmax2)
  % calculate deflection for minima or maxima
  for i=1:length(minmax1.latency)
    latency = minmax1.latency(i);
    % find preceding minimum
    tmp = latency - minmax2.latency;
    tmp(tmp<0)=Inf;
    [tmp,ind1] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp1 = 0;
    else
      amp1 = minmax2.amplitude(ind1);
    end;

    % find subsequent maxima or minima
    tmp = minmax2.latency - latency;
    tmp(tmp<0)=Inf;
    [tmp,ind2] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp2 = 0;
    else
      amp2 = minmax2.amplitude(ind2);
    end;

    minmax1.deflection(i) =...
      abs(minmax1.amplitude(i) - mean([amp1 amp2]));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

