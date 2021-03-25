function onset = ts_wform_onset(wforms,varargin)
%function onset = ts_wform_onset(wforms,[options])
%
% Purpose: calculate area under curve
%
% Required Input:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nrois,nconditions]
%
% Optional Parameters:
%   'onset_method': method used for determining threshold
%     allowed: 'quartiles' 'MAD' 'STDEV'
%     {default = 'quartiles'}
%   'onset_minimum': minimum time for onset (msec)
%     {default = 0}
%   'onset_kappa': tunable factor for scaling threshold
%     {default = 2.3}
%   'onset_baseline_flag': [0|1|2] subtract mean or median baseline from wforms
%     0: no baseline subtraction
%     1: subtract mean of each waveform baseline
%     2: subtract median of each waveform baseline
%     {default = 0}
%   'onset_baseline_range': time range (msec) to calculate baseline noise
%     {default = [-100,0]}
%   'onset_baseline_collapse_flag': [0|1|2] whether to collapse noise
%     across area or condition (for calculating threshold)
%     0: separate baseline for each area and condition
%     1: collapse baseline across conditions
%     2: collapse baseline across areas and conditions
%     {default = 0}
%   'onset_polarity': [-1,1] onset of negative or positive deflection
%     {default = 1}
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
%   onset: matrix of size [nrois,nconditions,2] with onset latencies
%     NaN if no onset is found
%
% Created:  02/08/12 by Don Hagler
% Rcnt Mod: 03/02/12 by Don Hagler
% Last Mod: 09/13/12 by Don Hagler
%


%% todo:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nconditions] or
%                              [ntpoints,nrois,nconditions]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'onset_method','quartiles',{'quartiles','MAD','STDEV'},...
  'onset_minimum',0,[],...
  'onset_kappa',2.3,[],...
  'onset_baseline_flag',0,[0 1 2],...
  'onset_baseline_range',[-100,0],[],...
  'onset_baseline_collapse_flag',0,[0 1 2],...
  'onset_polarity',1,[-1,1,1],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
};
onset = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);

if parms.onset_baseline_flag
  wforms = subtr_baseline(wforms,parms);
end;

if parms.onset_baseline_collapse_flag<2
  onset = get_onset(wforms,parms);
else
  onset = get_onset2(wforms,parms);
end;

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
  [tmp,b_ind0] = min(abs(parms.time - parms.onset_baseline_range(1)));
  [tmp,b_ind1] = min(abs(parms.time - parms.onset_baseline_range(2)));
  parms.onset_b_range = [b_ind0:b_ind1];
  [tmp,r_ind0] = min(abs(parms.time - parms.onset_minimum));
  r_ind1 = parms.ntpoints;
  parms.onset_r_range = [r_ind0:r_ind1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wforms = subtr_baseline(wforms,parms)
  nwf = parms.nrois*parms.nconditions;
  responses = reshape(wforms,[parms.ntpoints,nwf]);
  baseline = reshape(wforms(parms.onset_b_range,:,:),...
                      [numel(parms.onset_b_range),nwf]);
  if parms.onset_baseline_flag==1
    avg_baseline = ones(parms.ntpoints,1)*mean(baseline,1);
  elseif parms.onset_baseline_flag==2
    avg_baseline = ones(parms.ntpoints,1)*median(baseline,1);
  else
    avg_baseline = 0;
  end;
  responses = responses - avg_baseline;
  wforms = reshape(responses,[parms.ntpoints,parms.nrois,parms.nconditions]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function onset = get_onset(wforms,parms)
  onset = nan(parms.nrois,parms.nconditions);
  time = parms.time(parms.onset_r_range);
  if parms.onset_baseline_collapse_flag==2
    % find onset thresholds based on quartiles of baseline values
    baseline = mmil_rowvec(wforms(parms.onset_b_range,:,:));
    thresh = calc_thresh(baseline,parms);
  end;
  for r=1:parms.nrois
    if ~parms.onset_baseline_collapse_flag==1
      % find onset thresholds based on quartiles of baseline values
      baseline = mmil_rowvec(wforms(parms.onset_b_range,r,:));
      thresh = calc_thresh(baseline,parms);
    end;
    for c=1:parms.nconditions
      if ~parms.onset_baseline_collapse_flag==0
        % find onset thresholds based on quartiles of baseline values
        baseline = mmil_rowvec(wforms(parms.onset_b_range,r,c));
        thresh = calc_thresh(baseline,parms);
      end;
      response = squeeze(wforms(parms.onset_r_range,r,c));
      if parms.onset_polarity>0
        ind_onset = find(response>thresh,1,'first');
      else
        ind_onset = find(response<thresh,1,'first');
      end;
      if ~isempty(ind_onset)
        onset(r,c) = time(ind_onset);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function onset = get_onset2(wforms,parms)
  onset = [];
  time = parms.time(parms.onset_r_range);
  % find onset thresholds based on quartiles of baseline values
  baseline = mmil_rowvec(wforms(parms.onset_b_range,:,:));
  thresh = calc_thresh(baseline,parms);
  nwf = parms.nrois*parms.nconditions;
  response = reshape(wforms(parms.onset_r_range,:,:),...
    [numel(time),nwf]);
  if parms.onset_polarity>0
    ind_mask = find(response>thresh);
  else
    ind_mask = find(response<thresh);
  end;
  mask = zeros(size(response));
  mask(ind_mask) = 1;
  mask(end-1,:) = 0;
  mask(end,:) = 1; % make sure there is onset for every wf
  mask_diff = cat(1,zeros(1,nwf),diff(mask,1,1));
  ind_change = find(mask_diff~=0);
  if isempty(ind_change), return; end;
  [t_change,n_change] = ind2sub(size(response),ind_change);
  [~,ind_wf] = intersect(flipud(n_change),[1:nwf]);
  ind_wf = 1 + length(n_change) - ind_wf;
  onset = reshape(time(t_change(ind_wf)),[parms.nrois,parms.nconditions]);
  % set to NaN if onset is last time point
  onset(onset>=time(end-1)) = nan;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

function thresh = calc_thresh(baseline,parms)
  switch upper(parms.onset_method)
    case 'QUARTILES'
      q = prctile(baseline,[25,50,75]);
      m = q(2);
      q = q(3) - q(1);
    case 'MAD'
      % calculate median absolute deviation
      m = median(baseline);
      q = max(median(abs(baseline - m)) / 0.6745,eps);
    case 'STDEV'
      m = mean(baseline);
      q = std(baseline);
  end;
  thresh = m + parms.onset_polarity*parms.onset_kappa*q;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

