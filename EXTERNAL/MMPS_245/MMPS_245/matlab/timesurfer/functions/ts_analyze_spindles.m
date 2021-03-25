function results = ts_analyze_spindles(spindles,varargin)
%function results = ts_analyze_spindles(spindles,[options])
%
% Usage:
%  metrics = ts_analyze_spindles(spindles,'key1', value1,...);
%
% Required Input:
%   spindles: struct containing results from ts_find_spindles
%
% Optional parameters
%  'coher_flag': calculate cross-channel coherence in spindle frequency range
%    NOTE: this is slow!
%    {default = 1}
%  'spindle_freq_low': lower bound of spindle frequency range (in Hz)
%    {default = 9}
%  'spindle_freq_high': upper bound of spindle frequency range (in Hz)
%    {default = 17}
%  'corr_nchans_max': maximum number of channels to use for cross-channel
%     correlation calculation
%    {default = 500}
%  'recalc_prob_flag': [0|1] re-calculate spindle occurrence probabilities
%     for each channel based on amplitude and duration
%     NOTE: does not take into account the additional rejection criteria
%       used by ts_find_spindles (i.e. fourier_thresh, min_npeaks, min_peakrate)
%    {default = 0}
%  'peak_thresh': threshold for spindle peak detection
%    {default = 2}
%  'edge_thresh': threshold for spindle edge detection
%    results included in timing_fixed field
%    {default = 1}
%  'relative_edge_thresh_flag': [0|1] determine edge threshold
%     based on a fraction of peak amplitude
%     this option determines method used for rejection based on spindle duration
%    {default = 1}
%  'relative_edge_thresh': fraction of peak amplitude used for edge detection
%    results included in timing_relative field
%    {default = 0.2}
%  'reject_dur_flag': [0|1] reject spindles with short duration
%    for spindle probability metrics
%    {default = 1}
%  'dur_thresh': threshold for the minimal duration of the spindle
%    epoch (in seconds) in a given channel
%    {default = 0.2}
%  'min_epoch_delay': threshold for the minimum delay between epochs
%    epoch pairs with delays shorter than this will be combined
%    for purposes of calculating epoch delays and epoch durations
%    {default = 0.1}
%  'mask': mask vector with size = [1,ntpoints]
%     with values of 0 or 1 to exclude artifacts
%     used in analysis to exclude delays between spindles containing artifacts
%     {default = []}
%  'fft_norm_flag': [0|1] normalize FFT power spectra by frequency
%    {default = 1}
%  'fft_norm_exp': exponent applied to frequency when normalizing FFT spectra
%    {default = 1.2}
%  'verbose': [0|1] display status messages
%    {default = 0}
%
% Output:
%   results: struct containing various spindle metrics
%
% Created:  02/12/14 by Don Hagler
% Last mod: 03/03/16 by Don Hagler
%

% initially created as ts_calc_spindle_metrics

%% todo: calculate timing (dur, delay, interval) for each channel

%% todo: split into separate functions (files) for each metric type?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];

% parse input parameters
parms = check_input(spindles,varargin);

% initialize results struct
results = init_results(parms,spindles);

% calculate spindle co-occurrence probability across channels
if parms.verbose
  fprintf('%s: calculating spindle probabilities...\n',mfilename);
end;
results = calc_prob(spindles,parms,results,0);
if parms.recalc_prob_flag
  results.prob_orig = results.prob;
  results = calc_prob(spindles,parms,results,1);
end;

% calculate spindle amplitudes
if parms.verbose
  fprintf('%s: calculating spindle amplitudes...\n',mfilename);
end;
results = calc_amp(spindles,parms,results);

% calculate spindle onset and duration
if parms.verbose
  fprintf('%s: calculating spindle timing...\n',mfilename);
end;
% relative edge threshold
results = calc_timing(spindles,parms,results,1);
% fixed edge threshold
results = calc_timing(spindles,parms,results,0);

% calculate correlations between channels
if parms.nchans>1
  if parms.verbose
    fprintf('%s: calculating cross-channel correlations...\n',mfilename);
  end;
  results = calc_corr(spindles,parms,results);

  % calculate coherence
  if parms.coher_flag
    if parms.verbose
      fprintf('%s: calculating cross-channel coherence...\n',mfilename);
    end;
    results = calc_coher(spindles,parms,results);
  end;
end;

% calculate mean and coefficient of variation of RMS
if parms.verbose
  fprintf('%s: calculating RMS amplitudes...\n',mfilename);
end;
results = calc_rms(spindles,parms,results);

% calculate amplitude and phase from FFT
if parms.verbose
  fprintf('%s: calculating FFT amplitudes and phases...\n',mfilename);
end;
if parms.nchans>1
  results = calc_fft(spindles,parms,results,1);
end;
results = calc_fft(spindles,parms,results,0);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(spindles,options)
  parms = mmil_args2parms(options,{...
    'coher_flag',true,[false true],...
    'spindle_freq_low',9,[0.5,100],...
    'spindle_freq_high',17,[0.5,100],...
    'corr_nchans_max',500,[],...
    'recalc_prob_flag',false,[false true],...
    'peak_thresh',2,[0,Inf],...
    'edge_thresh',1,[0,Inf],...
    'relative_edge_thresh_flag',true,[false true],...
    'relative_edge_thresh',0.2,[0,1],...
    'reject_dur_flag',true,[false true],...
    'dur_thresh',0.2,[0,Inf],...
    'min_epoch_delay',0.1,[0,100],...
    'mask',[],[],...
    'fft_norm_flag',true,[false true],...
    'fft_norm_exp',1.2,[0,10],...
    'verbose',false,[false true],...
  ...
    'required_fields',{'sfreq','ind_spindle_rel','spindle_data'},[],...
    'fft_freq',[0:0.1:30],[],...
  });

  % check that input struct has required fields
  fnames = fieldnames(spindles);
  missing_fields = setdiff(parms.required_fields,fnames);
  if ~isempty(missing_fields)
    error('these required fields are missing from input spindles struct:\n%s\n',...
      sprintf('%s ',missing_fields{:}));   
  end;
  
  % extract some information for later use
  parms.sfreq = spindles.sfreq;
  parms.nepochs = size(spindles.spindle_data,1);
  parms.nchans = size(spindles.spindle_data,2);
  parms.ntpoints = size(spindles.spindle_data,3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms,spindles);
  results = [];
  fnames = {'sfreq','nepochs','nchans','ntpoints'};
  for i=1:length(fnames)
    results.(fnames{i}) = parms.(fnames{i});
  end;
  results.labels = spindles.labels;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_prob(spindles,parms,results,recalc_flag)
  if ~exist('recalc_flag','var'), recalc_flag = 0; end;
  spindle_chans = find_spindle_chans(spindles,parms,recalc_flag);
  results = update_results_prob(results,spindle_chans,recalc_flag,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle_chans = find_spindle_chans(spindles,parms,recalc_flag)
  if ~recalc_flag
    spindle_chans = spindles.spindle_chans;
  else
    spindle_chans = cell(parms.nepochs,1);
    for i=1:parms.nepochs
      % extract spindle power wave data between start and end of spindle
      ind1 = spindles.ind_spindle_rel(i,1);
      ind2 = max(spindles.ind_spindle_rel(i,2),ind1+1);
      data = squeeze(spindles.spindle_power_peak(i,:,ind1:ind2));
      peak_flags = zeros(parms.nchans,1);
      for j=1:parms.nchans
        chan_data = data(j,:);
        nt = length(chan_data);
        % find peaks
        [peak_vals,ind_peaks,~,~] = mmil_extrema(chan_data);
        % exclude edges
        ind_non_edge = find(ind_peaks>1 & ind_peaks<nt);
        if isempty(ind_non_edge), continue; end;
        % test whether peaks are large enough
        if any(peak_vals >= parms.peak_thresh)
          peak_flags(j) = 1;
        end;
      end;
      chans = find(peak_flags);
      spindle_chans{i} = chans';
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle_counts = calc_spindle_counts(spindle_chans,parms)
  spindle_counts = zeros(parms.nchans,parms.nchans);
  for i=1:parms.nepochs
    chans = spindle_chans{i};
    for j=1:length(chans)
      c1 = chans(j);
      spindle_counts(c1,c1) = spindle_counts(c1,c1) + 1;
      chans2 = setdiff(chans,c1);
      for k=1:length(chans2)
        c2 = chans2(k);
        spindle_counts(c1,c2) = spindle_counts(c1,c2) + 1;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = update_results_prob(results,spindle_chans,recalc_flag,parms)
  fname = 'prob';
  spindle_counts = calc_spindle_counts(spindle_chans,parms);
  dvals = diag(spindle_counts);
  spindle_freq = spindle_counts/parms.nepochs;
  spindle_freq_norm = bsxfun(@rdivide,spindle_counts,dvals);
  spindle_chans_mask = zeros(parms.nepochs,results.nchans);;
  for i=1:parms.nepochs
    chans = spindle_chans{i};
    spindle_chans_mask(i,chans) = 1;
  end;
  nchans = sum(spindle_chans_mask,2);
  
  results.(fname) = [];
  results.(fname).spindle_chans = spindle_chans;
  results.(fname).spindle_counts = spindle_counts;
  results.(fname).spindle_freq = spindle_freq;
  results.(fname).spindle_freq_norm = spindle_freq_norm;
  results.(fname).spindle_chans_mask = spindle_chans_mask;
  results.(fname).nchans = nchans;
  results.(fname).nchans_norm = nchans / parms.nchans;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_amp(spindles,parms,results)
  fname = 'amp';
  chans_mask = results.prob.spindle_chans_mask;
  % copy amplitudes from pre-calculated spindle_metrics
  amp = spindles.spindle_metrics.amp;
  amp_chan = nan(parms.nepochs,results.nchans);
  amp_chan(chans_mask>0) = amp(chans_mask>0);
  results.(fname) = [];
  results.(fname).amp = amp;
  results.(fname).amp_chan = amp_chan;
  % calculate mean and std of amp
  results.(fname).amp_mean = ts_nan_mean(amp,1);
  results.(fname).amp_std = ts_nan_std(amp,0,1);
  results.(fname).amp_chan_mean = ts_nan_mean(amp_chan,1);
  results.(fname).amp_chan_std = ts_nan_std(amp_chan,0,1);
  % calculate cross-channel correlation of amplitude
  if parms.nepochs
    results.(fname).amp_corr = corr(amp);
  end;
  % calculate maximum amplitude across channels
  results.(fname).amp_max = max(amp,[],2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_timing(spindles,parms,results,relative_edge_thresh_flag)
  if ~exist('relative_edge_thresh_flag','var')
    relative_edge_thresh_flag = 1;
  end;
  if relative_edge_thresh_flag
    fname = 'timing_relative';
  else
    fname = 'timing_fixed';
  end;
  results.(fname) = [];
  spindle_onset = nan(parms.nepochs,parms.nchans);
  spindle_offset = nan(parms.nepochs,parms.nchans);
  spindle_dur = nan(parms.nepochs,parms.nchans);
  spindle_peak = nan(parms.nepochs,parms.nchans);
  spindle_onset_std_chan = nan(parms.nepochs,1);
  spindle_dur_std_chan = nan(parms.nepochs,1);
  spindle_peak_std_chan = nan(parms.nepochs,1);
  for i=1:parms.nepochs
    % extract spindle power wave data between start and end of spindle
    ind1 = spindles.ind_spindle_rel(i,1);
    ind2 = max(spindles.ind_spindle_rel(i,2),ind1+1);
    data = spindles.spindle_power_edge(i,:,ind1:ind2);
    data = reshape(data,[size(data,2),size(data,3)]);
    chans = results.prob.spindle_chans{i};
    for j=1:length(chans)
      c = chans(j);
      chan_data = data(c,:);
      nt = length(chan_data);
      [tmp,ind_peak] = max(chan_data);
      if ind_peak==1 || ind_peak==nt
        % exclude peaks at edge of epoch if possible
        [peak_vals,ind_peaks,~,~] = mmil_extrema(chan_data);
        if length(ind_peaks)>1
          ind_non_edge = find(ind_peaks>1 & ind_peaks<nt);
          if ~isempty(ind_non_edge)
            peak_vals = peak_vals(ind_non_edge);
            ind_peaks = ind_peaks(ind_non_edge);
          end;
        end;
        % use largest remaining peak (already sorted by size)
        ind_peak = ind_peaks(1);
      end;
      % use relative or fixed threshold
      if relative_edge_thresh_flag
        edge_thresh = parms.relative_edge_thresh * chan_data(ind_peak);
      else
        edge_thresh = parms.edge_thresh;
      end;
      % find time points below threshold
      tmp_thresh = chan_data >= edge_thresh;
      tmp_diff = diff([0,tmp_thresh]);
      % find earliest crossing above threshold
      ind_cross = find(tmp_diff > 0);
      ind_cross = ind_cross(ind_cross<=ind_peak);
      if length(ind_cross>1)
        % exclude edge if possible
        ind_first = min(ind_cross(ind_cross>1));
      else
        ind_first = ind_cross;
      end;
      % find latest crossing below threshold
      ind_cross = find(tmp_diff < 0);
      ind_cross = ind_cross(ind_cross>=ind_peak);
      if length(ind_cross>1)
        % exclude edge if possible
        ind_last = max(ind_cross(ind_cross<nt));
      else
        ind_last = ind_cross;
      end;
      if isempty(ind_first), ind_first = 1; end;
      if isempty(ind_last), ind_last = nt; end;
      % calculate onset and duration
      spindle_onset(i,c) = ind_first/parms.sfreq;
      spindle_offset(i,c) = ind_last/parms.sfreq;
      spindle_dur(i,c) = (ind_last - ind_first + 1)/parms.sfreq;
      % calculate peak latency
      spindle_peak(i,c) = ind_peak/parms.sfreq;
    end;
    % calculate variation in spindle onsets between channels
    if length(chans)>1
      spindle_onset_std_chan(i) = ts_nan_std(spindle_onset(i,:));
      spindle_dur_std_chan(i) = ts_nan_std(spindle_dur(i,:));
      spindle_peak_std_chan(i) = ts_nan_std(spindle_peak(i,:));
    end;
  end;

  % reject spindles that are too short in duration
  if parms.reject_dur_flag &&...
       relative_edge_thresh_flag == parms.relative_edge_thresh_flag
    ind_reject = find(spindle_dur < parms.dur_thresh);
    if ~isempty(ind_reject)
      spindle_chans = results.prob.spindle_chans;
      spindle_onset(ind_reject) = nan;
      spindle_offset(ind_reject) = nan;
      spindle_dur(ind_reject) = nan;
      spindle_peak(ind_reject) = nan;
      % recalculate variation in spindle onset, dur, and peak between channels
      for i=1:parms.nepochs
        chans = find(~isnan(spindle_dur(i,:)));
        if length(chans)>1
          spindle_onset_std_chan(i) = ts_nan_std(spindle_onset(i,:));
          spindle_dur_std_chan(i) = ts_nan_std(spindle_dur(i,:));
          spindle_peak_std_chan(i) = ts_nan_std(spindle_peak(i,:));
        end;
        spindle_chans{i} = chans;
      end;
      results = update_results_prob(results,spindle_chans,0,parms);
    end;
  end;

  % mean across channels
  spindle_onset_mean_chan = ts_nan_mean(spindle_onset,2);
  spindle_dur_mean_chan = ts_nan_mean(spindle_dur,2);
  spindle_peak_mean_chan = ts_nan_mean(spindle_peak,2);
  % mean across epochs
  spindle_onset_mean_epoch = ts_nan_mean(spindle_onset,1);
  spindle_dur_mean_epoch = ts_nan_mean(spindle_dur,1);
  spindle_peak_mean_epoch = ts_nan_mean(spindle_peak,1);
  % variation across epochs
  spindle_onset_std_epoch = ts_nan_std(spindle_onset,0,1);
  spindle_dur_std_epoch = ts_nan_std(spindle_dur,0,1);
  spindle_peak_std_epoch = ts_nan_std(spindle_peak,0,1);
  % mean of variation across channels
  spindle_onset_std_chan_mean = ts_nan_mean(spindle_onset_std_chan);
  spindle_dur_std_chan_mean = ts_nan_mean(spindle_dur_std_chan);
  spindle_peak_std_chan_mean = ts_nan_mean(spindle_peak_std_chan);
  % mean of variation across epochs
  spindle_onset_std_epoch_mean = ts_nan_mean(spindle_onset_std_epoch);
  spindle_dur_std_epoch_mean = ts_nan_mean(spindle_dur_std_epoch);
  spindle_peak_std_epoch_mean = ts_nan_mean(spindle_peak_std_epoch);
  % mean across epochs and channels
  spindle_onset_mean = ts_nan_mean(spindle_onset(:));
  spindle_dur_mean = ts_nan_mean(spindle_dur(:));
  spindle_peak_mean = ts_nan_mean(spindle_peak(:));
  % variation across epochs and channels
  spindle_onset_std = ts_nan_std(spindle_onset(:));
  spindle_dur_std = ts_nan_std(spindle_dur(:));
  spindle_peak_std = ts_nan_std(spindle_peak(:));

  % calculate onset differences between channels
  spindle_onset_mean_diff = zeros(results.nchans,results.nchans);
  spindle_onset_std_diff = zeros(results.nchans,results.nchans);
  spindle_onset_count_diff = zeros(results.nchans,results.nchans);
  for j=1:results.nchans
    for k=1:results.nchans
      tmp_diff = spindle_onset(:,j) - spindle_onset(:,k);
      spindle_onset_mean_diff(j,k) = ts_nan_mean(tmp_diff,1);
      spindle_onset_std_diff(j,k) = ts_nan_std(tmp_diff,1);
      spindle_onset_count_diff(j,k) = length(find(~isnan(tmp_diff)));
    end;
  end;
  spindle_onset_zscore_diff = spindle_onset_mean_diff ./...
    (spindle_onset_std_diff./sqrt(spindle_onset_count_diff));
  spindle_onset_zscore_diff(isnan(spindle_onset_zscore_diff)) = 0;

  % calculate peak differences between channels
  spindle_peak_mean_diff = zeros(results.nchans,results.nchans);
  spindle_peak_std_diff = zeros(results.nchans,results.nchans);
  spindle_peak_count_diff = zeros(results.nchans,results.nchans);
  for j=1:results.nchans
    for k=1:results.nchans
      tmp_diff = spindle_peak(:,j) - spindle_peak(:,k);
      spindle_peak_mean_diff(j,k) = ts_nan_mean(tmp_diff,1);
      spindle_peak_std_diff(j,k) = ts_nan_std(tmp_diff,1);
      spindle_peak_count_diff(j,k) = length(find(~isnan(tmp_diff)));
    end;
  end;
  spindle_peak_zscore_diff = spindle_peak_mean_diff ./...
    (spindle_peak_std_diff./sqrt(spindle_peak_count_diff));
  spindle_peak_zscore_diff(isnan(spindle_peak_zscore_diff)) = 0;

  % calculate epoch durations and delays
  epoch_onset = spindles.ind_spindle(:,1)/spindles.sfreq;
  epoch_onset = epoch_onset + min(spindle_onset,[],2);
  epoch_offset = epoch_onset + max(spindle_offset,[],2);
  epoch_dur = epoch_offset - epoch_onset;
  epoch_delay = epoch_onset(2:end) - epoch_offset(1:end-1);
  nepochs = length(epoch_dur);
  ind_epoch_dur = [1:nepochs];
  ind_epoch_delay = [2:nepochs];

  % combine across epochs if delay is too short
  if parms.min_epoch_delay>0
    ind_short = find(epoch_delay < parms.min_epoch_delay);
    for j=1:length(ind_short)
      k = ind_short(j);
      if k==1, continue; end;
      epoch_dur(k-1) = epoch_dur(k-1) + epoch_delay(k) + epoch_dur(k);
    end;
    epoch_dur(ind_short) = [];
    ind_epoch_dur(ind_short) = [];
    epoch_delay(ind_short) = [];
    ind_epoch_delay(ind_short) = [];
  end;

  % calculate intervals between epoch centers
  epoch_interval = epoch_delay + ...
                   epoch_dur(1:end-1)/2 + ...
                   epoch_dur(2:end)/2;
  ind_epoch_interval = ind_epoch_delay;                   

  % exclude NaN's
  ind_nan = find(isnan(epoch_dur));
  if ~isempty(ind_nan)
    epoch_dur(ind_nan) = [];
    ind_epoch_dur(ind_nan) = [];
  end;
  ind_nan = find(isnan(epoch_delay) | isnan(epoch_interval));
  if ~isempty(ind_nan)
    epoch_delay(ind_nan) = [];
    ind_epoch_delay(ind_nan) = [];
    epoch_interval(ind_nan) = [];
    ind_epoch_interval(ind_nan) = [];
  end;  

  % exclude delays (and intervals) if they contain masked time points
  if ~isempty(parms.mask)
    ntmp = length(ind_epoch_delay);
    mask_flags = zeros(ntmp,1);
    for i=1:length(ind_epoch_delay)
      j = ind_epoch_delay(i);
      delay_onset = epoch_offset(j-1);
      delay_offset = epoch_onset(j);
      s0 = round(delay_onset * parms.sfreq);
      s1 = round(delay_offset * parms.sfreq);
      % check if mask includes these time points
      if any(parms.mask(s0:s1)==0)
        mask_flags(i) = 1;
      end;
    end;
    if any(mask_flags)
      fprintf('%s: excluding %d spindle delays because of mask\n',...
        mfilename,length(find(mask_flags==1)));
      ind_mask = find(mask_flags);
      ind_epoch_delay(ind_mask) = [];
      epoch_delay(ind_mask) = [];
      ind_epoch_interval(ind_mask) = [];
      epoch_interval(ind_mask) = [];
    end;
  end;

  % save onset results
  results.(fname).spindle_onset = spindle_onset;
  results.(fname).spindle_onset_mean = spindle_onset_mean;
  results.(fname).spindle_onset_std = spindle_onset_std;
  results.(fname).spindle_onset_mean_epoch = spindle_onset_mean_epoch;
  results.(fname).spindle_onset_std_epoch = spindle_onset_std_epoch;
  results.(fname).spindle_onset_std_epoch_mean = spindle_onset_std_epoch_mean;
  results.(fname).spindle_onset_std_chan_mean = spindle_onset_std_chan_mean;
  results.(fname).spindle_onset_mean_diff = spindle_onset_mean_diff;
  results.(fname).spindle_onset_std_diff = spindle_onset_std_diff;
  results.(fname).spindle_onset_count_diff = spindle_onset_count_diff;
  results.(fname).spindle_onset_zscore_diff = spindle_onset_zscore_diff;

  % save duration results
  results.(fname).spindle_dur = spindle_dur;
  results.(fname).spindle_dur_mean = spindle_dur_mean;
  results.(fname).spindle_dur_std = spindle_dur_std;
  results.(fname).spindle_dur_mean_epoch = spindle_dur_mean_epoch;
  results.(fname).spindle_dur_std_epoch = spindle_dur_std_epoch;
  results.(fname).spindle_dur_std_epoch_mean = spindle_dur_std_epoch_mean;
  results.(fname).spindle_dur_std_chan_mean = spindle_dur_std_chan_mean;

  % save peak results
  results.(fname).spindle_peak = spindle_peak;
  results.(fname).spindle_peak_mean = spindle_peak_mean;
  results.(fname).spindle_peak_std = spindle_peak_std;
  results.(fname).spindle_peak_mean_epoch = spindle_peak_mean_epoch;
  results.(fname).spindle_peak_std_epoch = spindle_peak_std_epoch;
  results.(fname).spindle_peak_std_epoch_mean = spindle_peak_std_epoch_mean;
  results.(fname).spindle_peak_std_chan_mean = spindle_peak_std_chan_mean;
  results.(fname).spindle_peak_mean_diff = spindle_peak_mean_diff;
  results.(fname).spindle_peak_std_diff = spindle_peak_std_diff;
  results.(fname).spindle_peak_count_diff = spindle_peak_count_diff;
  results.(fname).spindle_peak_zscore_diff = spindle_peak_zscore_diff;
  
  % save epoch dur and delay results
  results.(fname).ind_epoch_dur = ind_epoch_dur;
  results.(fname).epoch_dur = epoch_dur;
  results.(fname).ind_epoch_delay = ind_epoch_delay;
  results.(fname).epoch_delay = epoch_delay;
  results.(fname).ind_epoch_interval = ind_epoch_interval;
  results.(fname).epoch_interval = epoch_interval;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_corr(spindles,parms,results)
  % initialize output values
  results.corr = [];
  results.corr.epochs.R = zeros(parms.nepochs,parms.nchans,parms.nchans);
  results.corr.epochs.R2 = zeros(parms.nepochs,parms.nchans,parms.nchans);
  results.corr.epochs.r = zeros(parms.nepochs,1);
  results.corr.epochs.r2 = zeros(parms.nepochs,1);
  % limit number of channels used to calculate correlations  
  nchans = min(parms.nchans,parms.corr_nchans_max);
  ind_chans = round(linspace(1,parms.nchans,nchans));
  % get indices of upper half of matrix
  [X,Y] = meshgrid(1:nchans,1:nchans);
  mask_upper = (X>Y);
  % loop over epochs, calculate correation between channels
  for i=1:parms.nepochs
    data = select_spindle(spindles,i,parms);
    R = corr(data');
    R2 = R.^2;
    results.corr.epochs.r(i) = ts_nan_mean(abs(R(mask_upper)));
    results.corr.epochs.r2(i) = ts_nan_mean(R2(mask_upper));
    results.corr.epochs.R(i,:,:) = R;
    results.corr.epochs.R2(i,:,:) = R2;
  end;

  %% todo: use Fischer z transform before calculating mean
  %%       and then transform back afterwards?
  
  results.corr.r_mean = mean(results.corr.epochs.r,1);
  results.corr.r_std = std(results.corr.epochs.r,0,1);
  results.corr.r2_mean = mean(results.corr.epochs.r2,1);
  results.corr.r2_std = std(results.corr.epochs.r2,0,1);
  results.corr.R_mean = squeeze(mean(results.corr.epochs.R,1));
  results.corr.R_std = squeeze(std(results.corr.epochs.R,0,1));
  results.corr.R2_mean = squeeze(mean(results.corr.epochs.R2,1));
  results.corr.R2_std = squeeze(std(results.corr.epochs.R2,0,1));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = select_spindle(spindles,n,parms)
  ind_first = max(1,spindles.ind_spindle_rel(n,1));
  ind_last = min(parms.ntpoints,spindles.ind_spindle_rel(n,2));
  data = squeeze(spindles.spindle_data(n,:,ind_first:ind_last));
  if parms.nchans==1
    data = data';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_coher(spindles,parms,results)
  % initialize output values
  results.coher = [];
  results.coher.epochs = zeros(parms.nepochs,parms.nchans,parms.nchans);
  % loop over epochs, calculate cross-channel coherence
  for i=1:parms.nepochs
    data = select_spindle(spindles,i,parms);
    ntpoints = size(data,2);
    nfft = 2^nextpow2(ntpoints);
    coher = ones(parms.nchans,parms.nchans);
    for j=1:parms.nchans
      for k=j+1:parms.nchans
        [coh,freq] = mscohere(data(j,:),data(k,:),[],0,nfft,parms.sfreq);
        [tmp,ind_freq_low] = min(abs(freq - parms.spindle_freq_low));
        [tmp,ind_freq_high] = min(abs(freq - parms.spindle_freq_high));
        ind_spindle_freq = ind_freq_low:ind_freq_high;
        coher(j,k) = mean(coh(ind_spindle_freq));
        coher(k,j) = coher(j,k);
      end;
    end;
    results.coher.epochs(i,:,:) = coher;
  end;
  results.coher.coh_mean = squeeze(mean(results.coher.epochs,1));
  results.coher.coh_std = squeeze(std(results.coher.epochs,0,1));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_rms(spindles,parms,results)
  % initialize output values
  results.rms = [];
  results.rms.epochs = zeros(parms.nepochs,parms.nchans);
  % loop over epochs, calculate root mean square amplitude for each channel
  for i=1:parms.nepochs
    data = select_spindle(spindles,i,parms);
    results.rms.epochs(i,:) = std(data,0,2);
  end;
  results.rms.rms_mean = mean(results.rms.epochs,1);
  results.rms.rms_std = std(results.rms.epochs,0,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_fft(spindles,parms,results,chans_flag)
  if ~exist('chans_flag','var'), chans_flag = 0; end;
  if chans_flag
    fname = 'fft_chans';
  else
    fname = 'fft';
  end;
  % initialize output values
  results.(fname) = [];
  results.(fname).epochs = [];
  peak_freq_chan = nan(parms.nepochs,results.nchans);
  peak_amp_chan = nan(parms.nepochs,results.nchans);
  peak_freq = nan(parms.nepochs,1);
  peak_amp = nan(parms.nepochs,results.nchans);
  peak_phase = nan(parms.nepochs,results.nchans);
  % loop over epochs, calculate fft for each channel
  for i=1:parms.nepochs
    if chans_flag
      chans = results.prob.spindle_chans{i};
    else
      chans = [1:results.nchans];
    end;
    data = select_spindle(spindles,i,parms);
    data = data(chans,:);
    % set number of time points for fft
    ntpoints = size(data,2);
    nfft = 2^nextpow2(ntpoints);
    % calculate Fourier transform
    fft_cx = fft(data,nfft,2)/ntpoints;
    fft_cx = fft_cx(:,1:nfft/2+1);
    % calculate frequencies
    freq = parms.sfreq/2*linspace(0,1,nfft/2+1);
    % calculate power
    amp = 2*abs(fft_cx);
    % scale power by frequency
    if parms.fft_norm_flag
      amp = bsxfun(@times,amp,freq.^parms.fft_norm_exp);
    end;
    % calculate phase (in cycles, not radians)
    phase = angle(fft_cx)/(2*pi);
    % calculate amp and phase in spindle freq range
    [tmp,ind_freq_low] = min(abs(freq - parms.spindle_freq_low));
    [tmp,ind_freq_high] = min(abs(freq - parms.spindle_freq_high));
    ind_spindle_freq = ind_freq_low:ind_freq_high;
    spindle_freq = freq(ind_spindle_freq);
    spindle_amp = amp(:,ind_spindle_freq);
    spindle_phase = phase(:,ind_spindle_freq);
    spindle_amp_mean = mean(spindle_amp,1);
    % find peak within spindle freq range for each channel
    [tmp,ind_peak] = max(spindle_amp,[],2);
    peak_freq_chan(i,chans) = spindle_freq(ind_peak);
    peak_amp_chan(i,chans) = spindle_amp(ind_peak);
    % find peak within spindle freq range for average across channels
    [tmp,ind_peak] = max(spindle_amp_mean,[],2);
    peak_freq(i) = spindle_freq(ind_peak);
    peak_phase(i,chans) = spindle_phase(:,ind_peak);
    peak_amp(i,chans) = spindle_amp(:,ind_peak);
    % save results
    results.(fname).epochs(i).chans = chans;
    results.(fname).epochs(i).freq = freq;
    results.(fname).epochs(i).amp = amp;
    results.(fname).epochs(i).phase = phase;
    results.(fname).epochs(i).peak_freq_chan = peak_freq_chan(i,chans);
    results.(fname).epochs(i).peak_amp_chan = peak_amp_chan(i,chans);
    results.(fname).epochs(i).peak_freq = peak_freq(i);
    results.(fname).epochs(i).peak_amp = peak_amp(i,chans);
    results.(fname).epochs(i).peak_phase = peak_phase(i,chans);
  end;
  % calculate mean and std of peak frequency, amp, and phase
  results.(fname).peak_freq_chan_mean = ts_nan_mean(peak_freq_chan,1);
  results.(fname).peak_freq_chan_std = ts_nan_std(peak_freq_chan,0,1);
  results.(fname).peak_amp_chan_mean = ts_nan_mean(peak_amp_chan,1);
  results.(fname).peak_amp_chan_std = ts_nan_std(peak_amp_chan,0,1);
  results.(fname).peak_freq_mean = ts_nan_mean(peak_freq);
  results.(fname).peak_freq_std = ts_nan_std(peak_freq);
  results.(fname).peak_amp_mean = ts_nan_mean(peak_amp,1);
  results.(fname).peak_amp_std = ts_nan_std(peak_amp,0,1);
  results.(fname).peak_phase_mean = ts_nan_mean(peak_phase,1);
  results.(fname).peak_phase_std = ts_nan_std(peak_phase,0,1);
  % calculate cross-channel correlation of peak fft measures
  if ~chans_flag && parms.nepochs
    results.(fname).peak_freq_chan_corr = corr(peak_freq_chan);
    results.(fname).peak_amp_chan_corr = corr(peak_amp_chan);
    results.(fname).peak_amp_corr = corr(peak_amp);
    results.(fname).peak_phase_corr = corr(peak_phase);
  end;

  % calculate phase differences between channels
  spindle_phase_mean_diff = zeros(results.nchans,results.nchans);
  spindle_phase_std_diff = zeros(results.nchans,results.nchans);
  spindle_phase_count_diff = zeros(results.nchans,results.nchans);
  for j=1:results.nchans
    for k=1:results.nchans
      tmp_diff = peak_phase(:,j) - peak_phase(:,k);
      % remove extra phase wrap
      tmp_diff(tmp_diff < -0.5) = tmp_diff(tmp_diff < -0.5) + 1;
      tmp_diff(tmp_diff > 0.5) = tmp_diff(tmp_diff > 0.5) - 1;
      % convert to complex, calculate mean and std, then back to phase
      tmp_r = cos(tmp_diff*2*pi);
      tmp_i = sin(tmp_diff*2*pi);
      tmp_r_mean = ts_nan_mean(tmp_r,1);
      tmp_i_mean = ts_nan_mean(tmp_i,1);
      tmp_r_std = ts_nan_std(tmp_r,1);
      tmp_i_std = ts_nan_std(tmp_i,1);
      tmp_diff_mean = atan2(tmp_i_mean,tmp_r_mean)/(2*pi);
      tmp_diff_std = sqrt(tmp_r_std^2 + tmp_i_std^2)/(2*pi);
%      spindle_phase_mean_diff(j,k) = ts_nan_mean(tmp_diff,1);
%      spindle_phase_std_diff(j,k) = ts_nan_std(tmp_diff,1);
      spindle_phase_mean_diff(j,k) = tmp_diff_mean;
      spindle_phase_std_diff(j,k) = tmp_diff_std;
      spindle_phase_count_diff(j,k) = length(find(~isnan(tmp_diff)));
    end;
  end;
  spindle_phase_zscore_diff = spindle_phase_mean_diff ./...
    (spindle_phase_std_diff./sqrt(spindle_phase_count_diff));
  spindle_phase_zscore_diff(isnan(spindle_phase_zscore_diff)) = 0;
  % save results
  results.(fname).peak_phase_mean_diff = spindle_phase_mean_diff;
  results.(fname).peak_phase_std_diff = spindle_phase_std_diff;
  results.(fname).peak_phase_count_diff = spindle_phase_count_diff;
  results.(fname).peak_phase_zscore_diff = spindle_phase_zscore_diff;

  % calculate average Fourier spectra for each channel
  nfreq = length(parms.fft_freq);  
  epochs = results.(fname).epochs;
  power_sum = zeros(results.nchans,nfreq);
  power_ssq = zeros(results.nchans,nfreq);
  power_n = zeros(results.nchans,1);
  for i=1:results.nepochs
    chans = epochs(i).chans;
    if isempty(chans), continue; end;
    tmp_freq = epochs(i).freq;
    tmp_amp = epochs(i).amp;
    amp = spline(tmp_freq,tmp_amp,parms.fft_freq);
    power_sum(chans,:) = power_sum(chans,:) + amp;
    power_ssq(chans,:) = power_ssq(chans,:) + amp.^2;
    power_n(chans) = power_n(chans) + 1;
  end;
  n = repmat(power_n,[1,nfreq]);
  results.(fname).spectra_freq = parms.fft_freq;
  results.(fname).spectra_mean = power_sum ./ n;
  results.(fname).spectra_std = sqrt(n.*power_ssq - power_sum.^2)./(n.*(n-1));
  results.(fname).spectra_n = power_n;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

