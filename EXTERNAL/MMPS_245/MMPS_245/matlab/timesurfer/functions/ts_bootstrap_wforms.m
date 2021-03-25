function results = ts_bootstrap_wforms(wforms,varargin)
%function results = ts_bootstrap_wforms(wforms,[options])
%
% Purpose: calculate confidence intervals for average waveforms,
%   peak amplitude, peak latency, onset latency, and area under curve
%
% Required Input:
%   wforms: waveform matrix
%     size of wforms should be [ntpoints,nsubs] or [ntpoints,nsubs,2]
%     if 3rd dim has 2 elements, will treat as conditions to be compared
%
% Optional Resampling Parameters:
%   'niters': number of iterations
%     {default = 1000}
%   'randseed_flag': [0|1] generate different result every time
%     {default = 0}
%   'bias_corr_flag': bootstrap bias correction
%     {default = 1}
%
% Optional Measurement Parameters:
%   'peak_flag': [0|1] find main peak and calculate amplitude and latency
%     {default = 0}
%   'peak_pol': [-1,1] peak polarity
%     {default = 1}
%   'time': time vector used for peak, onset, and auc
%     {default = []}
%   'onset_flag': [0|1] find response onset
%     {default = 0}
%   'auc_flag': [0|1] calculate integrated power (area under curve)
%     {default = 0}
%   'wform_ci_flag': [0|1] calculate waveform confidence intervals
%     {default = 1}
%   'alpha': probability threshold
%     {default = 0.05}
%
% Output:
%   results: struct containing output
%
% Created:  02/06/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'niters',1000,[1 Inf],...
  'randseed_flag',false,[false true],...
  'bias_corr_flag',true,[false true],...
...
  'peak_flag',false,[false true],...
  'onset_flag',false,[false true],...
  'auc_flag',false,[false true],...
  'wform_ci_flag',true,[false true],...
  'alpha',0.05,[eps,1],...
... % for peaks, onset, and auc
  'sfreq',[],[],...
  't0',[],[],...
  't1',[],[],...
  'time',[],[],...
... % for peaks
  'smooth_sigma',[],[],...
  'peak_range',[],[],...
  'peak_pol',1,[-1,1,1],...
  'peak_mindiff',[],[],...
  'firstflag',[],[],...
  'resfact',[],[],...
... % for onset
  'onset_method',[],[],...
  'onset_minimum',[],[],...
  'onset_kappa',[],[],...
  'onset_baseline_range',[],[],...
  'onset_baseline_flag',[],[],...
  'onset_polarity',[],[],...
  'onset_baseline_collapse_flag',[],[],...
... % for auc
  'auc_range',[],[],...
  'auc_nbins',1,[],...
  'auc_baseline_flag',[],[],...
  'auc_baseline_range',[],[],...
  'powerflag',[],[],...
...
  'peaks_tags',{'smooth_sigma','peak_range','peak_pol','peak_mindiff',...
                'firstflag','resfact'...
                'sfreq','t0','t1','time'},[],...
  'onset_tags',{'onset_method','onset_minimum','onset_kappa',...
                'onset_baseline_range','onset_baseline_flag','onset_polarity',...
                'onset_baseline_collapse_flag','sfreq','t0','t1','time'},[],...
  'auc_tags',{'auc_range','auc_nbins','auc_baseline_flag',...
              'auc_baseline_range','powerflag','sfreq',...
              't0','t1','time'},[],...
};
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);
results.time = parms.time;

% calculate parametric statistics (mean and sem)
results = calc_parametric_stats(parms,results);

% resample waveforms
results = resamp_wforms(parms,results);

% calculate statistics from resampled waveforms
results = calc_resamp_stats(parms,results);

% calculate waveform confidence intervals
if parms.wform_ci_flag
  results = calc_wform_ci(parms,results);
end;

% calculate peak amplitude and latency with confidence intervals
if parms.peak_flag
  results = calc_peak_ci(parms,results);
end;

% calculate onset latency with confidence intervals
if parms.onset_flag
  results = calc_onset_ci(parms,results);
end;

% calculate area under curve confidence intervals
if parms.auc_flag
  results = calc_auc_ci(parms,results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  if numel(size(wforms))>3
    error('wforms has too many dimensions');
  end;
  if size(wforms,3)>2
    error('wforms has too many elements in 3rd dim');
  end;
  parms.ntpoints = size(wforms,1);
  parms.nsubs = size(wforms,2);
  if parms.nsubs<=1
    error('wforms must have more than one subjects (2nd dim)');
  end;
  if size(wforms,3)==1
    parms.wforms = wforms;
    parms.wforms1 = [];
    parms.wforms2 = [];
  elseif size(wforms,3)==2
    parms.wforms1 = wforms(:,:,1);
    parms.wforms2 = wforms(:,:,2);
    parms.wforms = parms.wforms1 - parms.wforms2;
  end;
  parms.dof = parms.nsubs - 1;
  parms.ci_tail = 100*parms.alpha/2;
  if parms.randseed_flag
    seedval = sum(100*clock);
  else
    seedval = 5489; % default
  end;
  stream = RandStream.create('mt19937ar','Seed',seedval);
  RandStream.setDefaultStream(stream);
  parms.pval_range = linspace(1/parms.niters,1,parms.niters);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_parametric_stats(parms,results)
  results.wforms_mean = mean(parms.wforms,2);
  results.wforms_sem =  std(parms.wforms,0,2)/sqrt(parms.nsubs);
  if ~isempty(parms.wforms1)
    results.wforms1_mean = mean(parms.wforms1,2);
    results.wforms1_sem =  std(parms.wforms1,0,2)/sqrt(parms.nsubs);
    results.wforms2_mean = mean(parms.wforms2,2);
    results.wforms2_sem =  std(parms.wforms2,0,2)/sqrt(parms.nsubs);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = resamp_wforms(parms,results)
  % resampling with replacement
  ind_resamp = randi(parms.nsubs,parms.nsubs,parms.niters);
  ind_resamp = reshape(ind_resamp,[1,parms.nsubs*parms.niters]);
  results.wforms_resamp = reshape(parms.wforms(:,ind_resamp),...
                      [parms.ntpoints,parms.nsubs,parms.niters]);
  if ~isempty(parms.wforms1)
    results.wforms1_resamp = reshape(parms.wforms1(:,ind_resamp),...
                         [parms.ntpoints,parms.nsubs,parms.niters]);
    results.wforms2_resamp = reshape(parms.wforms2(:,ind_resamp),...
                         [parms.ntpoints,parms.nsubs,parms.niters]);
  end;
  % jackknife resampling to calculate skew
  if parms.bias_corr_flag
    ind_jn = zeros(parms.nsubs-1,parms.nsubs);
    for s=1:parms.nsubs
      ind_jn(:,s) = setdiff([1:parms.nsubs],s);
    end;
    ind_jn = reshape(ind_jn,[1,(parms.nsubs-1)*parms.nsubs]);
    results.wforms_jackknife = reshape(parms.wforms(:,ind_jn),...
                        [parms.ntpoints,parms.nsubs-1,parms.nsubs]);
    if ~isempty(parms.wforms1)
      results.wforms1_jackknife = reshape(parms.wforms1(:,ind_jn),...
                           [parms.ntpoints,parms.nsubs-1,parms.nsubs]);
      results.wforms2_jackknife = reshape(parms.wforms2(:,ind_jn),...
                           [parms.ntpoints,parms.nsubs-1,parms.nsubs]);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_resamp_stats(parms,results)
  % calculate mean and sem for each iteration
  results.wforms_resamp_mean = squeeze(mean(results.wforms_resamp,2));
  if ~isempty(parms.wforms1)
    results.wforms1_resamp_mean = squeeze(mean(results.wforms1_resamp,2));
    results.wforms2_resamp_mean = squeeze(mean(results.wforms2_resamp,2));
  end;
  % calculate mean and stdev of mean estimates across iterations
  results.wforms_bs_mean = mean(results.wforms_resamp_mean,2);
  results.wforms_bs_sem = std(results.wforms_resamp_mean,0,2);
  if ~isempty(parms.wforms1)
    results.wforms1_bs_mean = mean(results.wforms1_resamp_mean,2);
    results.wforms1_bs_sem = std(results.wforms1_resamp_mean,0,2);
    results.wforms2_bs_mean = mean(results.wforms2_resamp_mean,2);
    results.wforms2_bs_sem = std(results.wforms2_resamp_mean,0,2);
  end;

  if parms.bias_corr_flag
    results.wforms_jackknife_mean = squeeze(mean(results.wforms_jackknife,2));
    if ~isempty(parms.wforms1)
      results.wforms1_jackknife_mean = squeeze(mean(results.wforms1_jackknife,2));
      results.wforms2_jackknife_mean = squeeze(mean(results.wforms2_jackknife,2));
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_wform_ci(parms,results)
  % calculate confidence intervals for mean wforms
  if parms.bias_corr_flag
    results.wforms_bs_ci = ts_calc_bs_ci_bias(results.wforms_resamp_mean,...
      results.wforms_jackknife_mean,results.wforms_mean,parms.alpha);
  else
    results.wforms_bs_ci = ...
      ts_calc_bs_ci(results.wforms_resamp_mean,parms.alpha);
  end;
  % calculate confidence intervals for mean wforms1 and wforms2
  if ~isempty(parms.wforms1)
    if parms.bias_corr_flag
      results.wforms1_bs_ci = ...
        ts_calc_bs_ci_bias(results.wforms1_resamp_mean,...
          results.wforms1_jackknife_mean,results.wforms1_mean,parms.alpha);
    else
      results.wforms1_bs_ci = ...
        ts_calc_bs_ci(results.wforms1_resamp_mean,parms.alpha);
    end;
    if parms.bias_corr_flag
      results.wforms2_bs_ci = ...
        ts_calc_bs_ci_bias(results.wforms2_resamp_mean,...
          results.wforms2_jackknife_mean,results.wforms2_mean,parms.alpha);
    else
      results.wforms2_bs_ci = ...
        ts_calc_bs_ci(results.wforms2_resamp_mean,parms.alpha);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_peak_ci(parms,results)
  tmp_results = [];
  if isempty(parms.wforms1)
    nr = 1;
  else
    nr = 2;
  end;
  for r=1:nr
    if isempty(parms.wforms1)
      wforms_tmp = reshape(results.wforms_resamp_mean,...
        [parms.ntpoints,1,parms.niters]);
    elseif r==1
      wforms_tmp = reshape(results.wforms1_resamp_mean,...
        [parms.ntpoints,1,parms.niters]);
    else
      wforms_tmp = reshape(results.wforms2_resamp_mean,...
        [parms.ntpoints,1,parms.niters]);
    end;
    args = mmil_parms2args(parms,parms.peaks_tags);
    peak = ts_wform_peaks(wforms_tmp,args{:});
    tmp_results(r).peak.amplitudes = peak.amplitude;
    tmp_results(r).peak.latencies = peak.latency;
    % calculate mean amplitude and latency
    tmp_results(r).peak.amplitude = ts_nan_mean(peak.amplitude);
    tmp_results(r).peak.latency = ts_nan_mean(peak.latency);
    % calculate bootstrap standard error of the mean of amplitude and latency
    tmp_results(r).peak.amplitude_sem = ts_nan_std(peak.amplitude,0,2);
    tmp_results(r).peak.latency_sem = ts_nan_std(peak.latency,0,2);
    % calculate confidence intervals for amplitude and latency
    if parms.bias_corr_flag
      % peaks from jackknife resampled mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = reshape(results.wforms_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      elseif r==1
        wforms_tmp = reshape(results.wforms1_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      else
        wforms_tmp = reshape(results.wforms2_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      end;  
      peak_jackknife = ts_wform_peaks(wforms_tmp,args{:});
      tmp_results(r).peak.amplitudes_jackknife = peak_jackknife.amplitude;
      tmp_results(r).peak.latencies_jackknife = peak_jackknife.latency;
      % peaks from mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = results.wforms_mean;
      elseif r==1
        wforms_tmp = results.wforms1_mean;
      else
        wforms_tmp = results.wforms2_mean;
      end;
      peak_mean = ts_wform_peaks(wforms_tmp,args{:});
      tmp_results(r).peak.amplitude_mean = peak_mean.amplitude;
      tmp_results(r).peak.latency_mean = peak_mean.latency;
      amp_ci = ts_calc_bs_ci_bias(peak.amplitude,...
        peak_jackknife.amplitude,peak_mean.amplitude,parms.alpha);
      lat_ci = ts_calc_bs_ci_bias(peak.latency,...
        peak_jackknife.latency,peak_mean.latency,parms.alpha);
    else
      amp_ci = ts_calc_bs_ci(peak.amplitude,parms.alpha);
      lat_ci = ts_calc_bs_ci(peak.latency,parms.alpha);
    end;
    tmp_results(r).peak.amplitude_ci = amp_ci;
    tmp_results(r).peak.latency_ci = lat_ci;
  end;

  if isempty(parms.wforms1)
    results.peak = tmp_results.peak;
  else
    results.peak1 = tmp_results(1).peak;
    results.peak2 = tmp_results(2).peak;
    % calculate differences in amplitudes and latencies
    results.peak = [];
    results.peak.amplitudes = results.peak1.amplitudes - results.peak2.amplitudes;
    results.peak.latencies = results.peak1.latencies - results.peak2.latencies;
    results.peak.amplitude = ts_nan_mean(results.peak.amplitudes);
    results.peak.latency = ts_nan_mean(results.peak.latencies);
    % calculate bootstrap standard error of the mean of amplitude and latency
    results.peak.amplitude_sem = ts_nan_std(results.peak.amplitude,0,2);
    results.peak.latency_sem = ts_nan_std(results.peak.latency,0,2);
    if parms.bias_corr_flag
      results.peak.amplitudes_jackknife =...
         results.peak1.amplitudes_jackknife - results.peak2.amplitudes_jackknife;
      results.peak.latencies_jackknife = ...
         results.peak1.latencies_jackknife - results.peak2.latencies_jackknife;
      results.peak.amplitude_mean =...
         results.peak1.amplitude_mean - results.peak2.amplitude_mean;
      results.peak.latency_mean = ...
         results.peak1.latency_mean - results.peak2.latency_mean;
      % calculate confidence interval for amplitude difference
      results.peak.amplitude_ci = ...
        ts_calc_bs_ci_bias(results.peak.amplitudes,...
          results.peak.amplitudes_jackknife,...
          results.peak.amplitude_mean,parms.alpha);
      % find p-value for null hypothesis for amplitude difference
      results.peak.amplitude_pval = ...
        ts_calc_bs_null_pval_bias(results.peak.amplitudes,...
          results.peak.amplitudes_jackknife,...
          results.peak.amplitude_mean,parms.pval_range);
      % calculate confidence interval for latency difference
      results.peak.latency_ci = ...
        ts_calc_bs_ci_bias(results.peak.latencies,...
          results.peak.latencies_jackknife,...
          results.peak.latency_mean,parms.alpha);
      % find p-value for null hypothesis for latency difference
      results.peak.latency_pval = ...
        ts_calc_bs_null_pval_bias(results.peak.latencies,...
          results.peak.latencies_jackknife,...
          results.peak.latency_mean,parms.pval_range);
    else
      % calculate confidence interval for amplitude difference
      results.peak.amplitude_ci = ...
        ts_calc_bs_ci(results.peak.amplitudes,parms.alpha);
      % find p-value for null hypothesis for amplitude difference
      results.peak.amplitude_pval = ...
        ts_calc_bs_null_pval(results.peak.amplitudes,parms.pval_range);
      % calculate confidence interval for latency difference
      results.peak.latency_ci = ...
        ts_calc_bs_ci(results.peak.latencies,parms.alpha);
      % find p-value for null hypothesis for latency difference
      results.peak.latency_pval = ...
        ts_calc_bs_null_pval(results.peak.latencies,parms.pval_range);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_onset_ci(parms,results)
  tmp_results = [];
  if isempty(parms.wforms1)
    nr = 1;
  else
    nr = 2;
  end;
  for r=1:nr
    if isempty(parms.wforms1)
      wforms_tmp = reshape(results.wforms_resamp_mean,[parms.ntpoints,1,parms.niters]);
    elseif r==1
      wforms_tmp = reshape(results.wforms1_resamp_mean,[parms.ntpoints,1,parms.niters]);
    else
      wforms_tmp = reshape(results.wforms2_resamp_mean,[parms.ntpoints,1,parms.niters]);
    end;
    args = mmil_parms2args(parms,parms.onset_tags);
    onset = ts_wform_onset(wforms_tmp,args{:});
    tmp_results(r).onset.latencies = onset;
    % calculate mean latency
    tmp_results(r).onset.latency = ts_nan_mean(onset);
    % calculate bootstrap standard error of the mean
    tmp_results(r).onset.latency_sem = ts_nan_std(onset,0,2);
    % calculate confidence intervals for onset latency
    if parms.bias_corr_flag
      % onset from jackknife resampled mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = reshape(results.wforms_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      elseif r==1
        wforms_tmp = reshape(results.wforms1_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      else
        wforms_tmp = reshape(results.wforms2_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      end;  
      onset_jackknife = ts_wform_onset(wforms_tmp,args{:});
      tmp_results(r).onset.latencies_jackknife = onset_jackknife;
      % onset from mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = results.wforms_mean;
      elseif r==1
        wforms_tmp = results.wforms1_mean;
      else
        wforms_tmp = results.wforms2_mean;
      end;
      onset_mean = ts_wform_onset(wforms_tmp,args{:});
      tmp_results(r).onset.latency_mean = onset_mean;
      lat_ci = ts_calc_bs_ci_bias(onset,...
        onset_jackknife,onset_mean,parms.alpha);
    else
      lat_ci = ts_calc_bs_ci(onset.latency,parms.alpha);
    end;
    tmp_results(r).onset.latency_ci = lat_ci;
  end;
  if isempty(parms.wforms1)
    results.onset = tmp_results.onset;
  else
    results.onset1 = tmp_results(1).onset;
    results.onset2 = tmp_results(2).onset;
    % calculate differences in onset latencies
    results.onset = [];
    results.onset.latencies = results.onset1.latencies - results.onset2.latencies;
    results.onset.latency = ts_nan_mean(results.onset.latencies);
    % calculate bootstrap standard error of the mean
    results.onset.latency_sem = ts_nan_std(results.onset.latencies,0,2);
    if parms.bias_corr_flag
      results.onset.latencies_jackknife = ...
         results.onset1.latencies_jackknife - results.onset2.latencies_jackknife;
      results.onset.latency_mean = ...
         results.onset1.latency_mean - results.onset2.latency_mean;
      % calculate confidence interval for latency difference
      results.onset.latency_ci = ...
        ts_calc_bs_ci_bias(results.onset.latencies,...
          results.onset.latencies_jackknife,...
          results.onset.latency_mean,parms.alpha);
      % find p-value for null hypothesis for latency difference
      results.onset.latency_pval = ...
        ts_calc_bs_null_pval_bias(results.onset.latencies,...
          results.onset.latencies_jackknife,...
          results.onset.latency_mean,parms.pval_range);
    else
      % calculate confidence interval for latency difference
      results.onset.latency_ci = ...
        ts_calc_bs_ci(results.onset.latencies,parms.alpha);
      % find p-value for null hypothesis for latency difference
      results.onset.latency_pval = ...
        ts_calc_bs_null_pval(results.onset.latencies,parms.pval_range);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_auc_ci(parms,results)
  tmp_results = [];
  if isempty(parms.wforms1)
    nr = 1;
  else
    nr = 2;
  end;
  for r=1:nr
    if isempty(parms.wforms1)
      wforms_tmp = reshape(results.wforms_resamp_mean,[parms.ntpoints,1,parms.niters]);
    elseif r==1
      wforms_tmp = reshape(results.wforms1_resamp_mean,[parms.ntpoints,1,parms.niters]);
    else
      wforms_tmp = reshape(results.wforms2_resamp_mean,[parms.ntpoints,1,parms.niters]);
    end;  
    args = mmil_parms2args(parms,parms.auc_tags);
    auc = ts_wform_auc(wforms_tmp,args{:});
    tmp_results(r).auc.vals = auc;
    % calculate mean auc
    tmp_results(r).auc.mean = ts_nan_mean(auc);
    % calculate bootstrap standard error of the mean
    tmp_results(r).auc.sem = ts_nan_std(auc,0,2);
    % calculate confidence intervals for auc
    if parms.bias_corr_flag
      % auc from jackknife resampled mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = reshape(results.wforms_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      elseif r==1
        wforms_tmp = reshape(results.wforms1_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      else
        wforms_tmp = reshape(results.wforms2_jackknife_mean,...
          [parms.ntpoints,1,parms.nsubs]);
      end;  
      auc_jackknife = ts_wform_auc(wforms_tmp,args{:});
      tmp_results(r).auc.vals_jackknife = auc_jackknife;
      % auc from mean wforms
      if isempty(parms.wforms1)
        wforms_tmp = results.wforms_mean;
      elseif r==1
        wforms_tmp = results.wforms1_mean;
      else
        wforms_tmp = results.wforms2_mean;
      end;
      auc_mean = ts_wform_auc(wforms_tmp,args{:});
      tmp_results(r).auc.val_mean = auc_mean;
      auc_ci = ts_calc_bs_ci_bias(auc,...
        auc_jackknife,auc_mean,parms.alpha);
    else
      auc_ci = ts_calc_bs_ci(auc,parms.alpha);
    end;
    tmp_results(r).auc.ci = auc_ci;
  end;
  if isempty(parms.wforms1)
    results.auc = tmp_results.auc;
  else
    results.auc1 = tmp_results(1).auc;
    results.auc2 = tmp_results(2).auc;
    % calculate differences in auc
    results.auc = [];
    results.auc.vals = results.auc1.vals - results.auc2.vals;
    results.auc.vals = results.auc1.vals - results.auc2.vals;
    % calculate mean difference in auc
    results.auc.mean = ts_nan_mean(results.auc.vals);
    % calculate bootstrap standard error of the mean
    results.auc.sem = ts_nan_std(results.auc.vals,0,2);
    if parms.bias_corr_flag
      results.auc.vals_jackknife =...
         results.auc1.vals_jackknife - results.auc2.vals_jackknife;
      results.auc.val_mean =...
         results.auc1.val_mean - results.auc2.val_mean;
      % calculate confidence interval for amplitude difference
      results.auc.ci = ...
        ts_calc_bs_ci_bias(results.auc.vals,...
          results.auc.vals_jackknife,...
          results.auc.val_mean,parms.alpha);
      % find p-value for null hypothesis for amplitude difference
      results.auc.pval = ...
        ts_calc_bs_null_pval_bias(results.auc.vals,...
          results.auc.vals_jackknife,...
          results.auc.val_mean,parms.pval_range);
    else
      % calculate confidence interval for amplitude difference
      results.auc.ci = ...
        ts_calc_bs_ci(results.auc.vals,parms.alpha);
      % find p-value for null hypothesis for amplitude difference
      results.auc.pval = ...
        ts_calc_bs_null_pval(results.auc.vals,parms.pval_range);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


