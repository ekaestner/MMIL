function results = ts_detect_spindles(data,varargin)
%function results = ts_detect_spindles(data,[options])
%
% Usage:
%  results = ts_detect_spindles(data,'key1', value1,...);
%
% Required Input:
%   data:  struct containing these fields (output of ts_prep_spindle_data):
%     sfreq: sampling frequency
%     labels: channel labels
%     time: time vector
%     data_raw: matrix of data with size = [nchans,ntpoints]
%     data_filt_narrow: matrix of bandpass-filtered data
%     data_filt_broad: matrix of bandpass-filtered data
%     spindle_power: matrix of estimated spindle power
%       (used for spindle edge detection)
%     spindle_power_smooth: matrix of estimated spindle power
%       with additional smoothing
%       (used for spindle peak detection)
%     low_power: matrix of estimated low band power
%     high_power: matrix of estimated low band power
%
% Optional input data:
%   'mask': mask vector with size = [1,ntpoints]
%     with values of 0 or 1 to exclude artifacts
%     {default = []}
%
% Optional parameters for spindle detection:
%  'peak_thresh': threshold for spindle peak detection
%    {default = 2}
%  'epoch_dur_init': epoch duration (sec)
%    for initial segmentation around the peak before edge detection
%    {default = 3}
%  'edge_thresh': threshold for spindle edge detection
%    {default = 1}
%  'relative_edge_thresh_flag': [0|1] determine edge threshold
%     based on a fraction of peak amplitude
%    {default = 1}
%  'relative_edge_thresh': fraction of peak amplitude
%     used instead of edge_thresh if relative_edge_thresh_flag = 1
%    {default = 0.5}
%  'dur_thresh': threshold for the minimum spindle duration in a channel (sec)
%    {default = 0.2}
%  'exclude_overlap_flag': [0|1] exclude smaller peaks
%                                within epochs of larger peaks
%    {default = 1}
%
% Optional parameters for spindle rejection:
%  'reject_chan_flag': reject individual channels from each spindle
%     otherwise, reject entire epoch
%    {default = 0}
%
% Optional parameters for spindle rejection based on non-spindle power:
%  'reject_low_flag': [0|1] reject epochs that have low_power
%    greater than low_band_thresh in any channel with that spindle
%    {default = 0}
%  'low_thresh': threshold for low_power rejection
%    {default = 5}
%  'reject_high_flag': [0|1] reject epochs that have high_power
%    greater than high_band_thresh in any channel with that spindle
%    {default = 0}
%  'high_thresh': threshold for high_power rejection
%    {default = 5}
%  'reject_band_any_chan_flag': [0|1] reject low or high band based
%     on ANY channel (instead of any channel with a spindle)
%    {default = 0}
%  'reject_band_onset_flag': [0|1] reject low or high band based
%     on overlap with spindle onset for any channel
%    {default = 0}
%
% Optional parameters for spindle rejection based on Fourier power ratio:
%  'reject_fourier_flag': [0|1] reject spindles with low power
%    in spindle frequency range or high power in non-spindle frequencies
%    {default = 0}
%  'spindle_freq_low': lower bound of spindle frequency range (in Hz)
%    {default = 10}
%  'spindle_freq_high': upper bound of spindle frequency range (in Hz)
%    {default = 16}
%  'nonspindle_freq_low': lower bound of non-spindle frequency range (in Hz)
%    {default = 5}
%  'nonspindle_freq_high': upper bound of non-spindle frequency range (in Hz)
%    {default = 8}
%  'fourier_thresh': threshold for the ratio between
%     spindle vs. non-spindle frequency power
%    {default = 0.5}
%
% Optional parameters for spindle rejection based on number of peaks:
%  'reject_npeaks_flag': [0|1] reject spindles with too few peaks
%    peak counting is based on data_filt_broad
%    {default = 0}
%  'min_npeaks': minimum number of peaks for an accepted spindle
%    {default = 3}
%  'reject_peakrate_flag': [0|1] reject spindles with too few peaks
%    peak counting is based on data_filt_broad
%    {default = 0}
%  'min_peakrate': minimum number of peaks per second for an accepted spindle
%    {default = 5}
%  'max_peakrate': maximum number of peaks per second for an accepted spindle
%    {default = 17}
%  'min_peak_deflection': fraction of channel-wise median abs deviation
%    if the difference between a local maximum and surrounding minima
%      is smaller than min_peak_deflection * med_abs
%      it is not considered a likely peak in the spindle oscillation
%    {default = 1}
%  'min_peak_ratio': exclude peaks from npeaks count if deflection
%     is smaller than min_peak_ratio * largest peak deflection for the epoch
%    {default = 0.25}
%
% Optional parameters for spindle rejection based on number of channels:
%  'reject_nchans_flag': [0|1] reject spindles with too few channels
%    {default = 0}
%  'min_nchans': minimum number of channels for an accepted spindle
%    {default = 2}
%
% Optional parameters for spindle extraction:
%  'trimmed_epochs_flag': [0|1] extracted spindles will have edges trimmed
%     epochs will have variable lengths
%    {default = 0}
%  'epoch_dur_extract': epoch duration (sec) for spindle extraction
%    unless trimmed_epochs_flag = 1
%    {default = 3}
%  'extract_reject_flag': [0|1] extract data for rejected spindles
%    {default = 0}
%  'extract_reject_types': cell array of types of rejects to be extracted
%    allowed: {'low','high','fourier','npeaks','peakrate','nchans'}
%    {default = {'low','high','fourier','npeaks','peakrate','nchans'}}
%     
% Output:
%   results: struct containing these fields:
%     nchans: number of channels of input data
%     ntpoints: number of time points of input data
%     sfreq: sampling frequency of data
%     spindle_data: matrix or cell array containing extracted data epochs
%       if trimmed_epochs_flag = 0, 
%         matrix with size [nepochs,nchans,ntpoints_epoch]
%         (fixed epoch duration)
%       if trimmed_epochs_flag = 1, cell array with size [nepochs,1]
%         containing data matrices of size  = [nchans,ntpoints_epoch]
%         (variable epoch duration)
%     spindle_data_filt_broad: matrix or cell aray containing extract epochs
%       with broad bandpass-filtered data
%     spindle_data_filt_narrow: matrix or cell aray containing extract epochs
%       with bandpass-filtered data
%     spindle_power_peak: matrix or cell aray containing extract epochs
%       with spindle power from peak wave
%     spindle_power_edge: matrix or cell aray containing extract epochs
%       with spindle power from edge wave
%     spindle_mask: matrix with size = [nchans,ntpoints]
%       containing zeros or ones
%       to indicate when and where spindles were detected
%     spindle_chans: cell array with size = [nepochs,1]
%       containing vectors with indices of channels with spindles detected
%     ind_spindle_chan: struct array with nchans elements
%       and three fields: first, middle, and last
%       indicating time samples of spindles detected for each channel
%     ind_spindle: matrix with size = [nepochs,2]
%       with first and last time points (samples) for each epoch
%       containing a spindle in one or more channels
%     ind_spindle_rel: matrix with size = [nepochs,2]
%       with first and last time points (samples) for each epoch (all channels)
%       relative to start of extracted epoch
%     -- optional fields created if extract_reject_flag:
%     reject_data: matrix or cell array containing extracted data epochs
%       if trimmed_epochs_flag = 0, matrix with size [nepochs,ntpoints_epoch]
%         (fixed epoch duration)
%       if trimmed_epochs_flag = 1, cell array with size [nepochs,1]
%         containing data matrices of size  = [nchans,ntpoints_epoch]
%         (variable epoch duration)
%     reject_power_peak: matrix or cell aray containing extract epochs
%       with rejected spindle power from peak wave
%     reject_power_edge: matrix or cell aray containing extract epochs
%       with rejected spindle power from edge wave
%     reject_mask: matrix with size = [nchans,ntpoints]
%       containing zeros or ones
%       to indicate when and where rejected spindles were detected
%     reject_chans: cell array with size = [nepochs,1]
%       containing vectors with indices of channels with rejected spindles detected
%     reject_types: cell array with size = [nepochs,1]
%       containing reasons for epoch rejection
%     ind_reject: matrix with size = [nepochs,2]
%       with first and last time points (samples) for each epoch
%       containing a rejected spindle in one or more channels
%     ind_reject_rel: matrix with size = [nepochs,2]
%       with first and last time points (samples) for each epoch (all channels)
%       relative to start of extracted epoch
%
% Created:  08/22/14 by Don Hagler
% Last mod: 08/10/15 by Don Hagler
%
% Based on code by August Tan and Xi Jiang
%

%% todo: exclude spindle epochs that include masked time points?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters
[parms,data] = check_input(data,varargin);

% initialize results struct
results = init_results(parms);

% calculate median absolute deviation for each channel
results.data_mad = mmil_wtd_mad(data.data_raw,[],2);
results.data_filt_broad_mad = mmil_wtd_mad(data.data_filt_broad,[],2);
results.data_filt_narrow_mad = mmil_wtd_mad(data.data_filt_narrow,[],2);

% find spindles for each channel
results.ind_spindle_chan = find_spindles(parms,data);

% label time points containing spindle epochs in each channel
results.spindle_mask_orig = ...
  create_spindle_mask_from_chans(parms,results.ind_spindle_chan);

% find spindle epochs across all channels
[results.ind_spindle,results.spindle_chans] = ...
  find_spindle_epochs(parms,results.spindle_mask_orig);

% reject spindle epochs
[results.ind_spindle,results.spindle_chans,...
 results.ind_reject,results.reject_chans,results.reject_types] = ...
    reject_spindle_epochs(parms,data,...
                          results.ind_spindle,results.spindle_chans);

if ~isempty(results.ind_reject)
  % update spindle_mask
  results.spindle_mask = ...
    create_spindle_mask(parms,results.ind_spindle,results.spindle_chans);
  results.spindle_mask = results.spindle_mask .* results.spindle_mask_orig;
  % create reject_mask
  results.reject_mask = ...
    create_spindle_mask(parms,results.ind_reject,results.reject_chans);
else
  results.spindle_mask = results.spindle_mask_orig;
end;

% extract spindle epochs from data matrix
[results.spindle_data,...
 results.spindle_data_filt_broad,results.spindle_data_filt_narrow,...
 results.spindle_power_peak,results.spindle_power_edge,...
 results.ind_spindle_rel,results.ind_spindle] = ...
   extract_spindles(parms,data,results.ind_spindle);

% calculate basic spindle metrics
results.spindle_metrics = ...
  calc_metrics(parms,data,results.ind_spindle,results.spindle_chans);

if parms.extract_reject_flag && ~isempty(results.reject_types)
  % determine which rejects to extract
  ind_extr = ...
    find(ismember(results.reject_types,parms.extract_reject_types));

  % save original rejects
  results.ind_reject_orig = results.ind_reject;
  results.reject_chans_orig = results.reject_chans;
  results.reject_types_orig = results.reject_types;

  % copy extracted rejects
  results.ind_reject = results.ind_reject(ind_extr,:);
  results.reject_chans = results.reject_chans(ind_extr);
  results.reject_types = results.reject_types(ind_extr);

  % extract rejected spindle epochs from data matrix
  [results.reject_data,...
   results.reject_data_filt_broad,results.reject_data_filt_narrow,...
   results.reject_power_peak,results.reject_power_edge,...
   results.ind_reject_rel,results.ind_reject] = ...
    extract_spindles(parms,data,results.ind_reject);

  % calculate basic spindle metrics for rejects
  results.reject_metrics = ...
    calc_metrics(parms,data,results.ind_reject,results.reject_chans);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = check_input(data,options)
  parms = mmil_args2parms(options,{...
  ... % spindle detection
    'mask',[],[],...
    'peak_thresh',2,[0,Inf],...
    'epoch_dur_init',3,[0.2,Inf],...
    'edge_thresh',1,[0,Inf],...
    'relative_edge_thresh_flag',true,[false true],...
    'relative_edge_thresh',0.5,[0,1],...
    'dur_thresh',0.2,[0,Inf],...
    'exclude_overlap_flag',true,[false true],...
  ... % spindle rejection
    'reject_chan_flag',false,[false true],...
  ... % spindle rejection based on non-spindle power
    'reject_low_flag',false,[false true],...
    'low_thresh',5,[0,Inf],...
    'reject_high_flag',false,[false true],...
    'high_thresh',5,[0,Inf],...
    'reject_band_any_chan_flag',false,[false true],...
    'reject_band_onset_flag',false,[false true],...
  ... % spindle rejection based on Fourier ratio
    'reject_fourier_flag',false,[false true],...
    'spindle_freq_low',10,[0.5,100],...
    'spindle_freq_high',16,[0.5,100],...
    'nonspindle_freq_low',5,[0,100],...
    'nonspindle_freq_high',8,[0,100],...
    'fourier_thresh',0.5,[0,Inf],...
  ... % spindle rejection based on number of peaks
    'reject_npeaks_flag',false,[false true],...
    'min_npeaks',3,[0,100],...
    'reject_peakrate_flag',false,[false true],...
    'min_peakrate',5,[0,100],...
    'max_peakrate',17,[0,100],...
    'min_peak_deflection',1,[0,Inf],...
    'min_peak_ratio',0.25,[0,Inf],...
  ... % spindle rejection based on number of channels
    'reject_nchans_flag',false,[false true],...
    'min_nchans',2,[2,Inf],...
  ... % spindle extraction
    'epoch_dur_extract',3,[0.2,Inf],...
    'trimmed_epochs_flag',false,[false true],...
    'extract_reject_flag',false,[false true],...
    'extract_reject_types',{'low','high','fourier','npeaks','peakrate','nchans'},...
                           {'low','high','fourier','npeaks','peakrate','nchans'},...
  ...
    'verbose',true,[false true],...
    'required_fields',{'sfreq','data_raw','data_filt_narrow','data_filt_broad',...
                       'spindle_power','spindle_power_smooth'},[],...
  });

  % check input data
  if ~isstruct(data)
    error('input data not a structure as expected');
  end;
  for i=1:length(parms.required_fields)
    fstr = parms.required_fields{i};
    if ~isfield(data,fstr)
      error('input data struct missing required field %s',fstr);
    end;
  end;

  % check data matrix
  if numel(size(data.data_raw))~=2
    error('data matrix must be 2-dimensional');
  end;
  if any(size(data.data_raw)~=size(data.spindle_power_smooth)) ||...
     any(size(data.data_raw)~=size(data.spindle_power))
    error('size of spindle_power_smooth and spindle_power must match data');
  end;
  [parms.nchans,parms.ntpoints] = size(data.data_raw);
  [data.nchans,data.ntpoints] = size(data.data_raw);

  % check mask
  if ~isempty(parms.mask)
    data.mask = mmil_rowvec(parms.mask);
    if length(data.mask)~=parms.ntpoints
      error('length of mask must match 2nd dim of data');
    end;
  else
    data.mask = ones(1,parms.ntpoints);
  end;

  % check extract_reject_types
  if ~iscell(parms.extract_reject_types)
    parms.extract_reject_types = {parms.extract_reject_types};
  end;

  % check chan labels
  if ~isfield(data,'labels')
    data.labels = cell(parms.nchans,1);
    for i=1:parms.nchans
      data.labels{i} = sprintf('chan%d',i);
    end;
  else
    if ~iscell(data.labels), data.labels = {data.labels}; end;
    if length(data.labels) ~= parms.nchans
      error('labels has wrong number of elements (%d not %d)\n',...
        length(data.labels),parms.nchans);
    end;
  end;
  parms.labels = data.labels;

  % set sfreq for parms
  parms.sfreq = data.sfreq;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results.labels = parms.labels;
  results.nchans = parms.nchans;
  results.ntpoints = parms.ntpoints;
  results.sfreq = parms.sfreq;
  results.spindle_data = [];
  results.spindle_data_filt_broad = [];
  results.spindle_data_filt_narrow = [];
  results.spindle_power_peak = [];
  results.spindle_power_edge = [];
  results.spindle_metrics = [];
  results.spindle_mask = [];
  results.spindle_chans = [];
  results.ind_spindle_chan = [];
  results.ind_spindle = [];
  results.ind_spindle_rel = [];
  if parms.extract_reject_flag
    results.reject_data = [];
    results.reject_data_filt_broad = [];
    results.reject_data_filt_narrow = [];
    results.reject_power_peak = [];
    results.reject_power_edge = [];
    results.reject_metrics = [];
    results.reject_mask = [];
    results.reject_chans = [];
    results.reject_types = [];
    results.ind_reject = [];
    results.ind_reject_rel = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_spindle_chan = find_spindles(parms,data)
  clear ind_spindle_chan;
  if parms.verbose
    fprintf('%s: finding spindles...\n',mfilename);
  end;
  for i=1:parms.nchans
    % get data time course for this channel
    data_chan = data.data_raw(i,:);
    % copy filtered data for counting peaks
    data_filt_chan = data.data_filt_broad(i,:);
    % get spindle power wave time courses
    peak_wave_chan = data.spindle_power_smooth(i,:).*data.mask;
    edge_wave_chan = data.spindle_power(i,:).*data.mask;
    % intialize spindle epochs
    ind_epoch = init_epochs(parms,peak_wave_chan);
    % detect edges to refine epoch boundaries
    ind_epoch = detect_edges(parms,ind_epoch,edge_wave_chan);
    % save spindles for this channel
    ind_spindle_chan(i) = ind_epoch;
    % exclude short spindles, duplicates, and resolve overlap
    ind_spindle_chan(i) = exclude_spindles(parms,ind_spindle_chan(i),...
        data_chan,data_filt_chan,peak_wave_chan);
  end
  if ~exist('ind_spindle_chan','var'), ind_spindle_chan = []; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = init_epochs(parms,peak_wave_chan)
  ind_epochs = [];
  % find maxima in peak wave timecourse
  [peak_vals,ind_epochs.middle,~,~] = mmil_extrema(peak_wave_chan);
  % exclude negative and smaller peaks, then sort by time point
  ind_epochs.middle = ...
    sort(ind_epochs.middle(peak_vals > parms.peak_thresh));
  % set first and last indices for epochs
  nsamples_half_epoch = round(0.5*parms.epoch_dur_init*parms.sfreq);
  ind_epochs.first = max(ind_epochs.middle - nsamples_half_epoch,1);
  ind_epochs.last  = min(ind_epochs.middle + nsamples_half_epoch,...
                         length(peak_wave_chan));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = detect_edges(parms,ind_epochs,edge_wave_chan)
  nepochs = length(ind_epochs.middle);
  for i=1:nepochs
    % get edge_wave for this epoch
    ind_middle = ind_epochs.middle(i);
    ind_first = ind_epochs.first(i);
    ind_last = ind_epochs.last(i);
    epoch_edge_wave = edge_wave_chan(ind_first:ind_last);
    % relative index of peak
    ind_peak = ind_middle - ind_first + 1;
    % find peak in edge_wave nearest to peak in peak_wave
    [peak_vals,ind_peaks,~,~] = mmil_extrema(epoch_edge_wave);
    [~,ind_nearest] = min(abs(ind_peaks-ind_peak));
    ind_peak = ind_peaks(ind_nearest);
    % use relative or fixed threshold
    if parms.relative_edge_thresh_flag
      edge_thresh = parms.relative_edge_thresh * epoch_edge_wave(ind_peak);
    else
      edge_thresh = parms.edge_thresh;
    end;
    % find time points below threshold
    ind_subthresh = find(epoch_edge_wave < edge_thresh);
    % find sub-threshold points nearest to peak
    ind_edge_first = max(ind_subthresh(ind_subthresh<ind_peak));
    ind_edge_last = min(ind_subthresh(ind_subthresh>ind_peak));
    % save new edge indices
    if ~isempty(ind_edge_first)
      ind_epochs.first(i) = ind_first + ind_edge_first - 1;
    end;
    if ~isempty(ind_edge_last)
      ind_epochs.last(i) = ind_first + ind_edge_last - 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = exclude_spindles(parms,...
                         ind_epochs,data,data_filt,peak_wave_chan)
  % exclude short duration epochs
  ind_epochs = exclude_short_spindles(parms,ind_epochs);
  % exclude duplicate epochs
  ind_epochs = exclude_duplicate_spindles(parms,ind_epochs,peak_wave_chan);
  % check for overlap, exclude smaller peaks within epochs of larger peaks
  if parms.exclude_overlap_flag
    ind_epochs = exclude_overlap(parms,ind_epochs,peak_wave_chan);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = exclude_short_spindles(parms,ind_epochs)
  ind_epochs_edges = cat(2,ind_epochs.first',ind_epochs.last');
  dur = calc_dur(parms,ind_epochs_edges);
  ind_keep = find(dur >= parms.dur_thresh);
  ind_epochs = update_spindles(ind_epochs,ind_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = exclude_duplicate_spindles(parms,...
                                                 ind_epochs,peak_wave_chan)
  nepochs = length(ind_epochs.middle);
  keep_flags = ones(1,nepochs);
  for i=1:nepochs
    if ~keep_flags(i), continue; end;
    % find other spindles with the same first and last
    ind_overlap = find(keep_flags &...
      ind_epochs.first == ind_epochs.first(i) &...
      ind_epochs.last == ind_epochs.last(i));
    if length(ind_overlap)<2, continue; end;
    % exclude smaller peaks
    peak_vals = peak_wave_chan(ind_epochs.middle(ind_overlap));
    [tmp,ind_max] = max(peak_vals);
    ind_excl = setdiff(ind_overlap,ind_overlap(ind_max));
    keep_flags(ind_excl) = 0;
  end;
  ind_keep = find(keep_flags);
  ind_epochs = update_spindles(ind_epochs,ind_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = exclude_overlap(parms,ind_epochs,peak_wave_chan)
  nepochs = length(ind_epochs.middle);
  keep_flags = ones(1,nepochs);
  for i=1:nepochs
    if ~keep_flags(i), continue; end;
    % find other spindles within this epoch
    ind_overlap = find(keep_flags &...
      ind_epochs.middle >= ind_epochs.first(i) &...
      ind_epochs.middle <= ind_epochs.last(i));
    if length(ind_overlap)<2, continue; end;
    % exclude smaller peaks
    peak_vals = peak_wave_chan(ind_epochs.middle(ind_overlap));
    [tmp,ind_max] = max(peak_vals);
    ind_excl = setdiff(ind_overlap,ind_overlap(ind_max));
    keep_flags(ind_excl) = 0;
  end;
  ind_keep = find(keep_flags);
  ind_epochs = update_spindles(ind_epochs,ind_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_epochs = update_spindles(ind_epochs,ind_keep)
  nepochs = length(ind_epochs.middle);
  if length(ind_keep) < nepochs
    ind_epochs.first = ind_epochs.first(ind_keep);
    ind_epochs.middle = ind_epochs.middle(ind_keep);
    ind_epochs.last = ind_epochs.last(ind_keep);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle_mask = create_spindle_mask_from_chans(parms,ind_spindle_chan)
  spindle_mask = zeros(parms.nchans,parms.ntpoints);
  for i=1:parms.nchans
    spindle_chan_mask = zeros(1,parms.ntpoints);
    nepochs = length(ind_spindle_chan(i).middle);
    for j=1:nepochs
      ind_epoch = ind_spindle_chan(i).first(j):ind_spindle_chan(i).last(j);
      spindle_chan_mask(ind_epoch) = 1;
    end 
    spindle_mask(i,:) = spindle_chan_mask;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spindle_mask = create_spindle_mask(parms,ind_spindle,spindle_chans)
  spindle_mask = zeros(parms.nchans,parms.ntpoints);
  nepochs = size(ind_spindle,1);
  for i=1:nepochs
    spindle_mask(spindle_chans{i},ind_spindle(i,1):ind_spindle(i,2)) = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans] = find_spindle_epochs(parms,spindle_mask)
  ind_spindle = []; spindle_chans = [];

  % mask of spindles across all channels
  spindle_vec = 1.0*(sum(spindle_mask,1)>0);
  % find spindle onsets and offsets
  spindle_diff = [0 diff(spindle_vec)];
  ind_first = find(spindle_diff == 1)';
  ind_last = find(spindle_diff == -1)';
  if spindle_vec(1)==1
    % first spindle at beginning of data
    ind_first = cat(1,1,ind_first);  
  end;
  if spindle_vec(end)==1
    % last spindle at end of data
    ind_last = cat(1,ind_last,parms.ntpoints);
  end;
  % check that onsets match offsets
  nepochs = length(ind_first);
  if nepochs ~= length(ind_last)
    error('number of spindle onsets (%d) does not match number of spindle offsets (%d)',...
      nepochs,length(ind_last));
  end;
  % set start and end times for each epoch
  ind_spindle = cat(2,ind_first,ind_last);
  % check for spindles too close to the edges
  if ~parms.trimmed_epochs_flag
    % calculate mid-point
    ind_middle = round(mean(ind_spindle,2));
    % calculate epoch start and end points
    ntpoints_epoch = round(parms.epoch_dur_extract*parms.sfreq);
    ind_first_epoch = ind_middle - round(0.5*ntpoints_epoch);
    ind_last_epoch = ind_first_epoch + ntpoints_epoch - 1;
    % exclude spindles too close to the edges
    ind_exclude = find(ind_first_epoch<0 | ind_last_epoch>parms.ntpoints);
    if ~isempty(ind_exclude)
      ind_spindle(ind_exclude,:) = [];
      nepochs = size(ind_spindle,1);
      ind_first = ind_spindle(:,1);
      ind_last = ind_spindle(:,2);
    end;
  end;
  % determine which channels have a detected spindle for each epoch
  spindle_chans = cell(nepochs,1);
  for i=1:nepochs
    ind_chans = [];
    for j=1:parms.nchans
      if any(spindle_mask(j,ind_first(i):ind_last(i)))
        ind_chans(end+1) = j;
      end;
    end;
    spindle_chans{i} = ind_chans;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    reject_spindle_epochs(parms,data,ind_spindle,spindle_chans)
  if parms.verbose
    fprintf('%s: rejecting spindles...\n',mfilename);
  end;
  % initialize output
  ind_reject = [];
  reject_chans = [];
  reject_types = [];
  % optionally reject spindles with suprathreshold low-band amplitudes
  if parms.reject_low_flag
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_low_band(parms,data,ind_spindle,spindle_chans,...
                                 ind_reject,reject_chans,reject_types);
  end;
  % optionally reject spindles with suprathreshold high-band amplitudes
  if parms.reject_high_flag
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_high_band(parms,data,ind_spindle,spindle_chans,...
                                  ind_reject,reject_chans,reject_types);
  end;
  % reject spindles without a significant peak in the spindle frequency range
  if parms.reject_fourier_flag
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_fourier_power(parms,data,ind_spindle,spindle_chans,...
                                      ind_reject,reject_chans,reject_types);
  end;
  % reject spindles based on number of peaks
  if parms.reject_npeaks_flag || parms.reject_peakrate_flag
    % reject spindles with too few peaks
    %              or with peak rate that is too low or too high
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_npeaks(parms,data,ind_spindle,spindle_chans,...
                               ind_reject,reject_chans,reject_types);
  end;
  % reject spindles based on number of channels
  if parms.reject_nchans_flag
    % reject spindles with too few channels
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_nchans(parms,data,ind_spindle,spindle_chans,...
                               ind_reject,reject_chans,reject_types);
  end;
  % sort rejects
  if ~isempty(ind_reject)
    [~,ind_sort] = sort(ind_reject(:,1));
    ind_reject = ind_reject(ind_sort,:);
    reject_chans = reject_chans(ind_sort);
    reject_types = reject_types(ind_sort);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    update_rejects(keep_flags,reject_type,...
                   ind_spindle,spindle_chans,...
                   ind_reject,reject_chans,reject_types)
  nepochs = length(ind_spindle);
  % save rejects
  ind_excl = find(~keep_flags);
  if length(ind_excl)>0
    ind_reject = cat(1,ind_reject,ind_spindle(ind_excl,:));
    reject_chans = cat(1,reject_chans,spindle_chans(ind_excl));
    reject_types = cat(1,reject_types,repmat({reject_type},[length(ind_excl),1]));
  end;
  % save remaining spindles
  ind_keep = find(keep_flags);
  if length(ind_keep) < nepochs
    ind_spindle = ind_spindle(ind_keep,:);
    spindle_chans = spindle_chans(ind_keep);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_low_band(parms,data,ind_spindle,spindle_chans,...
                                 ind_reject,reject_chans,reject_types)
  nepochs = size(ind_spindle,1);
  keep_flags = zeros(nepochs,1);
  for i=1:nepochs
    % option to exclude based on any channel or all spindle channels
    if parms.reject_band_any_chan_flag
      chans = 1:data.nchans;
    else
      chans = spindle_chans{i};
    end;
    ind0 = ind_spindle(i,1);
    ind1 = ind_spindle(i,2);
    % check for overlap between band and first half of spindle
    if parms.reject_band_onset_flag
      ind1 = min(parms.ntpoints,ind0 + round((ind1 - ind0)/2));
    end;
    low_power = data.low_power(chans,ind0:ind1);
    nchans = length(chans);
    excl_flags = (max(low_power,[],2) >= parms.low_thresh);
    % option to exclude individual channnels or entire epoch
    if parms.reject_chan_flag
      nexcl = length(find(excl_flags));
      if nexcl < nchans
        keep_flags(i) = 1;
        if nexcl > 0
          spindle_chans{i} = chans(excl_flags==0);
        end;
      end;
    elseif ~any(excl_flags)
      keep_flags(i) = 1;
    end;
  end;
  % exclude rejected epochs
  [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    update_rejects(keep_flags,'low',ind_spindle,spindle_chans,...
                   ind_reject,reject_chans,reject_types);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_high_band(parms,data,ind_spindle,spindle_chans,...
                                  ind_reject,reject_chans,reject_types)
  nepochs = size(ind_spindle,1);
  keep_flags = zeros(nepochs,1);
  for i=1:nepochs
    % option to exclude based on any channel or all spindle channels
    if parms.reject_band_any_chan_flag
      chans = 1:data.nchans;
    else
      chans = spindle_chans{i};
    end;
    ind0 = ind_spindle(i,1);
    ind1 = ind_spindle(i,2);
    % check for overlap between band and first half of spindle
    if parms.reject_band_onset_flag
      ind1 = min(parms.ntpoints,ind0 + round((ind1 - ind0)/2));
    end;
    high_power = data.high_power(chans,ind0:ind1);
    nchans = length(chans);
    excl_flags = (max(high_power,[],2) >= parms.high_thresh);
    % option to exclude individual channnels or entire epoch
    if parms.reject_chan_flag
      nexcl = length(find(excl_flags));
      if nexcl < nchans
        keep_flags(i) = 1;
        if nexcl > 0
          spindle_chans{i} = chans(excl_flags==0);
        end;
      end;
    elseif ~any(excl_flags)
      keep_flags(i) = 1;
    end;
  end;
  % exclude rejected epochs
  [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    update_rejects(keep_flags,'high',ind_spindle,spindle_chans,...
                   ind_reject,reject_chans,reject_types);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_fourier_power(parms,data,ind_spindle,spindle_chans,...
                                  ind_reject,reject_chans,reject_types)
  nepochs = size(ind_spindle,1);
  keep_flags = zeros(nepochs,1);
  % calculate ratio between spindle and non-spindle power
  F_ratio = calc_F_ratio(parms,data,ind_spindle,spindle_chans);
  % find epochs with high enough F_ratio
  for i=1:nepochs
    chans = spindle_chans{i};
    nchans = length(chans);
    excl_flags = (F_ratio(i,chans) < parms.fourier_thresh);
    % option to exclude individual channnels or entire epoch
    if parms.reject_chan_flag
      nexcl = length(find(excl_flags));
      if nexcl < nchans
        keep_flags(i) = 1;
        if nexcl > 0
          spindle_chans{i} = chans(excl_flags==0);
        end;
      end;
    elseif ~any(excl_flags)
      keep_flags(i) = 1;
    end;
  end;
  % exclude rejected epochs
  [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    update_rejects(keep_flags,'fourier',ind_spindle,spindle_chans,...
                   ind_reject,reject_chans,reject_types);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_ratio,spindle_freq,spindle_power,nonspindle_power] =...
      calc_F_ratio(parms,data,ind_spindle,spindle_chans)
  nepochs = size(ind_spindle,1);
  nchans = data.nchans;
  F_ratio = nan(nepochs,nchans);
  spindle_freq = nan(nepochs,nchans);
  spindle_power = nan(nepochs,nchans);
  nonspindle_power = nan(nepochs,nchans);
  for i=1:nepochs
    % set chans
    chans = spindle_chans{i};
    nchans = length(chans);
    % extract epoch data
    epoch_data = data.data_raw(chans,ind_spindle(i,1):ind_spindle(i,2));
    % remove linear trend
    epoch_data = detrend(epoch_data')';
    % calculate number of time points and fft samples
    ntpoints = size(epoch_data,2);
    nfft = 2^nextpow2(ntpoints);
    % calculate frequencies
    fft_freqs = parms.sfreq/2*linspace(0,1,nfft/2+1);
    % calculate indices for frequencies of interest
    [tmp,ind_freq_low] = min(abs(fft_freqs - parms.spindle_freq_low));
    [tmp,ind_freq_high] = min(abs(fft_freqs - parms.spindle_freq_high));
    ind_spindle_freqs = ind_freq_low:ind_freq_high;
    [tmp,ind_freq_low] = min(abs(fft_freqs - parms.nonspindle_freq_low));
    [tmp,ind_freq_high] = min(abs(fft_freqs - parms.nonspindle_freq_high));
    ind_nonspindle_freqs = ind_freq_low:ind_freq_high;
    for j=1:nchans
      c = chans(j);
      % calculate Fourier power
      fft_power = fft(epoch_data(j,:),nfft)/ntpoints;
      fft_power = 2*abs(fft_power(1:nfft/2+1));
      % find peak within spindle freq range for each channel      
      tmp_power = fft_power(ind_spindle_freqs);
      [tmp,ind_peak] = max(tmp_power,[],2);
      spindle_freq(i,c) = fft_freqs(ind_spindle_freqs(ind_peak));
      % max power within spindle frequency range
      spindle_power(i,c) = fft_power(ind_spindle_freqs(ind_peak));
      % max power within non-spindle frequency range
      nonspindle_power(i,c) = max(fft_power(ind_nonspindle_freqs));
    end;
  end;
  % calculate ratio of spindle power to non-spindle power
  F_ratio = spindle_power ./ nonspindle_power;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_npeaks(parms,data,ind_spindle,spindle_chans,...
                               ind_reject,reject_chans,reject_types)
  % calculate number of peaks for each spindle
  [npeaks,peakrate] = calc_npeaks(parms,data,ind_spindle,spindle_chans);
  % find epochs with high enough number of peaks
  if parms.reject_npeaks_flag
    nepochs = size(ind_spindle,1);
    keep_flags = zeros(nepochs,1);
    for i=1:nepochs
      chans = spindle_chans{i};
      nchans = length(chans);
      excl_flags = (npeaks(i,chans) < parms.min_npeaks);
      % option to exclude individual channnels or entire epoch
      if parms.reject_chan_flag
        nexcl = length(find(excl_flags));
        if nexcl < nchans
          keep_flags(i) = 1;
          if nexcl > 0
            spindle_chans{i} = chans(excl_flags==0);
          end;
        end;
      elseif ~any(excl_flags)
        keep_flags(i) = 1;
      end;
    end;
    % exclude rejected epochs
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      update_rejects(keep_flags,'npeaks',ind_spindle,spindle_chans,...
                     ind_reject,reject_chans,reject_types);
    % update peakrate
    ind_keep = find(keep_flags);
    npeaks = npeaks(ind_keep,:);
    peakrate = peakrate(ind_keep,:);
  end;
  % find epochs with high enough peak rate
  if parms.reject_peakrate_flag
    nepochs = size(ind_spindle,1);
    keep_flags = zeros(nepochs,1);
    for i=1:nepochs
      chans = spindle_chans{i};
      nchans = length(chans);
      excl_flags = (peakrate(i,chans) < parms.min_peakrate |...
                    peakrate(i,chans) > parms.max_peakrate);
      % option to exclude individual channnels or entire epoch
      if parms.reject_chan_flag
        nexcl = length(find(excl_flags));
        if nexcl < nchans
          keep_flags(i) = 1;
          if nexcl > 0
            spindle_chans{i} = chans(excl_flags==0);
          end;
        end;
      elseif ~any(excl_flags)
        keep_flags(i) = 1;
      end;
    end;
    % exclude rejected epochs
    [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      update_rejects(keep_flags,'peakrate',ind_spindle,spindle_chans,...
                     ind_reject,reject_chans,reject_types);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [npeaks,peakrate,peak_intv_mean,peak_intv_std] = ...
     calc_npeaks(parms,data,ind_spindle,spindle_chans)
  nepochs = size(ind_spindle,1);
  nchans = data.nchans;
  npeaks = nan(nepochs,nchans);
  peakrate = nan(nepochs,nchans);
  peak_intv_mean = nan(nepochs,nchans);
  peak_intv_std = nan(nepochs,nchans);
  data_mad = mmil_wtd_mad(data.data_filt_broad,[],2);
  for i=1:nepochs
    chans = spindle_chans{i};
    nchans = length(chans);
    % extract epoch data
    epoch_data = data.data_filt_broad(chans,ind_spindle(i,1):ind_spindle(i,2));
    % remove linear trend
    epoch_data = detrend(epoch_data')';
    % calculate duration
    dur = (1+diff(ind_spindle(i,:)))/parms.sfreq;
    for j=1:nchans
      c = chans(j);
      dmad = max(data_mad(c),eps);
      tmp_data = epoch_data(j,:);
      % find local maxima with large enough peak to peak amplitudes
      [maxtab,mintab] = mmil_peakdet(tmp_data,...
                                     parms.min_peak_deflection*dmad);
      tmp_minima = get_minmax(mintab);
      tmp_maxima = get_minmax(maxtab);
      tmp_minima = get_deflection(tmp_minima,tmp_maxima);
      tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
      % exclude peaks with deflection much smaller than max peak
      peak_vals = tmp_maxima.deflection;
      peak_lats = tmp_maxima.latency;
      max_peak_val = max(peak_vals);
      ind_keep = find(peak_vals >= parms.min_peak_ratio*max_peak_val);
      % save number of peaks
      npeaks(i,c) = length(ind_keep);
      % normalize number of peaks by duration to get peaks per second
      peakrate(i,c) = npeaks(i,c)/dur;
      % calculate mean and std of spacing between peaks
      peak_lats = peak_lats(ind_keep);
      peak_intv = diff(peak_lats);
      peak_intv_mean(i,c) = mean(peak_intv);
      peak_intv_std(i,c) = std(peak_intv);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax(mtab)
  minmax = struct('latency',[],'amplitude',[],'deflection',[]);
  j = 1;
  for i=1:size(mtab,1)
    latency = mtab(i,1);
    amplitude = mtab(i,2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
      reject_nchans(parms,data,ind_spindle,spindle_chans,...
                                  ind_reject,reject_chans,reject_types)
  nepochs = size(ind_spindle,1);
  keep_flags = zeros(nepochs,1);
  % find epochs with high enough number of channels
  for i=1:nepochs
    chans = spindle_chans{i};
    nchans = length(chans);
    if nchans >= parms.min_nchans
      keep_flags(i) = 1;
    end;
  end;
  % exclude rejected epochs
  [ind_spindle,spindle_chans,ind_reject,reject_chans,reject_types] = ...
    update_rejects(keep_flags,'nchans',ind_spindle,spindle_chans,...
                   ind_reject,reject_chans,reject_types);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spindle_data,...
          spindle_data_filt_broad,spindle_data_filt_narrow,...
          spindle_power_peak,spindle_power_edge,...
          ind_spindle_rel,ind_spindle] =...
          extract_spindles(parms,data,ind_spindle)
  nepochs = size(ind_spindle,1);
  ind_spindle_rel = zeros(nepochs,2);
  if parms.trimmed_epochs_flag
    spindle_data = cell(nepochs,1);
    spindle_data_filt_broad = cell(nepochs,1);
    spindle_data_filt_narrow = cell(nepochs,1);
    spindle_power_peak = cell(nepochs,1);
    spindle_power_edge = cell(nepochs,1);
  else
    ntpoints_epoch = round(parms.epoch_dur_extract*parms.sfreq);
    spindle_data = zeros(nepochs,parms.nchans,ntpoints_epoch);
    spindle_data_filt_broad = zeros(nepochs,parms.nchans,ntpoints_epoch);
    spindle_data_filt_narrow = zeros(nepochs,parms.nchans,ntpoints_epoch);
    spindle_power_peak = zeros(nepochs,parms.nchans,ntpoints_epoch);
    spindle_power_edge = zeros(nepochs,parms.nchans,ntpoints_epoch);
  end;
  if isempty(data.data_filt_broad)
    spindle_data_filt_broad = [];
  end;
  if isempty(data.data_filt_narrow)
    spindle_data_filt_narrow = [];
  end;
  for i=1:size(ind_spindle,1)
    ind_first = ind_spindle(i,1);
    ind_last = ind_spindle(i,2);
    if parms.trimmed_epochs_flag
      % calculate variable spindle length
      ntpoints_epoch = ind_last - ind_first + 1;
      % extract data
      spindle_data{i} = data.data_raw(:,ind_first:ind_last);
      % extract bandpass-filtered data
      if ~isempty(data.data_filt_broad)
        spindle_data_filt_broad{i} = data.data_filt_broad(:,ind_first:ind_last);
      end;
      if ~isempty(data.data_filt_narrow)
        spindle_data_filt_narrow{i} = data.data_filt_narrow(:,ind_first:ind_last);
      end;
      % extract spindle power
      spindle_power_peak{i} = data.spindle_power_smooth(:,ind_first:ind_last);
      spindle_power_edge{i} = data.spindle_power(:,ind_first:ind_last);
      % relative indices
      ind_spindle_rel(i,1) = 1;
      ind_spindle_rel(i,2) = ntpoints_epoch;
    else
      ind_middle = round(mean([ind_first,ind_last]));
      ind_first_epoch = ind_middle - round(0.5*ntpoints_epoch);
      ind_last_epoch = ind_first_epoch + ntpoints_epoch - 1;
      % check that first and last are within bounds of data
      ind_first_check = max(ind_first_epoch,1);
      ind_last_check  = min(ind_last_epoch,parms.ntpoints);
      ind_first_check_rel = ind_first_check - ind_first_epoch + 1;
      ind_last_check_rel = ind_last_check - ind_first_epoch + 1;
      % extract data
      spindle_data(i,:,ind_first_check_rel:ind_last_check_rel) =...
        data.data_raw(:,ind_first_check:ind_last_check);
      % extract bandpass-filtered data
      if ~isempty(data.data_filt_broad)
        spindle_data_filt_broad(i,:,ind_first_check_rel:ind_last_check_rel) =...
          data.data_filt_broad(:,ind_first_check:ind_last_check);
      end;
      if ~isempty(data.data_filt_narrow)
        spindle_data_filt_narrow(i,:,ind_first_check_rel:ind_last_check_rel) =...
          data.data_filt_narrow(:,ind_first_check:ind_last_check);
      end;
      % extract spindle power
      spindle_power_peak(i,:,ind_first_check_rel:ind_last_check_rel) =...
        data.spindle_power_smooth(:,ind_first_check:ind_last_check);
      spindle_power_edge(i,:,ind_first_check_rel:ind_last_check_rel) =...
        data.spindle_power(:,ind_first_check:ind_last_check);
      % calculate relative indices
      ind_spindle_rel(i,1) = ind_first - ind_first_epoch + 1;
      ind_spindle_rel(i,2) = ind_last - ind_first_epoch + 1;
      % make sure bounds are within extracted epoch
      if ind_first < ind_first_epoch
        ind_spindle(i,1) = ind_first_epoch;
        ind_spindle_rel(i,1) = 1;
      end;
      if ind_last > ind_last_epoch
        ind_spindle(i,2) = ind_last_epoch;
        ind_spindle_rel(i,2) = ntpoints_epoch;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metrics = calc_metrics(parms,data,ind_spindle,spindle_chans)
  metrics = [];
  nepochs = size(ind_spindle,1);
  if nepochs==0, return; end;
  if parms.verbose
    fprintf('%s: calculating metrics...\n',mfilename);
  end;
  metrics.amp = calc_amp(parms,data,ind_spindle);
  metrics.dur = calc_dur(parms,ind_spindle);
  [metrics.F_ratio,metrics.spindle_freq,...
   metrics.spindle_power,metrics.nonspindle_power] =...
      calc_F_ratio(parms,data,ind_spindle,spindle_chans);
  [metrics.npeaks,metrics.peakrate,...
   metrics.peak_intv_mean,metrics.peak_intv_std] = ...
      calc_npeaks(parms,data,ind_spindle,spindle_chans);
  % calculate max spindle band, low band, and high band power/amplitude
  metrics.sband_amp = calc_amp(parms,data,ind_spindle,'spindle_power');
  metrics.lband_amp = calc_amp(parms,data,ind_spindle,'low_power');
  metrics.hband_amp = calc_amp(parms,data,ind_spindle,'high_power');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function amp = calc_amp(parms,data,ind_spindle,fstr)
  if ~exist('fstr','var') || isempty(fstr)
    fstr = 'spindle_power_smooth';
  end;
  nepochs = size(ind_spindle,1);
  amp = nan(nepochs,data.nchans);
  for i=1:nepochs
    % extract epoch data
    epoch_data = data.(fstr)(:,ind_spindle(i,1):ind_spindle(i,2));
    amp(i,:) = max(epoch_data,[],2);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dur = calc_dur(parms,ind_epochs_edges)
  dur = (1+diff(ind_epochs_edges,[],2))/parms.sfreq;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

