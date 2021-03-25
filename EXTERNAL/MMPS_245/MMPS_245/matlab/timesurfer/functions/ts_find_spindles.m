function [spindle_results,spindle_data] = ts_find_spindles(data,sfreq,varargin)
%function [spindle_results,spindle_data] = ts_find_spindles(data,sfreq,[options])
%
% Usage:
%  [spindle_results,spindle_data] = ...
%     ts_find_spindles(data,sfreq,'key1', value1,...);
%
% Required Input:
%   data: data matrix with size = [nchans,ntpoints]
%   sfreq: sampling frequency (in Hz)
%
% Optional parameters for preparing data:
%  'nonspindle_corr_flag': [0|1] adjust spindle power time courses using
%     time courses of non-spindle power (low_band and high_band)
%     if offset_flag = 0, calculate ratio
%     if offset_flag = 1, calculate difference
%    {default = 1}
%  'power_flag': [0|1] how to estimate power in frequency bands
%     if 0, calculate amplitude (abs)
%     if 1, calculate power (^2)
%    {default = 0}
%  'presmooth_win': moving average window length (sec)
%     for smoothing band-pass filtered envelopes
%     used to create spindle_power
%     not used if wavelet_flag = 1
%    {default = 0.3}
%  'postsmooth_win': moving average window length (sec)
%     for smoothing estimated spindle likelihood time courses
%     used to create spindle_power_smooth from spindle_power
%     not used if wavelet_flag = 1
%    {default = 0.6}
%  'wavelet_flag': [0|1] estimate spindle amplitude through wavelet convolution
%    using each frequency within spindle_band
%    if power_flag = 1, convolution results will be squared
%    {default = 0}
%  'wavelet_dur': duration of the wavelet (sec)
%    {default = 1}
%  'wavelet_peak_width': wavelet width for peak detection (spindle_power_smooth)
%    {default = 0.6}
%  'wavelet_peak_win': moving average window length (sec)
%    of wavelet convolution for peak detection (spindle_power_smooth)
%    {default = 0.3}
%  'wavelet_edge_width': wavelet width for edge detection (spindle_power)
%    higher value = broader peak in freq. space
%    {default = 0.2}
%  'wavelet_edge_win': moving average window length (sec)
%    of wavelet convolution for edge detection (spindle_power)
%    {default = 0.2}
%  'offset_flag': [0|1] offset power time courses
%    by subtracting median (for each channel)
%    {default = 1}
%  'offset_chanavg_flag': [0|1] offset median averaged across channels
%    only applies if offset_flag = 1
%    {default = 0}
%  'norm_flag': [0|1] normalize power time courses
%    by dividing by median absolute deviation (for each channel)
%    {default = 1}
%  'norm_chanavg_flag': [0|1] normalize by MAD averaged across channels
%    only applies if norm_flag = 1
%    {default = 0}
%  'norm_postsmooth_flag': [0|1] re-normalize spindle amplitude after smoothing
%    {default = 0}
%  'low_band_weight': weighting to give to low_band power
%    in calculating non-spindle power
%    may range from 0 to 1; 0.5 means equal weighting with high_band
%    {default = 0.5}
%  'spindle_band': bounds of spindle frequency range (Hz)
%    {default = [10,16]}
%  'low_band': bounds of low-frequency non-spindle frequency range (Hz)
%    only used if nonspindle_corr_flag = 1
%    {default = [4,8]}
%  'high_band': bounds of high-frequency non-spindle frequency range (Hz)
%    {default = [18,30]}
%    only used if nonspindle_corr_flag = 1
%  'broad_band': bounds of broad-band frequency range (Hz)
%    {default = [2,25]}
%  'filt_tfrac': width of filter transition bands relative to cut-off frequency
%    {default = 0.3}
%  'butter_flag': use Butterworth filter instead of FFT-based filter
%    {default = 0}
%  'butter_order': order of Butterworth filter
%    {default = 3}
%
% Optional parameters for spindle detection:
%  'mask': mask vector with size = [1,ntpoints]
%     with values of 0 or 1 to exclude artifacts
%     {default = []}
%  'peak_thresh': threshold for spindle peak detection
%    applied to spindle_power_smooth (see ts_prep_spindle_data)
%    {default = 3}
%  'epoch_dur_init': epoch duration (sec)
%    for initial segmentation around the peak before edge detection
%    {default = 4}
%  'edge_thresh': threshold for spindle edge detection
%    applied to spindle_power (see ts_prep_spindle_data)
%    {default = 3}
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
%  'reject_low_flag': [0|1] reject epochs that have low_power
%    greater than low_band_thresh in ANY channel
%    {default = 0}
%  'low_thresh': threshold for low_power rejection
%    {default = 5}
%  'reject_high_flag': [0|1] reject epochs that have high_power
%    greater than high_band_thresh in ANY channel
%    {default = 0}
%  'high_thresh': threshold for high_power rejection
%    {default = 5}
%  'reject_band_any_chan_flag': [0|1] reject low or high band based
%     on ANY channel (instead of any channel with a spindle)
%    {default = 0}
%  'reject_band_onset_flag': [0|1] reject low or high band based
%     on overlap with spindle onset for any channel
%    {default = 0}
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
%  'reject_npeaks_flag': [0|1] reject spindles with too few peaks
%    {default = 0}
%  'min_npeaks': minimum number of peaks per second for an accepted spindle
%    {default = 3}
%  'reject_peakrate_flag': [0|1] reject spindles with too few peaks
%    peak counting is based on data_filt_broad if supplied
%      or if not, either data_filt_narrow if supplied or data
%    {default = 0}
%  'min_peakrate': minimum number of peaks per second for an accepted spindle
%    {default = 5}
%  'max_peakrate': maximum number of peaks per second for an accepted spindle
%    {default = 17}
%  'min_peak_deflection': fraction of channel-wise median abs deviation
%    if the difference between a local maximum and surrounding minima
%      is smaller than min_peak_deflection * MAD
%      it is not considered a likely peak in the spindle oscillation
%    {default = 1}
%  'min_peak_ratio': exclude peaks from npeaks count if deflection
%     is smaller than min_peak_ratio * largest peak deflection for the epoch
%    {default = 0.25}
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
%  'labels': cell array of channel labels
%    if empty will set to channel1, channel2, etc.
%    {default = []}
%
% Output:
%   spindle_results: struct containing these fields:
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
%   spindle_data: struct containing these fields:
%     sfreq: sampling frequency
%     labels: channel labels
%     time: time vector
%     data_raw: input data
%     data_filt_narrow: matrix of bandpass-filtered data (size of data)
%       using spindle_band
%     data_filt_broad: matrix of bandpass-filtered data (size of data)
%       using broad_band
%     spindle_power: matrix of estimated spindle power (or amplitude)
%     spindle_power_smooth: matrix of estimated spindle power (or amplitude)
%       with additional smoothing (postsmooth_win)
%     low_power: matrix of estimated low band power (or amplitude)
%     high_power: matrix of estimated low band power (or amplitude)
%     nonspindle_power: matrix of estimated non-spindle power (or amplitude)
%       weighted average of low_power and high_power (see low_band_weight)
%       included only if nonspindle_corr_flag = 1
%     spindle_power_orig: matrix of estimated non-spindle power (or amplitude)
%       before normalization with non-spindle power
%       included only if nonspindle_corr_flag = 1
%
% Created:  08/22/14 by Don Hagler
% Last mod: 08/10/15 by Don Hagler
%

%% todo: exclude spindle epochs that include masked time points?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spindle_results = []; spindle_data = [];
if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters
parms = check_input(data,varargin);

% prepare data filtering and calculating spindle power
args = mmil_parms2args(parms,parms.prep_tags);
spindle_data = ts_prep_spindle_data(data,sfreq,args{:});

% detect spindles
args = mmil_parms2args(parms,parms.detect_tags);
spindle_results = ts_detect_spindles(spindle_data,args{:});

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(data,options)
  parms = mmil_args2parms(options,{...
  ... % data preparation
    'nonspindle_corr_flag',true,[false true],...
    'power_flag',false,[false true],...
    'presmooth_win',0.3,[0,10],...
    'postsmooth_win',0.6,[0,10],...
    'wavelet_flag',false,[false true],...
    'wavelet_dur',1,[0 Inf],...
    'wavelet_peak_width',0.6,[0.01,100],...
    'wavelet_peak_win',0.3,[0,10],...
    'wavelet_edge_width',0.2,[0.01,100],...
    'wavelet_edge_win',0.2,[0,10],...
    'offset_flag',true,[false true],...
    'offset_chanavg_flag',false,[false true],...
    'norm_flag',true,[false true],...
    'norm_chanavg_flag',false,[false true],...
    'norm_postsmooth_flag',false,[false true],...
    'low_band_weight',0.5,[0,1],...
    'spindle_band',[10,16],[0.5,100],...
    'low_band',[4,8],[0.5,100],...
    'high_band',[18,30],[0.5,100],...
    'broad_band',[2,25],[0.5,100],...
    'filt_tfrac',0.3,[0.01,1],...
    'butter_flag',false,[false true],...
    'butter_order',4,[1,10],...
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
    'labels',[],[],...
...
    'verbose',false,[false true],...
    'prep_tags',{'nonspindle_corr_flag' 'power_flag' 'presmooth_win' ...
                 'postsmooth_win' 'wavelet_flag' 'wavelet_dur'...
                 'wavelet_peak_width' 'wavelet_peak_win'...
                 'wavelet_edge_width' 'wavelet_edge_win'...
                 'offset_flag' 'offset_chanavg_flag'...
                 'norm_flag' 'norm_chanavg_flag' 'norm_postsmooth_flag' ...
                 'low_band_weight' 'spindle_band' 'low_band' 'high_band' ...
                 'broad_band' 'filt_tfrac' ...
                 'butter_flag' 'butter_order' ...
                 'labels' ...
                 'verbose'},[],...
    'detect_tags',{'mask' ...
                   'peak_thresh' 'epoch_dur_init' 'edge_thresh'...
                   'relative_edge_thresh_flag' 'relative_edge_thresh'...
                   'dur_thresh' 'exclude_overlap_flag'...
                   'reject_chan_flag'...
                   'reject_low_flag' 'low_thresh'...
                   'reject_high_flag' 'high_thresh'...
                   'reject_band_any_chan_flag' 'reject_band_onset_flag'...
                   'reject_fourier_flag' 'spindle_freq_low'...
                   'spindle_freq_high' 'nonspindle_freq_low'...
                   'nonspindle_freq_high' 'fourier_thresh'...
                   'reject_npeaks_flag' 'min_npeaks'...
                   'reject_peakrate_flag' 'min_peakrate' 'max_peakrate'...
                   'min_peak_deflection' 'min_peak_ratio'...
                   'reject_nchans_flag' 'min_nchans'...
                   'epoch_dur_extract' 'trimmed_epochs_flag'...
                   'extract_reject_flag' 'extract_reject_types' ...
                   'verbose'},[] ...
  });

  % check data matrix
  if numel(size(data))~=2
    error('data matrix must be 2-dimensional');
  end;
  [parms.nchans,parms.ntpoints] = size(data);

  % check mask
  if ~isempty(parms.mask)
    parms.mask = mmil_rowvec(parms.mask);
    if length(parms.mask)~=parms.ntpoints
      error('length of mask must match 2nd dim of data');
    end;
  else
    parms.mask = ones(1,parms.ntpoints);
  end;

  % check extract_reject_types
  if ~iscell(parms.extract_reject_types)
    parms.extract_reject_types = {parms.extract_reject_types};
  end;

  % check chan labels
  if isempty(parms.labels)
    parms.labels = cell(parms.nchans,1);
    for i=1:parms.nchans
      parms.labels{i} = sprintf('chan%d',i);
    end;
  else
    if ~iscell(parms.labels), parms.labels = {parms.labels}; end;
    if length(parms.labels) ~= parms.nchans
      error('labels has wrong number of elements (%d not %d)\n',...
        length(parms.labels),parms.nchans);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
