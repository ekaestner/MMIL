function results = ts_prep_spindle_data(data,sfreq,varargin)
%function results = ts_prep_spindle_data(data,sfreq,[options])
%
% Usage:
%  results = ts_find_spindles(data,sfreq,'key1', value1,...);
%
% Required Input:
%   data: data matrix with size = [nchans,ntpoints]
%   sfreq: sampling frequency (in Hz)
%
% Optional parameters:
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
%
% Optional parameters for offset and normalization:
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
%
% Optional parameters for filtering:
%  'single_flag': [0|1] convert data to single to save space
%    {default = 0}
%  'notch_flag': [0|1] apply notch filters to input data
%    {default = 0}
%  'notch_cf': notch filter cut-out frequency (Hz) (may be vector)
%    {default = 60}
%  'notch_tb': notch filter transition band (Hz) (may be vector)
%    {default = 5}
%  'spindle_band': bounds of spindle frequency range (Hz)
%    {default = [10,16]}
%  'low_band': bounds of low-frequency non-spindle frequency range (Hz)
%    {default = [4,8]}
%  'high_band': bounds of high-frequency non-spindle frequency range (Hz)
%    {default = [18,25]}
%  'broad_band': bounds of broad-band frequency range (Hz)
%    {default = [4,25]}
%  'filt_tfrac': width of filter transition bands relative to cut-off frequency
%    {default = 0.3}
%  'butter_flag': use Butterworth filter instead of FFT-based filter
%    {default = 0}
%  'butter_order': order of Butterworth filter
%    {default = 3}
%
% Optional parameters to place information in output struct
%   'labels': cell array of channel labels
%     if empty will set to channel1, channel2, etc.
%     {default = []}
%   'time': time vector
%     if empty will set based on sfreq and ntpoints ([0:ntpoints-1]/sfreq)
%     {default = []}
%
% Output:
%   results: struct containing these fields:
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
% Created:  01/27/15 by Don Hagler
% Last mod: 04/09/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters
parms = check_input(data,sfreq,varargin);

% intialize results struct
results = init_results(parms);

% convert to single to save space
if parms.single_flag
  data = single(data);
end;

% apply notch filter(s) to data
if parms.notch_flag
  if parms.verbose
    fprintf('%s: notch filtering data...\n',mfilename);
  end;
  results.data_raw = notch_filter_data(parms,data);
else
  results.data_raw = data;
end;
clear data;

% apply narrow band-pass filter to data
if parms.verbose
  fprintf('%s: filtering data with narrow band...\n',mfilename);
end;
results.data_filt_narrow = filter_data(parms,...
                             results.data_raw,parms.spindle_band);

% apply broad band-pass filter to data
if parms.verbose
  fprintf('%s: filtering data with broad band...\n',mfilename);
end;
results.data_filt_broad = filter_data(parms,results.data_raw,parms.broad_band);

% calculate envelope of spindle power (or amplitude)
if parms.verbose
  fprintf('%s: calculating spindle power...\n',mfilename);
end;
if parms.wavelet_flag
  results.spindle_power = convolve_wavelets(parms,...
                            results.data_raw,parms.spindle_band,0);
else
  results.spindle_power = calc_envelope(parms,results.data_filt_narrow);
end;

% calculate low band power (or amplitude)
if parms.verbose
  fprintf('%s: calculating low-band power...\n',mfilename);
end;
if parms.wavelet_flag
  results.low_power = convolve_wavelets(parms,...
                        results.data_raw,parms.low_band,0);
else
  data_filt_low = filter_data(parms,results.data_raw,parms.low_band);
  results.low_power = calc_envelope(parms,data_filt_low);
end;

% calculate high band power (or amplitude)
if parms.verbose
  fprintf('%s: calculating high-band power...\n',mfilename);
end;
if parms.wavelet_flag
  results.high_power = convolve_wavelets(parms,...
                         results.data_raw,parms.high_band,0);
else
  data_filt_high = filter_data(parms,results.data_raw,parms.high_band);
  results.high_power = calc_envelope(parms,data_filt_high);
end;

% adjust spindle_power for nonspindle_power
if parms.nonspindle_corr_flag
  if parms.verbose
    fprintf('%s: adjusting for non-spindle power...\n',mfilename);
  end;
  results.nonspindle_power = zeros(size(results.spindle_power));
  if parms.low_band_weight>0
    results.nonspindle_power = results.nonspindle_power +...
                               parms.low_band_weight * results.low_power;
  end;
  if parms.low_band_weight<1
    results.nonspindle_power = results.nonspindle_power +...
                               (1-parms.low_band_weight) * results.high_power;
  end;
  % adjust spindle_power by nonspindle_power
  results.spindle_power_orig = results.spindle_power;
  if parms.offset_flag
    results.spindle_power = results.spindle_power - results.nonspindle_power;
  else
    results.spindle_power = results.spindle_power ./...
                            (results.nonspindle_power + eps);
  end;
end;

% smooth spindle power
if parms.wavelet_flag
  % convolve with wider wavelets and more smoothing
  if parms.verbose
    fprintf('%s: calculating smooth spindle power...\n',mfilename);
  end;
  results.spindle_power_smooth = ...
    convolve_wavelets(parms,results.data_raw,parms.spindle_band,1);
  if parms.nonspindle_corr_flag
    if parms.verbose
      fprintf('%s: adjusting for non-spindle power...\n',mfilename);
    end;
    results.nonspindle_power_smooth = zeros(size(results.spindle_power_smooth));
    if parms.low_band_weight>0
      results.low_power_smooth = convolve_wavelets(parms,...
                                   results.data_raw,parms.low_band,1);
      results.nonspindle_power_smooth = results.nonspindle_power_smooth +...
                               parms.low_band_weight * results.low_power_smooth;
    end;
    if parms.low_band_weight<1
      results.high_power_smooth = convolve_wavelets(parms,...
                                    results.data_raw,parms.high_band,1);
      results.nonspindle_power_smooth = results.nonspindle_power_smooth +...
                          (1-parms.low_band_weight) * results.high_power_smooth;
    end;
    % adjust spindle_power by nonspindle_power
    results.spindle_power_smooth_orig = results.spindle_power_smooth;
    if parms.offset_flag
      results.spindle_power_smooth =...
        results.spindle_power_smooth - results.nonspindle_power_smooth;
    else
      results.spindle_power_smooth = results.spindle_power_smooth ./...
                                     (results.nonspindle_power_smooth + eps);
    end;
  end;
else
  % post-smooth
  if parms.verbose
    fprintf('%s: smoothing spindle power...\n',mfilename);
  end;
  if parms.postsmooth_win>0
    results.spindle_power_smooth = smooth_data(parms,results.spindle_power,...
                                               parms.postsmooth_win);
    if parms.norm_postsmooth_flag &&...
       (parms.offset_flag || parms.norm_flag)
      results.spindle_power_smooth = ...
        norm_data(parms,results.spindle_power_smooth);
    end;
  else
    results.spindle_power_smooth = results.spindle_power;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(data,sfreq,options)
  parms = mmil_args2parms(options,{...
    'sfreq',sfreq,[0 Inf],...
  ...
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
  ... % offset and normalization
    'offset_flag',true,[false true],...
    'offset_chanavg_flag',false,[false true],...
    'norm_flag',true,[false true],...
    'norm_chanavg_flag',false,[false true],...
    'norm_postsmooth_flag',false,[false true],...
    'low_band_weight',0.5,[0,1],...
  ... % filtering
    'single_flag',true,[false true],...
    'notch_flag',false,[false true],...
    'notch_cf',60,[0.5,1000],...
    'notch_tb',5,[0.01,100],...
    'spindle_band',[10,16],[0.5,100],...
    'low_band',[4,8],[0.5,100],...
    'high_band',[18,25],[0.5,100],...
    'broad_band',[4,25],[0.5,100],...
    'filt_tfrac',0.3,[0.01,1],...
    'butter_flag',false,[false true],...
    'butter_order',4,[1,10],...
  ... % info
    'labels',[],[],...
    'time',[],[],...
  ...
    'verbose',false,[false true],...
  });

  % check data matrix
  if numel(size(data))~=2
    error('data matrix must be 2-dimensional');
  end;
  [parms.nchans,parms.ntpoints] = size(data);

  % check filtering bands
  band_names = {'spindle_band','low_band','high_band','broad_band'};
  for i=1:length(band_names)
    band_name = band_names{i};
    if numel(parms.(band_name))~=2
      error('%s must have 2 elements',band_name);
    end;
  end;  

  % check notch filter settings
  if parms.notch_flag
    if length(parms.notch_cf)>1 && length(parms.notch_tb)==1
      parms.notch_tb = parms.notch_tb * ones(size(parms.notch_cf));
    end;
    if length(parms.notch_cf) ~= length(parms.notch_tb)
      error('notch_cf and notch_tb have mismatched elements');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  results.sfreq = parms.sfreq;
  if isempty(parms.labels)
    results.labels = cell(parms.nchans,1);
    for i=1:parms.nchans
      results.labels{i} = sprintf('channel%d',i);
    end;
  else
    results.labels = parms.labels;
  end;
  if isempty(parms.time)
    results.time = [0:parms.ntpoints-1]/parms.sfreq;
  else
    results.time = parms.time;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = notch_filter_data(parms,data)
  for i=1:length(parms.notch_cf)
    cf = parms.notch_cf(i);
    if parms.butter_flag
      [b,a] = butter(parms.butter_order,cf/(parms.sfreq/2),'stop');
      data = (filtfilt(b,a,data'))';
    else
      tb = parms.notch_tb(i);
      data = ts_freq_filt(data',parms.sfreq,cf,tb,'notch',1)';
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = filter_data(parms,data,cf_vec)
  tb_vec = parms.filt_tfrac*cf_vec;
  for i=1:parms.nchans
    tmp_data = data(i,:)';
    % remove mean in case there is a problem with detrending
    tmp_data = detrend(tmp_data,'constant');
    % try to detrend (may fail in some cases)
    warning('off');
    tmp_data = detrend(tmp_data);
    warning('on');
    if parms.butter_flag
      [b,a] = butter(parms.butter_order,cf_vec/(parms.sfreq/2));
      tmp_data = filtfilt(b,a,double(tmp_data));
    else
      tmp_data = ts_freq_filt(tmp_data,...
        parms.sfreq,cf_vec,tb_vec,'bandpass')';
    end;
    data(i,:) = tmp_data';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = convolve_wavelets(parms,data,freq_band,peak_flag)
  % convolve data
  tparms = [];
  tparms.wavelet_freqs = [freq_band(1):freq_band(2)];
  tparms.wavelet_dur = parms.wavelet_dur;
  if peak_flag
    tparms.wavelet_width = parms.wavelet_peak_width;
    tparms.wavelet_win = parms.wavelet_peak_win;
  else
    tparms.wavelet_width = parms.wavelet_edge_width;
    tparms.wavelet_win = parms.wavelet_edge_win;
  end;
  tparms.offset_flag = 0;
  tparms.norm_flag = 0;
  args = mmil_parms2args(tparms);
  data = ts_convolve_spindles(data,parms.sfreq,args{:});
  % square result
  if parms.power_flag
    data = data.^2;
  end;
  % offset and normalize  
  if parms.offset_flag || parms.norm_flag
    data = norm_data(parms,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = calc_envelope(parms,data)
  if parms.power_flag
    data = data.^2;
  else
    data = abs(data);
  end;
  data = smooth_data(parms,data,parms.presmooth_win);
  if parms.offset_flag || parms.norm_flag
    data = norm_data(parms,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = norm_data(parms,data)
  % calculate median absolute deviation and median
  [wmad,wmed] =  mmil_wtd_mad(data,[],2);
  % subtract median
  if parms.offset_flag
    if parms.offset_chanavg_flag
      data = data - mean(wmed);
    else
      data = bsxfun(@minus,data,wmed);
    end;
  end;
  % normalize by mad
  if parms.norm_flag
    if parms.norm_chanavg_flag
      data = data / mean(wmad);
    else
      data = bsxfun(@rdivide,data,wmad);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = smooth_data(parms,data,win)
  if win==0, return; end;
  ntpoints = size(data,2);
  nwin = round(win*parms.sfreq);
  avg_window = tukeywin(nwin)'/nwin;
  for i=1:parms.nchans
    data_chan = data(i,:);
    data_tmp = zeros(1,ntpoints+2*nwin);
    data_tmp(nwin+1:nwin+ntpoints) = data_chan;
    data_tmp = conv(data_tmp,avg_window,'same');
    data_chan = data_tmp(nwin+1:nwin+ntpoints);
    data(i,:) = data_chan;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

