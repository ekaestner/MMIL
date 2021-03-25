function avg_data = ts_postprocess_avg(avg_data,varargin);
% avg_data_post = ts_postprocess_avg(avg_data,[options]);
%
% Usage:
%  avg_data_post = ts_postprocess_avg(avg_data, 'key1', value1,...);
%
% Required input:
%  avg_data - averaged data structure (output from avg_fif_data.m)
%
% Optional parameters:
%  'stim_delay': offset applied to time zero (msec)
%    (will be subtracted from time vector)
%    can be used to adjust for stimulus onset delay after trigger
%    {default = 0}
%  'bandpass_flag': [0|1] bandpass frequency filter
%    {default = 0}
%  'bandpass_low_cf': low cutoff frequency (high-pass filter) (Hz)
%    {default = 0}
%  'bandpass_low_tb': low cutoff transition band (Hz)
%    {default = 0}
%  'bandpass_high_cf': high cutoff frequency (low-pass filter) (Hz)
%    {default = 100}
%  'bandpass_high_tb': high cutoff transition band (Hz)
%    {default = 0}
%  'notch_flag': [1|0] notch fft filter before averaging
%    {default = 0}
%  'notch_cf': notch center frequency (notch filter) (Hz)
%    {default = 0}
%  'notch_tb': notch transition band (Hz)
%    {default = 0}
%  'dsfact': downsampling factor -- must be an integer
%    data is downsampled to lower sampling frequency
%    e.g. original sampling freq = 1000, dsfact = 4,
%        resulting sampling freq = 250
%    {default = 1} (no downsampling)
%  'detrend_flag': [0|1] linear trend removal
%    {default = 1}
%  'baseline_flag': [0|1] subtract mean of baseline period
%    {default = 1}
%  'baseline_start': start time of baseline period (msec)
%    relative to time zero
%    {default = -Inf}
%  'baseline_end': end time of baseline period (msec)
%    relative to time zero
%    {default = 0}
%  'badchans': vector of bad channel numbers
%    {default = []}
%  'badchanfile': name of text file containing bad channel labels
%    {default = []}
%  'rm_badchans_flag': [0|1] complete removal of badchans
%    from data matrices -- otherwise just set to zero
%    {default = 0}
%
% Output:
%   avg_data - structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        typestring   (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
%        type         (int)
%        kind         (int)
%        lognum       (int)
%      coor_trans     (struct)
%        device2head  (4x4 matrix of doubles) (or empty)
%        mri2head     (4x4 matrix of doubles) (or empty)
%      averages       (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples)
%      noise          (struct)
%        num_trials   (int)
%        scale_fact   (double)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%
% Note on order of operations:
%  When specified (i.e. bandpass=1), filtering is performed first.
%     If filtering is done, detrending is done first regardless of whether
%     detrend_flag=1.
%  Downsampling is then performed (if dsfact>1),
%    followed by detrending (if detrend_flag=1).
%    and finally baseline correction (if baseline_flag=1),
%
%  Created:  04/26/06  by Don Hagler
%  Last Mod: 12/29/13  by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'stim_delay',0,[],...
  'bandpass_flag',false,[false,true],...
  'bandpass_low_cf',0,[0,Inf],...
  'bandpass_low_tb',0,[0,Inf],...
  'bandpass_high_cf',100,[0,Inf],...
  'bandpass_high_tb',0,[0,Inf],...
  'notch_flag',false,[false,true],...
  'notch_cf',0,[0,Inf],...
  'notch_tb',0,[0,Inf],...
  'dsfact',1,[1,Inf],...
  'detrend_flag',true,[false,true],...
  'baseline_flag',true,[false,true],...
  'baseline_start',-Inf,[-Inf,Inf],...
  'baseline_end',0,[-Inf,Inf],...
  'badchans',[],[],...
  'badchanfile',[],[],...
  'rm_badchans_flag',false,[false true],... 
...
  'datatype','single',{'single','double'},...
  'pad_nsamples',1000,[0,Inf],...
  'verbose',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reject bad parameters
if ~mmil_isint(parms.dsfact)
  error('dsfact must be an integer >= 1');
end;
if parms.bandpass_flag
  if parms.bandpass_low_cf > parms.bandpass_high_cf
    error('low cutoff freq must be less than high cutoff freq');
  end;
end;

% read badchan file
labels = {avg_data.sensor_info.label}; 
if ~isempty(parms.badchanfile)
  badchan_i = ts_read_txt_badchans(parms.badchanfile,labels);
else
  badchan_i = [];
end;
parms.badchans = unique([parms.badchans,badchan_i]);

% check bad chans are in bounds
if ~isempty(parms.badchans)
  parms.badchans = parms.badchans(...
    union(find(parms.badchans>0),find(parms.badchans<=avg_data.num_sensors)));
end;
% make sure there are some good channels left
goodchans = setdiff([1:avg_data.num_sensors],parms.badchans);
if isempty(goodchans)
  error('no good channels specified');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nconds = length(avg_data.averages);

% remove or zero bad chans
if ~isempty(parms.badchans)
  if parms.rm_badchans_flag
    if parms.verbose
      fprintf('%s: removing badchans...\n',mfilename);
    end;
    avg_data.sensor_info = avg_data.sensor_info(goodchans);
    avg_data.num_sensors = length(goodchans);
    for j=1:nconds
      avg_data.averages(j).data = avg_data.averages(j).data(goodchans,:);
    end;
    avg_data.noise.covar = avg_data.noise.covar(goodchans,:);
    avg_data.noise.covar = avg_data.noise.covar(:,goodchans);
  else
    fprintf('%s: zeroing badchans...\n',mfilename);
    for b=1:length(parms.badchans)
      k = parms.badchans(b);
      avg_data.sensor_info(k).badchan = 1;
    end;
    for j=1:nconds
      avg_data.averages(j).data(parms.badchans,:) = 0;
    end;
    avg_data.noise.covar(parms.badchans,:) = 0;
    avg_data.noise.covar(:,parms.badchans) = 0;
  end;
end;

% filter
if parms.bandpass_flag || parms.notch_flag
  if parms.verbose
    fprintf('%s: applying bandpass filter...\n',mfilename);
  end;
  for j=1:nconds
    matsize = size(avg_data.averages(j).data);
    % remove linear trend
    avg_data.averages(j).data = detrend(avg_data.averages(j).data')';
    % zero pad data to avoid filter artifacts
    num_samples = matsize(2) + 2*parms.pad_nsamples;
    tmp_data = zeros(matsize(1),num_samples,parms.datatype);
    tmp_data(:,parms.pad_nsamples+1:num_samples-parms.pad_nsamples) =...
      avg_data.averages(j).data;
    if parms.bandpass_flag && parms.notch_flag
      tmp_data = ts_freq_filt(tmp_data',avg_data.sfreq,...
          [parms.bandpass_low_cf,parms.bandpass_high_cf,parms.notch_cf],...
          [parms.bandpass_low_tb,parms.bandpass_high_tb,parms.notch_tb],...
          'bandpassnotch')';
    elseif parms.bandpass_flag
        tmp_data = ts_freq_filt(tmp_data',avg_data.sfreq,...
          [parms.bandpass_low_cf,parms.bandpass_high_cf],...
          [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass')';
    elseif parms.notch_flag
        tmp_data = ts_freq_filt(tmp_data',avg_data.sfreq,...
          parms.notch_cf,parms.notch_tb,...
          'notch')';
    end;
    if isempty(tmp_data)
      error('error filtering data');
    end;
    % unpad data
    avg_data.averages(j).data = ...
      tmp_data(:,parms.pad_nsamples+1:num_samples-parms.pad_nsamples);
  end;
end

% adjust times with stim delay
if parms.stim_delay
  for j=1:nconds
    avg_data.averages(j).time = ...
      avg_data.averages(j).time - parms.stim_delay/1000;
  end;
end

% downsample
if parms.dsfact~=1
  if parms.verbose
    fprintf('%s: downsampling...\n',mfilename);
  end;
  for j=1:nconds
    avg_data.averages(j).data = resample(double(avg_data.averages(j).data'),1,parms.dsfact,0)';
    avg_data.averages(j).stdev = resample(double(avg_data.averages(j).stdev'),1,parms.dsfact,0)';
    avg_data.averages(j).time = resample(double(avg_data.averages(j).time),1,parms.dsfact,0);
    if strcmp(parms.datatype,'single')
      avg_data.averages(j).data = single(avg_data.averages(j).data);
      avg_data.averages(j).stdev = single(avg_data.averages(j).stdev);
      avg_data.averages(j).time = single(avg_data.averages(j).time);
    end;
  end;
  avg_data.sfreq = avg_data.sfreq/parms.dsfact;
end;

% linear trend removal
if parms.detrend_flag
  if parms.verbose
    fprintf('%s: removing linear trend...\n',mfilename);
  end;
  for j=1:nconds
    avg_data.averages(j).data = detrend(avg_data.averages(j).data')';
  end;
end;

% subtract baseline
if parms.baseline_flag
  if parms.verbose
    fprintf('%s: subtracting baseline...\n',mfilename);
  end;
  for j=1:nconds
    % could assume that all conditions have the same timing,
    % but why not be flexible?
    stim_delay = -avg_data.averages(j).time(1)*1000;
    num_samples = length(avg_data.averages(1).time);
    baseline_start_samp = round((parms.baseline_start + stim_delay) * avg_data.sfreq/1000);
    baseline_end_samp   = round((parms.baseline_end + stim_delay) * avg_data.sfreq/1000);
    baseline_start_samp = min(max(1,baseline_start_samp),num_samples);
    baseline_end_samp = min(max(1,baseline_end_samp),num_samples);
    mean_baseline = ...
      mean(avg_data.averages(j).data(:,baseline_start_samp:baseline_end_samp),2);
    avg_data.averages(j).data = ...
      avg_data.averages(j).data - mean_baseline*ones(1,num_samples);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.verbose
  fprintf('%s: finished\n',mfilename);
end;

