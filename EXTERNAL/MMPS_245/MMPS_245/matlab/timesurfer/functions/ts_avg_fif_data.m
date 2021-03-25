function output_data = ts_avg_fif_data(datafile,varargin);
% output_data = ts_avg_fif_data(datafile,[options]);
%
% Usage:
%  output_data = ts_avg_fif_data(datafile, 'key1', value1,...);
%
% Required input:
%  datafile - full or relative path name of Neuromag fif raw data file
%
% Optional parameters:
%  'evnts': [struct] event structure with the following fields:
%     evnts.type      = string (e.g. trigger, skip, reject, manual)
%     evnts.latency   = expressed in samples, first sample of file is 1
%     evnts.condition = numeric event code
%     evnts.duration  = expressed in samples
%  'hdr': [struct] hdr structure with the following fields:
%     hdr.sfreq           sampling frequency
%     hdr.nChans          number of channels
%     hdr.nSamples        number of samples per buffer
%     hdr.tfirst          first time point in file
%     hdr.tlast           last time point in file
%     hdr.nBuffs          number of buffers in file (not including skips)
%     hdr.sensor          structure containing sensor information
%
%   N.B. if evnts or hdr structures are not supplied, will read events
%     and header info from raw data file (this can be time consuming)
%     see ts_read_fif_events and ts_read_fif_header
%
%  'prestim_dur':  duration of prestimulus period (msec)
%     { default: 100 }
%  'poststim_dur': duration of poststimulus period (msec)
%     { default: 400 }
%  'stim_delay': duration of stimulus onset delay after trigger (msec)
%     { default: 0 }
%  'badchans': vector of bad channel indices
%     { default: [] }
%  'badchanfile': name of text file containing bad channel labels
%    {default: []}
%  'bandpass_flag': [1|0] Toggle bandpass frequency filter
%     { default: 0 }
%  'bandpass_low_cf': low cutoff frequency (high-pass filter) (Hz)
%     { default: 0 }
%  'bandpass_low_tb': low cutoff transition band (Hz)
%     { default: 0 }
%  'bandpass_high_cf': high cutoff frequency (low-pass filter) (Hz)
%     { default: 100 }
%  'bandpass_high_tb': high cutoff transition band (Hz)
%     { default: 0 }
%  'notch_flag': [1|0] Toggle notch fft filter before averaging
%     {default: 0}
%  'notch_cf':  notch center frequency (notch filter) (Hz)
%     {default: 0}
%  'notch_tb':  notch transition band (Hz)
%     {default: 0}
%  'dsfact':  downsampling factor -- must be an integer
%     data is downsampled to lower sampling frequency
%     e.g. original sampling freq = 1000, dsfact = 4,
%         resulting sampling freq = 250
%     { default: 1 (no downsampling) }
%  'detrend_flag': [1|0] Toggle linear trend removal
%    { default: 1 }
%  'baseline_flag': [1|0] Toggle subtract mean of baseline period
%    { default: 1 }
%  'baseline_start': start time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     { default: -Inf } (start at beginning of prestimulus period)
%  'baseline_end'   - end time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     { default: 0 } (end at trigger onset)
%  'ncov_ex_evnts': vector of event codes that should not
%     be used in calculating the noise covariance matrix
%     { default: [] }
%  'reject_mag':  automatic rejection threshold for magnetometer channels (fT)
%     if 0, rejection based on magnetometers is disabled
%     { default: 6000 }
%  'reject_grad': automatic rejection threshold for gradiometer channels (fT/cm)
%     if 0, rejection based on gradiometers is disabled
%     { default: 3000 }
%  'reject_eeg': automatic rejection threshold for eeg channels (uV)
%     if 0, rejection based on eeg is disabled
%     { default: 0 }
%  'reject_eog': automatic rejection threshold for eog channel (uV)
%     if 0, rejection based on eog is disabled
%     { default: 200 }
%  'readtrans_flag':  [1|0] Toggle read device2head transform from fif file
%     if not found in fif file, will cause core dump
%     { default: 1 }
%  'max_num_trials': maximum number of trials per condition
%    {default: Inf} (infinite)
%  'save_epochs_flag': [1|0] instead of calculating averages, return structure
%     containing epochs (single trial time courses)
%     { default: 0 }
%
% Output:
%   avg_data - structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        type         (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
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
%        stdev        (matrix of doubles) (num_sensors x num samples)
%      noise          (struct)
%        num_trials   (int)
%        num_samples  (int)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%   epoch_data - structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        type         (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
%      coor_trans     (struct)
%        device2head  (4x4 matrix of doubles) (or empty)
%        mri2head     (4x4 matrix of doubles) (or empty)
%      epochs         (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples x num_trials)
%      noise          (struct)
%        num_trials   (int)
%        num_samples  (int)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%
% Note on order of operations:
%  When specified (e.g. bandpass_flag=1), filtering is performed first.
%     If filtering is done, detrending is done first regardless of whether
%     detrend_flag=1.
%  Downsampling is then performed (if dsfact>1),
%    followed by detrending (if detrend_flag=1).
%    and finally baseline correction (if baseline_flag=1),
%  Automatic artifact rejection is then done, followed by averaging.
%
%  created:  04/18/06   by Don Hagler
%  last mod: 11/04/14   by Don Hagler
%

%% todo: replace use of fiff access toolbox with MNE matlab toolbox
%%  e.g. use ts_MNE_loadfif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'evnts',[],[],...
  'hdr',[],[],...
  'prestim_dur',100,[0,Inf],...
  'poststim_dur',400,[0,Inf],...
  'stim_delay',0,[],...
  'badchans',[],[],...
  'detrend_flag',true,[false,true],...
  'baseline_flag',true,[false,true],...
  'baseline_start',-Inf,[-Inf,Inf],...
  'baseline_end',0,[-Inf,Inf],...
  'reject_mag',6000,[],...
  'reject_grad',3000,[],...
  'reject_eeg',0,[],...
  'reject_eog',200,[],...
  'readtrans_flag',true,[false,true],...
  'bandpass_flag',false,[false,true],...
  'bandpass_low_cf',0,[0,Inf],...
  'bandpass_low_tb',0,[0,Inf],...
  'bandpass_high_cf',100,[0,Inf],...
  'bandpass_high_tb',0,[0,Inf],...
  'notch_flag',false,[false,true],...
  'notch_cf',0,[0,Inf],...
  'notch_tb',0,[0,Inf],...
  'dsfact',1,[1,Inf],...
  'ncov_ex_evnts',[],[],...
  'max_num_trials',Inf,[],...
  'save_epochs_flag',false,[false,true],...
  'badchanfile',[],[],...
...
  'prescale_mag',10^15,[],...  % to fT
  'prescale_grad',10^13,[],... % to fT/cm
  'prescale_eeg',10^6,[],...   % to uV
  'prescale_eog',10^6,[],...   % to uV
  'datatype','single',{'single','double'},...
  'max_noise_samples',5000,[],...
  'concat_ncov_flag',false,[false true],...% concatenate noise before calculating noise covar
});                                        % or calculate noise covar for each trial and then avg

output_data = [];

if(~exist(datafile,'file'))
  error('datafile %s not found',datafile);
end

SKIP_LENGTH = 100; % do not change -- must be consistent with ts_read_fif_data

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

% optionally read events and header (if not supplied on command line)
if isempty(parms.evnts)
  % read events
  fprintf('%s: reading events for %s...\n',mfilename,datafile);
  [parms.hdr,parms.evnts]=ts_read_fif_events(datafile);
elseif isempty(parms.hdr)
  % read header
  fprintf('%s: reading header for %s...\n',mfilename,datafile);
  parms.hdr = ts_read_fif_header(datafile,1);
end

% read badchan file
labels = parms.hdr.sensors.label; 
if ~isempty(parms.badchanfile)
  badchan_i = ts_read_txt_badchans(parms.badchanfile,labels);
else
  badchan_i = [];
end;
parms.badchans = unique([parms.badchans,badchan_i]);

% check bad chans are in bounds
if ~isempty(parms.badchans)
  parms.badchans = parms.badchans(...
    union(find(parms.badchans>0),find(parms.badchans<=parms.hdr.nChans)));
end;
% make sure there are some good channels left
goodchans = setdiff([1:parms.hdr.nChans],parms.badchans);
if isempty(goodchans)
  error('no good channels specified');
end;

% convert from msec to samples and introduce stim_delay
prestim_samp  = round((parms.prestim_dur - parms.stim_delay) * parms.hdr.sfreq/1000);
poststim_samp = round((parms.poststim_dur + parms.stim_delay) * parms.hdr.sfreq/1000);
time_orig = [-prestim_samp:poststim_samp]/parms.hdr.sfreq - parms.stim_delay/1000;

% account for downsampling
sfreq_ds = parms.hdr.sfreq/parms.dsfact;
time = downsample(time_orig,parms.dsfact);
num_samples = length(time);
if strcmp(parms.datatype,'single')
  time = single(time);
end;

% find sample closest in time to that desired for baseline and noise
[tmp,baseline_start_samp] = min(abs(time-parms.baseline_start/1000));
[tmp,baseline_end_samp] = min(abs(time-parms.baseline_end/1000));

% make sure events have proper fields
if isempty(parms.evnts)
  error('evnts structure is empty');
end;
if ~isfield(parms.evnts, 'type') | ~isfield(parms.evnts, 'latency') | ...
   ~isfield(parms.evnts, 'duration') | ~isfield(parms.evnts, 'condition')
  error('events structure does not have correct fields');
end;
% sort evnts by latency
[slat, indx] = sort([parms.evnts.latency]);
parms.evnts = parms.evnts(indx);
% separate evnts by type
skip_i = find(strcmp('skip',{parms.evnts.type}));
reject_i = find(strcmp('reject',{parms.evnts.type}));
other_i = setdiff(1:length(parms.evnts),union(skip_i,reject_i));
skip_evnts = parms.evnts(skip_i);
reject_evnts = parms.evnts(reject_i);
parms.evnts = parms.evnts(other_i);
skips_beg   = cell2mat({skip_evnts.latency});
skips_end   = skips_beg + cell2mat({skip_evnts.duration}) - 1;
rejects_beg = cell2mat({reject_evnts.latency});
rejects_end = rejects_beg + cell2mat({reject_evnts.duration}) - 1;
evnts_lats  = cell2mat({parms.evnts.latency});
% make sure event conditions are numeric
for evnum=1:length(parms.evnts)
  if ~isnumeric(parms.evnts(evnum).condition) | size(parms.evnts(evnum).condition)~=1
    error('non-numeric event code found for event %d',evnum);
  end;
end;
% make list of unique event codes
event_codes = unique(cell2mat({parms.evnts.condition}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output_data structure;
%output_data.subject_id = ''; % let user set this or somehow get from fif?
%output_data.date = ''; % how to get this from fif?
output_data.num_sensors = parms.hdr.nChans;
output_data.sfreq = sfreq_ds;
output_data.sensor_info = [];
for k = 1:parms.hdr.nChans
  output_data.sensor_info(k).typestring = parms.hdr.sensors.typestring{k};
  output_data.sensor_info(k).label = parms.hdr.sensors.label{k};
  output_data.sensor_info(k).loc = parms.hdr.sensors.loc{k};
  output_data.sensor_info(k).badchan = 0;
  output_data.sensor_info(k).type = parms.hdr.sensors.type(k);
  output_data.sensor_info(k).kind = parms.hdr.sensors.kind(k);
  output_data.sensor_info(k).lognum = parms.hdr.sensors.lognum(k);
end;
for b = 1:length(parms.badchans)
  k = parms.badchans(b);
  output_data.sensor_info(k).badchan = 1;
end;  
output_data.coor_trans = [];
if parms.readtrans_flag
  output_data.coor_trans.device2head = loadtrans(datafile); % Uutela's fiff access toolbox meg_pd
else
  output_data.coor_trans.device2head = [];
end;
output_data.coor_trans.mri2head = [];
% initialize averages/epochs
if parms.save_epochs_flag
  for j=1:length(event_codes),
    output_data.epochs(j).event_code=event_codes(j);
    output_data.epochs(j).num_trials=0;
    output_data.epochs(j).num_rejects.mag=0;
    output_data.epochs(j).num_rejects.grad=0;
    output_data.epochs(j).num_rejects.eeg=0;
    output_data.epochs(j).num_rejects.eog=0;
    output_data.epochs(j).num_rejects.manual=0;
    output_data.epochs(j).num_rejects.skip=0;
    output_data.epochs(j).time=time;
    output_data.epochs(j).data=[];
  end
else
  data = zeros(parms.hdr.nChans,num_samples,parms.datatype);
  for j=1:length(event_codes),
    output_data.averages(j).event_code=event_codes(j);
    output_data.averages(j).num_trials=0;
    output_data.averages(j).num_rejects.mag=0;
    output_data.averages(j).num_rejects.grad=0;
    output_data.averages(j).num_rejects.eeg=0;
    output_data.averages(j).num_rejects.eog=0;
    output_data.averages(j).num_rejects.manual=0;
    output_data.averages(j).num_rejects.skip=0;
    output_data.averages(j).time=time;
    output_data.averages(j).data=data;
    output_data.averages(j).stdev=data;
  end
end;
output_data.noise = [];
output_data.noise.num_trials = 0;
output_data.noise.num_samples = 0;
output_data.noise.covar = zeros(parms.hdr.nChans,parms.hdr.nChans,parms.datatype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine which channels are which type
mag_i = find(strcmp('mag',{parms.hdr.sensors.typestring{:}}));
grad_i = find(strncmp('grad',{parms.hdr.sensors.typestring{:}},4));
eeg_i = find(strcmp('eeg',{parms.hdr.sensors.typestring{:}}));
eog_i = find(strcmp('eog',{parms.hdr.sensors.typestring{:}}));
% generate rejection thresholds for each channel
reject_thresh = zeros(parms.hdr.nChans,1,parms.datatype);
reject_thresh(mag_i)  = parms.reject_mag/parms.prescale_mag;
reject_thresh(grad_i) = parms.reject_grad/parms.prescale_grad;
reject_thresh(eeg_i)  = parms.reject_eeg/parms.prescale_eeg;
reject_thresh(eog_i)  = parms.reject_eog/parms.prescale_eog;
reject_thresh(find(reject_thresh<=0)) = inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read data
fprintf('%s: reading data from file %s...\n',mfilename,datafile);
noise_matrix = [];
evnum = 0;
while evnum < length(parms.evnts)
  evnum = evnum + 1;
  j = find(event_codes==parms.evnts(evnum).condition);
  begsamp = parms.evnts(evnum).latency - prestim_samp;
  endsamp = parms.evnts(evnum).latency + poststim_samp;

  % read data around the stimulus epoch to reduce ringing from the filter
  num_stimsamps = endsamp - begsamp + 1;
  if parms.bandpass_flag
    read_begsamp = begsamp - num_stimsamps;
    read_endsamp = endsamp + num_stimsamps;
  else
    read_begsamp = begsamp;
    read_endsamp = endsamp;
  end;

  if(read_begsamp < 1)
    continue;
  end;
  if(read_endsamp >= parms.hdr.tlast*parms.hdr.sfreq)
    fprintf('%s: reached end of file\n',mfilename);
    break;
  end;

  % read to end of buffer that read_endsamp is in
  % must treat skips and data buffers differently
  skip_i = find(parms.hdr.skips<=read_endsamp);
  num_samps = read_endsamp - length(skip_i)*SKIP_LENGTH;
  % calculate number of regular (non-skip) buffers
  num_buffs = ceil(num_samps/parms.hdr.nSamples);
  % round up to end of buffer
  num_samps = num_buffs*parms.hdr.nSamples;
  % add back in skips
  num_samps = num_samps + length(skip_i)*SKIP_LENGTH;
  read_endsamp = num_samps;

  % check for skip
  skip_i = find((skips_beg >= read_begsamp & skips_beg <= read_endsamp) | ...
                (skips_end >= read_begsamp & skips_end <= read_endsamp));
  if ~isempty(skip_i)
    if parms.save_epochs_flag
      output_data.epochs(j).num_rejects.skip=...
        output_data.epochs(j).num_rejects.skip + 1;
    else
      output_data.averages(j).num_rejects.skip=...
        output_data.averages(j).num_rejects.skip + 1;
    end;
    continue;
  end;

  % check for manual reject
  reject_i = find((rejects_beg >= begsamp & rejects_beg <= endsamp) | ...
                (rejects_end >= begsamp & rejects_end <= endsamp));
  if ~isempty(reject_i)
    if parms.save_epochs_flag
      output_data.epochs(j).num_rejects.manual=...
        output_data.epochs(j).num_rejects.manual + 1;
    else
      output_data.averages(j).num_rejects.manual=...
        output_data.averages(j).num_rejects.manual + 1;
    end;
    continue;
  end;

  % read selected samples from file
  [data,eofstatus] = ...
    ts_read_fif_data(datafile,parms.hdr,read_begsamp,read_endsamp);
  if(eofstatus)
    fprintf('%s: reached end of file\n',mfilename);
    break;
  end

  % make tmpevnts containing events within read data
  % (purpose is to reduce unnecessary rereading of buffers)
  tmp_evnts_i = find((evnts_lats > read_begsamp & evnts_lats < read_endsamp));
  tmp_evnts = [];
  for tmp_evnum=1:length(tmp_evnts_i)
    reject_event_flag = 0;
    begsamp = parms.evnts(tmp_evnts_i(tmp_evnum)).latency - prestim_samp;
    endsamp = parms.evnts(tmp_evnts_i(tmp_evnum)).latency + poststim_samp;
    if parms.bandpass_flag
      % make sure enough time will be read to get pre- and post-read buffers
      tmp_read_begsamp = begsamp - num_stimsamps;
      tmp_read_endsamp = endsamp + num_stimsamps;
    else
      tmp_read_begsamp = begsamp;
      tmp_read_endsamp = endsamp;
    end;
    if(tmp_read_begsamp >= read_begsamp & tmp_read_endsamp <= read_endsamp)
      tmp_evnts = [tmp_evnts parms.evnts(tmp_evnts_i(tmp_evnum))];
    end;
  end;

  % bandpass and/or notch filter
  if parms.bandpass_flag || parms.notch_flag
    data = detrend(data')';
    if parms.bandpass_flag && parms.notch_flag
      data = ts_freq_filt(data',parms.hdr.sfreq,...
        [parms.bandpass_low_cf,parms.bandpass_high_cf,parms.notch_cf],...
        [parms.bandpass_low_tb,parms.bandpass_high_tb,parms.notch_tb],...
        'bandpassnotch')';
    elseif parms.bandpass_flag
      data = ts_freq_filt(data',parms.hdr.sfreq,...
        [parms.bandpass_low_cf,parms.bandpass_high_cf],...
        [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass')';
    elseif parms.notch_flag
      data = ts_freq_filt(data',parms.hdr.sfreq,...
        parms.notch_cf,parms.notch_tb,...
        'notch')';
    end;
    if isempty(data)
      error('error filtering data');
    end;
  end;

  % zero bad channels
  data(parms.badchans,:) = 0;

  % cycle through all tmpevnts
  for tmp_evnum=1:length(tmp_evnts)
    j = find(event_codes==tmp_evnts(tmp_evnum).condition);
    begsamp = tmp_evnts(tmp_evnum).latency - prestim_samp;
    endsamp = tmp_evnts(tmp_evnum).latency + poststim_samp;

    if parms.save_epochs_flag
      if output_data.epochs(j).num_trials>=parms.max_num_trials
        break;
      end;
    elseif output_data.averages(j).num_trials>=parms.max_num_trials
      break;
    end;

    if tmp_evnum>1
      % check for manual reject
      reject_i = find((rejects_beg >= begsamp & rejects_beg <= endsamp) | ...
                    (rejects_end >= begsamp & rejects_end <= endsamp));
      if ~isempty(reject_i)
        reject_event_flag = 1;
        if parms.save_epochs_flag
          output_data.epochs(j).num_rejects.manual=...
            output_data.epochs(j).num_rejects.manual + 1;
        else
          output_data.averages(j).num_rejects.manual=...
            output_data.averages(j).num_rejects.manual + 1;
        end;
        continue;
      end;
    end;

    % resize to desired stimulus epoch
    endsamp = endsamp - read_begsamp + 1;
    begsamp = begsamp - read_begsamp + 1;
    tmp_data = data(:,begsamp:endsamp);

    % downsample
    if parms.dsfact~=1
      tmp_data = downsample(tmp_data',parms.dsfact)';
      tmp_nsamp = size(tmp_data,2);
      if tmp_nsamp<num_samples
        tmp_data(:,tmp_nsamp+1:num_samples)=...
          zeros(1,num_samples-tmp_nsamp,parms.datatype);
      end;
    end;

    % linear trend removal
    if parms.detrend_flag
      tmp_data = detrend(tmp_data')';
    end;

    % baseline correction
    if parms.baseline_flag
      mean_baseline = mean(tmp_data(:,baseline_start_samp:baseline_end_samp),2);
      tmp_data = tmp_data - mean_baseline*ones(1,num_samples);
    end;

    % automatic artifact rejection
    [reject_chans,reject_samples] = ...
      find(abs(tmp_data)>...
        reject_thresh*ones(1,num_samples));

    % classify rejects
    if ~isempty(reject_chans)
      reject_event_flag = 1;
      k=reject_chans(1);
      type = regexprep(parms.hdr.sensors.typestring{k},'grad\d','grad');
      if parms.save_epochs_flag
        output_data.epochs(j).num_rejects.(type)=...
          output_data.epochs(j).num_rejects.(type) + 1;
      else
        output_data.averages(j).num_rejects.(type)=...
          output_data.averages(j).num_rejects.(type) + 1;
      end;
    end;

    if ~reject_event_flag % add to data matrices only if not rejected
      if parms.save_epochs_flag
        % add data to epoch struct
        output_data.epochs(j).num_trials = output_data.epochs(j).num_trials + 1;
        output_data.epochs(j).data(:,:,output_data.epochs(j).num_trials) = tmp_data;
      else
        % add data to average struct
        output_data.averages(j).data = output_data.averages(j).data + tmp_data;
        output_data.averages(j).stdev = output_data.averages(j).stdev + tmp_data.*tmp_data;
        output_data.averages(j).num_trials = output_data.averages(j).num_trials + 1;
      end;

      % add data to noise matrix
      ex=find(parms.ncov_ex_evnts==tmp_evnts(tmp_evnum).condition);
      if isempty(ex)
        noise = tmp_data(:,baseline_start_samp:baseline_end_samp);
        nsamps = size(noise,2);
        noise = noise - mean(noise,2)*ones(1,nsamps);
        if parms.concat_ncov_flag
          if size(noise_matrix,2) < parms.max_noise_samples
            noise_matrix = [noise_matrix noise];
            output_data.noise.num_trials = output_data.noise.num_trials + 1;
            output_data.noise.num_samples = output_data.noise.num_samples + nsamps;
          end;
        else
          output_data.noise.covar = output_data.noise.covar + cov(noise');
          output_data.noise.num_trials = output_data.noise.num_trials + 1;
          output_data.noise.num_samples = output_data.noise.num_samples + nsamps;
        end;
      end;
    end;    
    
    if tmp_evnum>1 % more than one event in tmp_evnts
      % add one to evnum to keep track that we are processing this one now
      evnum = evnum+1; 
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~parms.save_epochs_flag
  % calculate averages from sums
  for j=1:length(event_codes)
    N = output_data.averages(j).num_trials;
    if (N>0)
      % calculate average
      output_data.averages(j).data = ...
        output_data.averages(j).data / N;
    end;
    if (N>1)
      % calculate standard deviation
      output_data.averages(j).stdev = ...
        sqrt((output_data.averages(j).stdev - N*output_data.averages(j).data)/(N-1));
    end;
  end;
end;

% calculate noise covariance matrix
if parms.concat_ncov_flag
  if ~isempty(noise_matrix)
    output_data.noise.covar = cov(noise_matrix');
  end;
else
  output_data.noise.covar = output_data.noise.covar/output_data.noise.num_trials;
end;

fprintf('%s: finished\n',mfilename);
