function [data_out,badbaselines] = ts_freqanalysis_makebaseline (data_in,varargin)
% function [data_out,badbaselines] = ts_freqanalysis_makebaseline (data_in,varargin)
% 
% Usage: [data_out,badbaselines] = ts_freqanalysis_makebaseline (data_in,'option',value...
%
% This function performs frequency analysis on the given epoch_data set and
% then returns a baseline_data data set across all conditions.
%
% Required Input:
%
% epoch_data - a valid Time Surfer epoch data set.
%
% Optional Input:
%
% reject_flag    - perform baseline artifact rejection (default = 1)
% reject_exclude - remove bad trials in one of the following ways:
%                  'inall'           - remove any bad trial across all
%                                      channels
%                  'groupchannels'   - remove bad trials by grouping
%                                      channels (by name) and any trial
%                                      bad in the group is
%                                      excluded across that group
%                  'bychannel'(default)- remove bad trials only in the
%                                      channel that it was bad in
% bl_threshold   - threshold for rejection. theshold x 1st decile of power values (default = 3)
% baseline_start - start time point of the baseline in ms (default = -80 ms)
% baseline_end   - end time point of the baseline in ms (default = -5 ms)
%
%    Preprocessing Options (done prior to time-frequency analysis)
%
%      detrend  - 'yes' or 'no' {default = 'yes'}
%
%    Options for each method
%
%    WLTCONVOL - performs time-frequency analysis on any time series trial data
%                using the 'wavelet method' based on Morlet wavelets.
% 
%      output   - power ('pow') or power and cross-spectra ('powandcsd')
%                 {default = pow}
%      width    - width of the wavelet determines the temporal and spectral
%                 resolution of the analysis {defualt = half of foi}
%      gwidth   - length of the used wavelets in standard deviations
%                 of the implicit Gaussian kernel and should be choosen
%                 >= 3; {default = 3}
%
%    TFR -  computes time-frequency representations of single-trial
%           data using a convolution in the time-domain with Morlet's wavelets.
%
%      width - 'width' of wavelets expressed in cycles (default = 7)
%      downsample - ratio for downsampling which occurs after convolution
%
%    MTMFFT - performs frequency analysis on any time series
%             trial data using the 'multitaper method' (MTM) based on discrete
%             prolate spheroidal sequences (Slepian sequences) as tapers. Alternatively,
%             you can use conventional tapers (e.g. Hanning).
%
%      output    - power ('pow'), power and cross-spectra ('powandcsd'), complex 
%                  Fourier spectra ('fourier') {default = 'power'}
%      taper     - dpss, hanning, ect.. (see WINDOW) {default = 'dpss'}
%      tapsmofrq - number, the amount of spectral smoothing through
%                  multi-tapering. Note that 4 Hz smoothing means
%                  plus-minus 4 Hz, i.e. a 8 Hz smoothing box. {default = 4}
%      keeptapers- 'yes' or 'no' - return individual tapers or average
%      pad       - length in ms to which the data can be padded out
%                  {default = 'maxperlen'}
%
%    MTMCONVOL - performs time-frequency analysis on any time series trial data
%                using the 'multitaper method' (MTM) based on Slepian sequences as tapers. 
%                Alternatively, you can use conventional tapers (e.g. Hanning).
%       
%      In addition to the options for MTMFFT:
%  
%      t_ftimwin - length of time window in ms (vector the length of number
%                  of foi's)
%
%    MTMWELCH - performs frequency analysis on any time series
%               trial data using the 'multitaper method' (MTM) based on discrete
%               prolate spheroidal sequences (Slepian sequences) as tapers. Alternatively,
%               you can use conventional tapers (e.g. Hanning).
%  
%               Besides multitapering, this function uses Welch's averaged, modified
%               periodogram method. The data is divided into a number of sections with
%               overlap, each section is windowed with the specified taper(s) and the
%               powerspectra are computed and averaged over the sections in each trial.
%
%      In addition to the options for MTMCONVOL:
%
%      toi       - the times on which the analysis windows should be
%                  centered (in ms)
% Output:
%
% baseline_data  - a baseline data set that is averaged across all (good)
%                  trials of all conditions
%
%   NOTE: the num_trials field may be an array dependent on the the
%   reject_exclude parameter
%
% Optional Output:
%
% badbaselines   - a cell array (the size of the number of channels) that
%                  has a list of bad trials found for each channel

%% Check Inputs 

if nargin < 1, help(mfilename); end

error(nargoutchk(1,2,nargout,'string'));

parms = mmil_args2parms(varargin,...
                        {...
                         'method','wltconvol',[],...
                         'baseline_events',[],[],...
                         'baseline_start',-80,[],...
                         'baseline_end',-5,[],...
                         'foi',[1:1:60],[],...
                         'channels',      [],   [],...
                         'chantype',      [], {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                         'channelcmb',    {'all' 'all'},[],...
                         'logfile',[],[],...
                         'logfid',1,[],...
												 'trials_flag', 0, sort([false true])...
                       },...
                       false);

%% Rejection parms
rejparms = mmil_args2parms(varargin,...
                        {...
                         'reject_flag',true,sort([false true]),...
                         'reject_exclude','bychannel',{'inall','groupchannels','bychannel'},...
                         'bl_threshold',10,[],...
                         'chan_threshold',80,[],...
                         },...
                         false);

% temporarily force these flag values.  Jason Sherfey / 31-Oct-2008
parms.trials_flag = 0;
rejparms.reject_flag = 1;

rejargs  = mmil_parms2args(rejparms);                    

  
errors = ts_checkdata(data_in);
if ~isempty(errors)
 mmil_error(parms, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:}));
end;

if ~isfield(data_in,'epochs')
    mmil_error(parms,'Input data must be epoch_data.');
end


if (isempty(parms.baseline_events))
   parms.conditions  = num2cell(1:length(data_in.epochs));
else
    if ~iscell(parms.baseline_events), parms.baseline_events = num2cell(parms.baseline_events); end;
    if (~isempty(setdiff([parms.baseline_events{:}], [data_in.epochs.event_code])))
     mmil_error(parms, 'Event code doesn''t exist in epoch data: %d.', ...
       min(setdiff([parms.baseline_events{:}], [data_in.epochs.event_code])));
   else
     [a,parms.conditions]=intersect([data_in.epochs.event_code], [parms.baseline_events{:}]);
   end
end

if (~iscell(parms.baseline_events)), parms.baseline_events = num2cell(parms.baseline_events); end;
if (~iscell(parms.conditions)),  parms.conditions  = num2cell(parms.conditions);  end;
 
% Check channel selections
 if (~isempty(parms.chantype) && ~isempty(parms.channels))
   mmil_error(parms, 'Cannot specify BOTH chantype AND channels param values.');

   % Choose all channels by default
 elseif (isempty(parms.chantype) && isempty(parms.channels))
   parms.channels = 1:data_in.num_sensors;

   % Convert all inputs directly to channel indices
 elseif (isempty(parms.channels))
   switch parms.chantype
     case {'mag' 'grad1' 'grad2' 'eeg', 'other'}
       parms.channels = find(strcmp(parms.chantype,{data.sensor_info.typestring}));
     case {'grad'}
       parms.channels = find(strncmp(parms.chantype,{data.sensor_info.typestring},...
         length(parms.chantype)));
     case 'meg'
       [a,parms.channels] = find(ismember({data.sensor_info.typestring}, ...
         {'mag', 'grad1', 'grad2'}));
     case 'all'
       parms.channels = setdiff(1:data.num_sensors, find(strcmp('other', {data.sensor_info.typestring})));
   end;
 end;

parms.toi = [parms.baseline_start:(1/data_in.sfreq)*1000:parms.baseline_end];
                                        % ^---sec -> ms

%% Frequency Analysis Options 

switch (parms.method)
  case 'wltconvol'
     faparms = mmil_args2parms (varargin, ...
                            { 'output' , 'pow', {'pow','powandcsd'},...
                              'width',[],[],...% parms.foi/2, [],...
                              'gwidth', [],[]...%3, []...
                            },...
                           false);
     faparms.method = parms.method;
     faparms.toi = parms.toi / 1000;                % convert ms to sec for FieldTrip
     faparms.foi = parms.foi;
     faparms.channelcmb = parms.channelcmb;
  case 'tfr'
    faparms = mmil_args2parms (varargin, ...
          { 'width', 7, [],...
            'downsample',1, []...
          },...
      false);
    faparms.waveletwidth = width;
    faparms = rmfield(faparms,'width');
    faparms.method = parms.method;
    faparms.foi = parms.foi;
  case 'mtmfft'
    faparms = mmil_args2parms (varargin, ...
          { 'output', 'pow', {'pow','powandcsd','fourier'},...
            'taper', 'dpss', {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin',...
                              'chebwin','flattopwin','gausswin','hamming','hann','kaiser','nuttallwin',...
                              'parzenwin','rectwin','tukeywin','triang','dpss'},...
            'tapsmofrq', 4, [],...
            'keeptapers','no',{'yes','no'},...
            'pad','maxperlen',[] ...
          },...
      false);
    faparms.method     = parms.method;
    faparms.channelcmb = parms.channelcmb;
    faparms.foilim     = [parms.foi(1) parms.foi(end)];    
  case 'mtmconvol'
    faparms = mmil_args2parms (varargin, ...
          { 'output', 'pow', {'pow','powandcsd','fourier'},...
            'taper', 'dpss', {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin',...
                              'chebwin','flattopwin','gausswin','hamming','hann','kaiser','nuttallwin',...
                              'parzenwin','rectwin','tukeywin','triang','dpss'},...
            't_ftiwin',[],[],...
            'tapsmofrq', 4, [],...
            'keeptapers','no',{'yes','no'},...
            'pad','maxperlen',[] ...
          },...
      false);
    faparms.method     = parms.method;
    faparms.channelcmb = parms.channelcmb;
    faparms.foi        = parms.foi;    
    faparms.toi        = parms.toi / 1000;
    faparms.t_ftiwin   = faparms.t_ftiwin / 1000;
  case 'mtmwelch'
    faparms = mmil_args2parms (varargin, ...
          { 'output', 'pow', {'pow','powandcsd','fourier'},...
            'taper', 'dpss', {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin',...
                              'chebwin','flattopwin','gausswin','hamming','hann','kaiser','nuttallwin',...
                              'parzenwin','rectwin','tukeywin','triang','dpss'},...
            't_ftiwin',[],[],...
            'tapsmofrq', 4, [],...
            'keeptapers','no',{'yes','no'},...
            'pad','maxperlen',[] ...
          },...
      false);
    faparms.method     = parms.method;
    faparms.channelcmb = parms.channelcmb;
    faparms.foi        = parms.foi;    
    faparms.toi        = parms.toi / 1000;
    faparms.t_ftiwin   = faparms.t_ftiwin / 1000;
  otherwise
    mmil_error(parms, 'Method not recognized: %s', parms.tool);
end
faparms.trials_flag = parms.trials_flag;
faparms.channel = 'all';

if data_in.sfreq < 2*(max(parms.foi))
    mmil_error(parms,'The sampling frequency (%d) of the data is not high enough to compute frequencies as high as %d Hz.',...
        data_in.sfreq,max(parm.foi));
end  

%for i = 1:length(data_in.epochs)
%    mintimeneeded = faparms.toi(1) - ((1/faparms.foi(1))/2);
%    maxtimeneeded = ((1/faparms.foi(1))/2) + faparms.toi(end);
%    if mintimeneeded < data_in.epochs(i).time(1) || ...
%       maxtimeneeded > data_in.epochs(i).time(end)        
%       
%        mmil_logstr(parms,'WARNING: There is not enough padding of the data to calculate such low frequencies.');
%        mmil_logstr(parms,'For %d Hz and toi of %0.5f to %0.5f seconds, the data should range from %0.5f to %0.5f seconds.',...
%              faparms.foi(1),faparms.toi(1),faparms.toi(end),mintimeneeded,maxtimeneeded);
%        minallowabletime = max(data_in.epochs(i).time(1) + ((1/faparms.foi(1))/2),faparms.toi(1));
%        maxallowabletime = min(data_in.epochs(i).time(end) - ((1/faparms.foi(1))/2),faparms.toi(end));
%        faparms.toi = minallowabletime:1/data_in.sfreq:maxallowabletime;         
%        mmil_logstr(parms,'Set time of interest to: %0.5f to %0.5f seconds.',minallowabletime,maxallowabletime);
%        
%    end
%end

%% Parse out preprocessing parameters

ppparms = mmil_args2parms (varargin, ...
          { 
           'detrend','yes',{'yes','no'},...
           'lnfilter','yes',{'yes','no'},...
           'lnfreq',60,[],...
					 'blc','no',{'yes','no'},...
					 'blcwindow',[],[],...
          },...
          false);
ppparms.keeptrials = 'yes';
ppparms.feedback   = 'no' ;

    
%% Make the epochs into one condition

mmil_logstr(parms,'Combining trials across conditions %s.',num2str(cell2mat(parms.conditions)));
old_epoch_data    = data_in.epochs;
data_in.epochs = [];
for e = 1:length(parms.conditions)
    if e == 1
        data_in.epochs            = old_epoch_data(parms.conditions{e});
        data_in.epochs.event_code = 1;
    else % for rest combine across
        data_in.epochs.data               = cat(3,data_in.epochs.data,old_epoch_data(parms.conditions{e}).data);
        data_in.epochs.num_trials         = data_in.epochs.num_trials         + old_epoch_data(parms.conditions{e}).num_trials;
        data_in.epochs.num_rejects.mag    = data_in.epochs.num_rejects.mag    + old_epoch_data(parms.conditions{e}).num_rejects.mag;
        data_in.epochs.num_rejects.grad   = data_in.epochs.num_rejects.grad   + old_epoch_data(parms.conditions{e}).num_rejects.grad;
        data_in.epochs.num_rejects.eeg    = data_in.epochs.num_rejects.eeg    + old_epoch_data(parms.conditions{e}).num_rejects.eeg;
        data_in.epochs.num_rejects.eog    = data_in.epochs.num_rejects.eog    + old_epoch_data(parms.conditions{e}).num_rejects.eog;
        data_in.epochs.num_rejects.manual = data_in.epochs.num_rejects.manual + old_epoch_data(parms.conditions{e}).num_rejects.manual;
        data_in.epochs.num_rejects.skip   = data_in.epochs.num_rejects.skip   + old_epoch_data(parms.conditions{e}).num_rejects.skip;
    end
end 
clear old_epoch_data

%% Frequency Analysis

mmil_logstr(parms,'Converting to fieldtrip.');
mmil_logstr(parms, 'Running t/f analysis, method=%s, on %d frequencies and %d times.',...
     faparms.method, length(faparms.foi), length(faparms.toi));
ft_epochs = ts_data2fieldtrip(data_in,...
    'condition',1,...
    'channels',parms.channels,...
    'dimord','chan_time');
warning off MATLAB:rankDeficientMatrix
ft_timelock = timelockanalysis(ppparms,ft_epochs);
clear ft_epochs;
% Keep track of all the channel indexes for the future
[dum,in_ts,in_ft] = intersect({data_in.sensor_info.label},ft_timelock.label);
[dum,ex_ts]       = setdiff({data_in.sensor_info.label},ft_timelock.label);
[in_ft,idx]       = sort(in_ft);
in_ts             = in_ts(idx);
if rejparms.reject_flag
    %% Process by channel to find the bad trials
    mmil_logstr(parms,'Processing individual channels to find artifacts.');
    faparms.keeptrials = 'yes';
    %% Modify the data_in set for channel by channel analysis
    epoch_template             = rmfield(data_in,'epochs');
    epoch_template.epochs      = rmfield(data_in.epochs,'data');
    epoch_template.epochs.data = [];
    epoch_template.sensor_info = [];
    badtrials    = [];
    tic
    for ch = 1:length(ft_timelock.label)
        mmil_logstr(parms,'Channel %s (%d/%d).',ft_timelock.label{ch},ch,length(ft_timelock.label));
        faparms.channel            = ch;
        warning off MATLAB:NonIntegerInput
        ft_timefreq                = freqanalysis(faparms,ft_timelock,1);
        epoch_template.sensor_info = data_in.sensor_info(in_ts(ch));; 
        data_out_temp(ch)          = ts_fieldtrip2timefreq (ft_timefreq,epoch_template);
        clear ft_timefreq
        %% Artifact Rejection %%
        [data_out_temp(ch),badtrials_bychan{ch}] = ts_freqanalysis_bl_reject (data_out_temp(ch),rejargs{:});
        if (length(badtrials_bychan{ch}) / data_in.epochs.num_trials) >= (rejparms.chan_threshold/100)
                data_in.sensor_info(in_ts(ch)).badchan = 1;
                mmil_logstr(parms,'Marking channel %s as having a bad baseline period.',data_out_temp(ch).sensor_info.label);
        end;
        if strcmpi(rejparms.reject_exclude,'inall')
            %% If it is bad channel don't let its trials ruin the whole set
            if data_in.sensor_info(in_ts(ch)).badchan ~= 1
               badtrials = unique([badtrials badtrials_bychan{ch}]);            
            end
            %% Break early if number of remaining trials reaches 0
            if (length(badtrials) >= data_in.epochs.num_trials)
                mmil_error(parms,'There are no good baselines, consider increasing the ''threshold'' parameter.');
            end
            mmil_logstr(parms,'Total bad trials: %.0f.',length(badtrials));
        end
        epoch_template.sensor_info = [];
    end
    mmil_logstr(parms,'Elapsed time %.0f seconds.',toc);   
    switch rejparms.reject_exclude
        case 'inall'
            clear data_out_temp
            data_out = remove_across_all(data_in,1,badtrials,ppparms,faparms);
        case 'groupchannels'
            data_out = remove_by_group  (data_in,1,badtrials_bychan,data_out_temp,ppparms,faparms,in_ts,in_ft);
        case 'bychannel'
            data_out = remove_by_channel(data_in,1,data_out_temp,in_ts,in_ft,badtrials_bychan);
        otherwise
    end
    tmp = badtrials_bychan;
    badtrials_bychan = [];
    for c =1:length(in_ts)
        badtrials_bychan{in_ts(c)} = tmp{in_ft(c)};
    end
    badbaselines = badtrials_bychan;
else
    faparms.keeptrials  = 'no';
    warning off MATLAB:NonIntegerInput
    ft_timefreq         = freqanalysis(faparms,ft_timelock,parms.trials_flag);
    clear ft_timelock
    data_out            = ts_fieldtrip2timefreq(ft_timefreq,data_in);
    clear ft_timefreq
    for c = 1:length(data_in.sensor_info)
        badbaselines{c} = [];
    end
end

errors = ts_checkdata(data_out);
if ~isempty(errors)
  mmil_error(parms, 'Errors in baseline data structure: %s.', sprintf('\t%s\n', errors{:}));
end;

warning on all

%% SUB FUNCTIONS

function data_out = remove_across_all (data_in,condition,badtrials,ppparms,faparms)
% reruns the freqency analysis on the data set removing the badtrials
% across all the channels
mmil_logstr(parms,'Removing %s of %s trials across all channels.',...
    num2str(length(badtrials)),num2str(data_in.epochs.num_trials));
% Remove the trials from the orginal data
tmp_data            = data_in.epochs.data;
data_in.epochs.data = [];
data_in.epochs.data = tmp_data(:,:,...
    setdiff([1:size(tmp_data,3)],badtrials));
data_in.epochs.num_rejects.eeg = data_in.epochs.num_rejects.eeg + length(badtrials);
data_in.epochs.num_trials      = data_in.epochs.num_trials      - length(badtrials);
mmil_logstr(parms,'Reprocessing across whole data set without bad trials.');
ft_epochs = ts_data2fieldtrip(data_in,...
    'condition',1,...
    'channels',parms.channels,...
    'dimord','chan_time');
warning off MATLAB:rankDeficientMatrix
ft_timelock = timelockanalysis(ppparms,ft_epochs);
clear ft_epochs;
faparms.keeptrials = 'no';
faparms.channel    = 'all';
warning off MATLAB:NonIntegerInput
ft_timefreq                = freqanalysis(faparms,ft_timelock,faparms.trials_flag);
data_out                   = ts_fieldtrip2timefreq(ft_timefreq,data_in);
clear ft_timefreq

function data_out = remove_by_group (data_in,condition,badtrials_bychan,data_out_temp,ppparms,faparms,in_ts,in_ft)
% remove trials across groups of channels
data_out                     = data_out_temp(1);
data_out.sensor_info         = data_in.sensor_info;
data_out.timefreq.num_trials = repmat(data_out.timefreq.num_trials,[1 length(data_out.sensor_info)]);
data_out.timefreq.data       = nan(length(data_in.sensor_info),...
                                   size(data_out_tmp(1).timefreq.power,2),...
                                   size(data_out_tmp(1).timefreq.power,3));
data_out.timefreq.power      = nan(length(data_in.sensor_info),...
                                   size(data_out_tmp(1).timefreq.data,2),...
                                   size(data_out_tmp(1).timefreq.power,3));
grouped                      = group_channels(data_in.sensor_info);
for g = 1:length(grouped)
   [dum,group_indxs]    = intersect(in_ts,grouped(g).indices);
   badtrialsbygroup     = unique([badtrials_bychan{in_ft(group_indxs)}]); 
   tmp_data             = data_in;
   tmp_data.epochs.data = data_in.epochs.data(:,:,setdiff([1;size(data_in.epochs.data,3)],badtrialsbygroup));
   mmil_logstr(parms,'Reprocessing %s channel(s) without %d bad trials across the group.',grouped(g).name,length(badtrialsbygroup));
   %% Only process the selected channels in the group
   ft_epochs   = ts_data2fieldtrip (tmp_data,...
                                     'condition',1,...
                                     'channels',grouped(g).indices,...
                                     'dimord','chan_time');
   warning off MATLAB:rankDeficientMatrix
   ft_timelock = timelockanalysis(ppparms,ft_epochs);
   clear ft_epochs;
   faparms.keeptrials = 'no';
   faparms.channel = 'all';
   warning off MATLAB:NonIntegerInput
   ft_timefreq  = freqanalysis(faparms,ft_timelock,faparms.trials_flag);
   clear ft_timelock
   tmp_timefreq = ts_fieldtrip2timefreq(ft_timefreq,data_in);
   clear ft_timefreq
   [a,includedchans] =intersect({data_out.sensor_info.label}, ft_timefreq.label);
   includedchans     = sort(includedchans);
   %% Sort it into the proper locations
   data_out.timefreq.power(includedchans,:,:)  = tmp_timefreq.timefreq.power(includedchans,:,:);
   data_out.timefreq.data (includedchans,:,:)  = tmp_timefreq.timefreq.data(includedchans,:,:);
   data_out.timefreq.num_trials(includedchans) = data_in.epochs.num_trials - length(badtrialsbygroup);
   clear tmp_data tmp_timefreq
end
[data_out.timefreq.num_trials(find([data_out.sensor_info.badchan]))] = deal(0);

function data_out = remove_by_channel(data_in,condition,data_out_tmp,in_ts,in_ft,badtrials)
% remove trials by channel - just uses the data_out_tmp
data_out                     = data_out_tmp(1);
data_out.sensor_info         = data_in.sensor_info;
data_out.timefreq.num_trials = repmat(data_out.timefreq.num_trials,[1 length(data_out.sensor_info)]);
data_out.timefreq.data       = nan(length(data_in.sensor_info),...
                                   size(data_out_tmp(1).timefreq.power,2),...
                                   size(data_out_tmp(1).timefreq.power,3));
data_out.timefreq.power      = nan(length(data_in.sensor_info),...
                                   size(data_out_tmp(1).timefreq.data,2),...
                                   size(data_out_tmp(1).timefreq.power,3));
for ch = 1:length(in_ft)
    data_out.timefreq.data (in_ts(ch),:,:)  = data_out_tmp(in_ft(ch)).timefreq.data(1,:,:);
    data_out.timefreq.power(in_ts(ch),:,:)  = data_out_tmp(in_ft(ch)).timefreq.power(1,:,:);
    data_out.timefreq.num_trials(in_ts(ch)) = data_out.timefreq.num_trials(in_ts(ch)) - length(badtrials{in_ft(ch)});
end
[data_out.timefreq.num_trials(find([data_out.sensor_info.badchan]))] = deal(0);

function grouped = group_channels(sensor_info)
% Group the channels
groups = unique(regexp([sensor_info.label],'[A-Za-z]*(?=[\d*])','match'));
for g = 1:length(groups)
    exp = sprintf('%s(?=\\d)',groups{g});
    grouped(g).name     = groups{g};
    grouped(g).indices = [];
    for ch = 1:length(sensor_info)
        [indx1,indx2]=regexp(sensor_info(ch).label,exp,'start','end');
        if isempty(indx1), indx1 = 0; end
        if isempty(indx2), indx2 = 0; end
        if (indx1 == 1) && (indx2 == length(groups{g}))
            grouped(g).indices(end+1) = ch;
        end
    end
end
if isempty(grouped)
    %% Because there are no groups this will essentially be
    %% running by channel inform user but run anyway
    mmil_logstr(parms,'There were no groups of channels found for rejection. Run with ''exclude'' set to ''bychannel'' to save time.');
    grouped = [];
    for g = 1:length(sensor_info)
        grouped(g).name    = sensor_info(g).label;
        grouped(g).indices = g;
    end
else
    leftovers = setdiff([sensor_info.lognum],[grouped.indices]);
    for l = length(leftovers)
        grouped(end+1).name  = sensor_info(leftovers(l)).label;
        grouped(end).indices = leftovers(l);
    end
end
mmil_logstr(parms,'Found %d groups of channels to reject by.',length(grouped));








