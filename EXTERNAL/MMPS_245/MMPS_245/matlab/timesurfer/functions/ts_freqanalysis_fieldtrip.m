function [data_out,badtrials,baseline_data,badbaselines] = ts_freqanalysis_fieldtrip (data_in,varargin)
% function [data_out,badtrials,baseline_data,badbaselines] = ts_freqanalysis_fieldtrip (data_in,varargin)
% 
% Usage: [data_out,badtrials] = ts_freqanalysis_fieldtrip (data_in,'option',value...
%
% This function performs frequency analysis on the given epoch_data set and
% then returns timefreq_data structure
%
% Required Input:
%
% epoch_data - a valid Time Surfer epoch data set.
%
% Optional Input:
%
% reject_flag    - perform baseline artifact rejection (default = 1)
% reject_exclude - remove bad trials in one of the following ways:
%                  'inall' (default) - remove any bad trial across all
%                                      channels
%                  'groupchannels'   - remove bad trials by grouping
%                                      channels (by name) and any trial
%                                      bad in the group is
%                                      excluded across that group
%                  'bychannel'       - remove bad trials only in the
%                                      channel that it was bad in
% trial_threshold- threshold for rejection. z-score compared to baseline
% chan_threshold - % of trials bad causes channel to be marked as bad
% baseline_data  - baseline data set created by ts_freqanalysis_makebaseline
%    OR
% baseline_threshold - calculate baseline data during the run and use this
%                      threshold for rejecting baesline trials - see
%                      ts_freqanalysis_makebaseline for more information
% baseline_start  - starting time point for baseline (in ms)
% baseline_end    - ending time point for baseline (in ms)
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
%      sf
%      st
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
% timefreq_data  - a baseline data set that is averaged across all (good)
%                  trials of all conditions
%
%   NOTE: the num_trials field may be an array dependent on the the
%   reject_exclude parameter
%
% Optional Output:
%
% badtrials      - a cell array (the size of the number of channels) that
%                  has a list of bad trials found for each channel
% baseline_data  - if created 
% badabaselines  - if created
%
% Created: 					 by Rajan Patel
% Rcnt Mod:	09/13/08 by Jason Sherfey 
% Last Mod: 09/15/12 by Don Hagler
%

%
% Revision 05-Sep-2008 by Jason Sherfey
%		Added parms.trials_flag and save_trials() to save individual trial data per channel with rejection.
%
% Revision 13-Sep-2008 by Jason Sherfey
%  	Added option to save trials without rejection
%
% Revision 19-Sep-2008 by Jason Sherfey
% 	Modified to work with Andrei's modified versions of freqanalysis & freqanalysis_wltconvol_cmplx
%		12:03 - modified for the case: (reject_flag=0, trials_flag=0)
%  	12:29 - commented out the resetting toi if not enough padding

%% Check Inputs 

if nargin < 1, help(mfilename); end

error(nargoutchk(1,4,nargout,'string'));

parms = mmil_args2parms(varargin,...
                        {...
                         'savesingles_flag',0,[0 1],...
                         'savecomplex_flag',1,[0 1],...
                         'method','wltconvol',[],...
                         'foi','all',[],...
                         'sf',[],[],...
                         'toi',[],[],...
                         'toilim',[],[],...
                         'event_codes',[],[],...                        
                         'events',[],[],...
                         'conditions',[],[],...         
                         'channels',      [],   [],...
                         'chantype',      [], {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                         'channelcmb',    {'all' 'all'},[],...
                         'verbose',1,{0,1},...
                         'logfile',[],[],...
                         'logfid',1,[],...
                         'freq_path',[],[],...
                         'datafile',[],[],...
												 'datapath',[],[],...
                         'filename',[],[],...
												 'f_name','',[],...
												 'trials_flag', 0, sort([false true]),...
                         'save_flag', 1, sort([false true]),...
												 'makebaseline_flag', 0, sort([false true]),...
												 'visualreject_flag', 0, sort([false true])...
                         'overwrite',0,{0,1},...
                       },...
                       false);

if isempty(parms.event_codes) && ~isempty(parms.events)
  parms.event_codes = parms.events;
end
data_in  = ts_checkdata_header(data_in,'precision','double','events',parms.event_codes);
% the next line needs to be tested for bugs
data_in  = ts_data_selection(data_in,'channels',parms.channels,'chantype',parms.chantype); 

if ~isempty(parms.filename)
  [pathstr,name,ext] = fileparts(parms.filename{1});
  parms.freq_path = pathstr;
  parms.f_name = parms.filename{1};
elseif isempty(parms.f_name)
  parms.f_name = parms.datafile; 
end
if ~iscell(parms.f_name), parms.f_name = {parms.f_name}; end
% if isempty(parms.toi), parms.toi = data_in.epochs(1).time*1000; end

errors = ts_checkdata(data_in);
if ~isempty(errors)
 mmil_error(parms, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:}));
end;

if ~isfield(data_in,'epochs')
    mmil_error(parms,'Input data must be epoch_data.');
end

% Set foi if necessary
foi = [2:12 14:2:24 25:5:55 70:10:200];
sf  = [1 1 2*ones(1,9) 3*ones(1,6) 5*ones(1,7) 10*ones(1,14)];
if ischar(parms.foi)
  if strcmpi(parms.foi,'megshort')
    parms.foi = foi(3:24);
    parms.sf  = sf (3:24);
  elseif strcmpi(parms.foi,'meg') || strcmpi(parms.foi,'meglong') ...
      || (~isempty(parms.chantype) && ismember(parms.chantype,{'mag','grad','grad1','grad2','meg'}))
    parms.foi = foi(1:24);
    parms.sf  = sf (1:24);
  elseif strcmpi(parms.foi,'ieegshort')
    parms.foi = foi(3:end);
    parms.sf  = sf (3:end);
  elseif strcmpi(parms.foi,'all') || strcmpi(parms.foi,'ieeg') || strcmpi(parms.foi,'ieeglong') ...
      || (~isempty(parms.chantype) && ismember(parms.chantype,{'eeg','ieeg'})) ...
      || (~isempty(parms.chantype) && ismember('eeg',{data_in.sensor_info.typestring})) ...
      || (~isempty(parms.chantype) && ismember('ieeg',{data_in.sensor_info.typestring}))
    parms.foi = foi;
    parms.sf  = sf;
  else
    parms.foi = foi;
    parms.sf  = sf;
  end
elseif isempty(parms.sf) || ~all(ismember(parms.sf,sf)) || ~all(ismember(parms.foi,foi)) || length(parms.sf)~=length(parms.foi)    
  mmil_logstr(parms,'Not using the TimeSurfer-recommended frequencies of interest and spectral resolution!');
%   warning('Not using the TimeSurfer-recommended frequencies of interest and spectral resolution!');
end
if isempty(parms.foi) || ischar(parms.foi)
  parms.foi = 1:60;
end

% Check channel selections
 if (~isempty(parms.chantype) && ~isempty(parms.channels))
%    mmil_error(parms, 'Cannot specify BOTH chantype AND channels param values.');
    parms.channels = [];
 end
   % Choose all channels by default
 if (isempty(parms.chantype) && isempty(parms.channels))
   parms.channels = 1:data_in.num_sensors;

   % Convert all inputs directly to channel indexes
 elseif (isempty(parms.channels))
   switch parms.chantype
     case {'mag' 'grad1' 'grad2' 'eeg', 'other'}
       parms.channels = find(strcmp(parms.chantype,{data_in.sensor_info.typestring}));
     case {'grad'}
       parms.channels = find(strncmp(parms.chantype,{data_in.sensor_info.typestring},...
         length(parms.chantype)));
     case 'meg'
       [a,parms.channels] = find(ismember({data_in.sensor_info.typestring}, ...
         {'mag', 'grad1', 'grad2'}));
     case 'all'
       parms.channels = setdiff(1:data_in.num_sensors, find(strcmp('other', {data_in.sensor_info.typestring})));
   end;
 end;

 if ~isempty(parms.event_codes) && isnumeric(parms.event_codes) && ~iscell(parms.event_codes)
   parms.event_codes = {parms.event_codes};
 end
 % Make event codes a cell array
 if (~isempty(parms.event_codes) && ~isempty(parms.conditions))
   mmil_error(parms, 'Cannot specify BOTH event_codes AND conditions.');

   %specified none
 elseif (isempty(parms.event_codes) && isempty(parms.conditions))
   parms.event_codes = { data_in.epochs.event_code };
   parms.conditions  = num2cell(1:length(data_in.epochs));

   % specified conditions
 elseif (isempty(parms.event_codes))
   if (min(parms.conditions) < 1 || max(parms.conditions) > length(data_in.epochs))
     mmil_error(parms, 'Conditions are out of range; nConditions=%d', length(data_in.epochs));
   else
     parms.event_codes = { data_in.epochs(parms.conditions).event_code };
   end;

   % specified event_codes
 else
   if (~isempty(setdiff([parms.event_codes{:}], [data_in.epochs.event_code])))
     mmil_error(parms, 'Event code doesn''t exist in epoch data: %d.', ...
       min(setdiff([parms.event_codes{:}], [data_in.epochs.event_code])));
   else
     [a,parms.conditions]=intersect([data_in.epochs.event_code], [parms.event_codes{:}]);
     parms.conditions    = num2cell(parms.conditions);
   end;
 end;

 % Make sure both conditions and event_codes are cell arrays
 if (~iscell(parms.event_codes)), parms.event_codes = num2cell(parms.event_codes); end;
 if (~iscell(parms.conditions)),  parms.conditions  = num2cell(parms.conditions);  end;
 
 if ~isempty(parms.toilim)
   parms.toi = parms.toilim;
 end
 if isempty(parms.toi)
     for i = 1:length(parms.conditions)
        if i==1
          first_time  = data_in.epochs(parms.conditions{i}).time(1);
          last_time   = data_in.epochs(parms.conditions{i}).time(end); 
        else
         first_time = min(first_time,data_in.epochs(parms.conditions{i}).time(1));
         last_time  = max(data_in.epochs(parms.conditions{i}).time(end));
        end
     end
     parms.toi = [first_time:1/data_in.sfreq:last_time]*1000;
 elseif length(parms.toi)==2
   parms.toi = [parms.toi(1):1/data_in.sfreq:parms.toi(2)]*1000;
 end

%% Frequency Analysis Options 

switch (parms.method)
  case 'wltconvol'
     faparms = mmil_args2parms (varargin, ...
                            { 'output' , 'pow', {'pow','powandcsd'},...
                              'width',[],[],...% parms.foi/2, [],...
                              'gwidth',[],[],...%3, []...
                              'st',[],[]...
                            },...
                           false);
     faparms.sf = parms.sf;
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

% if data_in.sfreq < 2*(max(parms.foi))
%     mmil_error(parms,'The sampling frequency (%d) of the data is not high enough to compute frequencies as high as %d Hz.',...
%         data_in.sfreq,max(parms.foi));
% end  

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
           'detrend','no',{'yes','no'},...
           'lnfilter','no',{'yes','no'},...
           'lpfilter','no',{'yes','no'},...
           'lpfreq',200,[],...
           'lnfreq',60,[],...
           'hpfilter','no',[],...
           'hpfreq',[],[],...
           'blc','no',[],...
           'blcwindow','all',[],...
          },...
          false);
if isnumeric(ppparms.blc)
  if ppparms.blc == 1
    ppparms.blc = 'yes'; 
  elseif ppparms.blc == 0
    ppparms.blc = 'no'; 
  end
end
ppparms.keeptrials = 'yes';
ppparms.feedback   = 'no' ;

ppparms.removemean = 'no';
%ppparms.channel    = {data_in.sensor_info(parms.channels).label};

%% Rejection parms

rejparms = mmil_args2parms(varargin,...
                        {...
                         'reject_flag',false,sort([false true]),...
                         'reject_exclude','bychannel',{'inall','groupchannels','bychannel'},...
                         'trial_threshold',50,[],...
                         'chan_threshold',80,[],...
                         'baseline_data',[],[],...
                         'bl_threshold',10,[],...
                         'baseline_start',-80,[],...
                         'baseline_end',-5,[],...
                         'baseline_events',[],[],...
                         },...
                         false);         

if 			rejparms.reject_flag		&&		parms.makebaseline_flag		&&		~isempty(rejparms.baseline_data)
	mmil_logstr(parms,'Creating new baseline data for use with rejection.');
	mmil_logstr(parms,'Note: Baseline data provided is not required since new baseline data will be created.');
elseif	~rejparms.reject_flag		&&		parms.makebaseline_flag		&&		~isempty(rejparms.baseline_data)
	mmil_logstr(parms,'Note: Creating new baseline data, but no rejection will be performed.');
	mmil_logstr(parms,'Note: Baseline data provided is not required since no rejection will be performed.');	
elseif	~rejparms.reject_flag		&&		~parms.makebaseline_flag		&&		~isempty(rejparms.baseline_data)
%	mmil_logstr(parms,'Note: Baseline data provided is not required since no rejection will be performed.');
elseif	rejparms.reject_flag		&&		~parms.makebaseline_flag		&&		~isempty(rejparms.baseline_data)
%	mmil_logstr(parms,'Baseline data provided will be used for rejection.');
elseif	rejparms.reject_flag		&&		parms.makebaseline_flag		&&		isempty(rejparms.baseline_data)
	mmil_logstr(parms,'Creating baseline data for use with rejection.');
elseif	~rejparms.reject_flag		&&		parms.makebaseline_flag		&&		isempty(rejparms.baseline_data)
	mmil_logstr(parms,'Note: Creating new baseline data, but no rejection will be performed.');
elseif	rejparms.reject_flag		&&		~parms.makebaseline_flag		&&		isempty(rejparms.baseline_data)
%	mmil_logstr(parms,'WARNING: Setting makebaseline to "yes" since rejection will be performed.');
	mmil_logstr(parms,'Creating baseline data for use with rejection.');	
	parms.makebaseline_flag = 1;
	if nargout < 2, 
		mmil_error(parms, 'Not enough output arguments.  See help ts_freqanalysis_fieldtrip for more information.'); 
	end
end

if parms.makebaseline_flag
    mmil_logstr(parms,'Creating baseline data set for use with rejection.');
		faparms.channel    = {data_in.sensor_info(parms.channels).label};
    faargs = mmil_parms2args(faparms); 
    [baseline_data,badbaselines] = ts_freqanalysis_makebaseline(data_in,...
                                                                faargs{:},...
                                                                'reject_flag',rejparms.reject_flag,...
                                                                'reject_exclude',rejparms.reject_exclude,...
                                                                'bl_threshold',rejparms.bl_threshold,...
                                                                'chan_threshold',rejparms.chan_threshold,...
                                                                'baseline_start',rejparms.baseline_start,...
                                                                'baseline_end',rejparms.baseline_end,...
                                                                'baseline_events',rejparms.baseline_events);
	
elseif rejparms.reject_flag
    baseline_data = rejparms.baseline_data;
else
    baseline_data = [];
    badbaselines  = [];
end

if isempty(rejparms.chan_threshold) || (rejparms.chan_threshold == 0)
    rejparms.chan_threshold = 100;
end

rejparms      = rmfield(rejparms,'baseline_data');
rejargs       = mmil_parms2args(rejparms);

orig_sensor_info = getfield(data_in,'sensor_info');

if ~isfield(parms,'filename')
  parms.filename = {};
end
%% Frequency Analysis
for i = 1:length(parms.conditions)
    mmil_logstr(parms,'Converting to fieldtrip.');
    mmil_logstr(parms, 'Running t/f analysis, method=%s, on %d frequencies and %d times.',...
          faparms.method, length(faparms.foi), length(faparms.toi));
		% guarantee no new channels will be set to badchan (so all channels will be analyzed each condition)
		tmp_sensor_info = data_in.sensor_info; 
		data_in.sensor_info = orig_sensor_info;
    if ~strcmp(class(data_in.epochs(parms.conditions{i}).data),'double')
      data_in.epochs(parms.conditions{i}).data = double(data_in.epochs(parms.conditions{i}).data);
    end
    ft_epochs = ts_data2fieldtrip(data_in,'condition',parms.conditions{i},'channels',parms.channels,'dimord','chan_time');
		data_in.sensor_info = tmp_sensor_info; 
		clear tmp_sensor_info;
    
    warning off MATLAB:rankDeficientMatrix
    ft_timelock = timelockanalysis(ppparms,ft_epochs);

		% Visual rejection
		if (parms.visualreject_flag)
    	cfg             = [];
	    cfg.channel     = num2cell(parms.channels);
	    cfg.keepchannel = 'yes';
	    cfg.method      = 'channel';
	    rej_data        = ft_timelock;
	    rej_ok          = 2;
	    while rej_ok == 2
	        rej_data = rejectvisual(cfg,rej_data);
	        rej_ok   = menu('Happy with artifact rejection?', 'Yes','No');
	    end;
	    ft_timelock = rej_data;
		end;

    %% Keep track of all the channel indexes for the future
    [dum,in_ts,in_ft] = intersect({data_in.sensor_info.label},ft_timelock.label);
    [dum,ex_ts]       = setdiff({data_in.sensor_info.label},ft_timelock.label);
    [in_ft,idx]       = sort(in_ft);
    in_ts             = in_ts(idx);

		clear ft_epochs;
    if rejparms.reject_flag
        %% Process by channel to find the bad trials
        mmil_logstr(parms,'Processing individual channels to find artifacts.');
        faparms.keeptrials = 'yes';        
        %% Modify the data_in set for channel by channel analysis
        epoch_template             = rmfield(data_in,'epochs');
        epoch_template.epochs      = rmfield(data_in.epochs(parms.conditions{i}),'data');
        epoch_template.epochs.data = [];
        epoch_template.sensor_info = [];
        %% setup baseline_data for by channel analysis
        bl_power                   = baseline_data.timefreq.power;
        bl_data                    = baseline_data.timefreq.data;
        if length(baseline_data.timefreq.num_trials) ~= 1
            bl_trials = baseline_data.timefreq.num_trials;
        else
            bl_trials = repmat(baseline_data.timefreq.num_trials,[1 length(baseline_data.sensor_info)]);
        end
        baseline_data.timefreq.num_trials = [];
        baseline_data.timefreq.power      = [];
        baseline_data.timefreq.data       = [];
        badtrials                         = [];
        tic
        for ch = 1:length(ft_timelock.label)				
            % mmil_logstr(parms,'Channel %s (%d/%d).',ft_timelock.label{ch},ch,length(ft_timelock.label));
            faparms.channel            = ch;
            warning off MATLAB:NonIntegerInput
            ft_timefreq                = freqanalysis(faparms,ft_timelock,1);
            epoch_template.sensor_info = data_in.sensor_info(in_ts(ch));
            data_out_temp(ch)          = ts_fieldtrip2timefreq (ft_timefreq,epoch_template);
						clear ft_timefreq; 						
						% if statement added by Jason Sherfey on 05-Sep-2008
						if parms.trials_flag && strcmpi(rejparms.reject_exclude,'bychannel'), tfchan = data_out_temp(ch); end

            %% Artifact Rejection %%
            % use only baseline data for this channel
            baseline_data.timefreq.power          = bl_power (in_ts(ch),:,:);
            baseline_data.timefreq.data           = bl_data  (in_ts(ch),:,:);
            baseline_data.timefreq.num_trials     = bl_trials(in_ts(ch));            
            %fprintf('Baseline channel: %s, Timefreq Channel: %s.\n',baseline_data.sensor_info(in_ts(ch)).label,ft_timelock.label{ch});
            [data_out_temp(ch),badtrials_bychan{ch}] = ts_freqanalysis_reject (data_out_temp(ch),baseline_data,rejargs{:});
            %% Reject a Channel if has over specified (chanthresh option) % of bad trials
            if (length(badtrials_bychan{ch}) / data_in.epochs(parms.conditions{i}).num_trials) >= (rejparms.chan_threshold/100)
                data_in.sensor_info(in_ts(ch)).badchan = 1;
                mmil_logstr(parms,'Marking channel %s as bad.',data_out_temp(ch).sensor_info.label);
            end;
            %badtrials = unique([badtrials badtrials_bychan{ch}]);
            %mmil_logstr(parms,'Total bad trials: %.0f.',length(badtrials));
            if strcmpi(rejparms.reject_exclude,'inall')
                %% If it is bad channel don't let its trials ruin the whole set
                if data_in.sensor_info(in_ts(ch)).badchan ~= 1
                    badtrials = unique([badtrials badtrials_bychan{ch}]);
                end
                badtrials = unique([badtrials badtrials_bychan{ch}]);
                %% Break early if number of remaining trials reaches 0
                if (length(badtrials) >= data_in.epochs(parms.conditions{i}).num_trials)
                    mmil_error(parms,'There are no good trials, consider increasing ''trial_threshold'' parameter.');
                end
                mmil_logstr(parms,'Total bad trials: %.0f.',length(badtrials));
            end
						% if statement added by Jason Sherfey on 05-Sep-2008
						if parms.trials_flag && strcmpi(rejparms.reject_exclude,'bychannel')
							tfchan = remove_by_channel(data_in,parms.conditions{i},tfchan,in_ts,in_ft,badtrials_bychan,faparms);
							parms.filename{end+1} = save_trials(tfchan,parms,faparms,rejparms);
							clear tfchan;
						end

            epoch_template.sensor_info = [];
        end	
        mmil_logstr(parms,'Elapsed time %.0f seconds.',toc);

        switch rejparms.reject_exclude
            case 'inall'
                clear data_out_temp
                data_out(i) = remove_across_all(data_in,parms.conditions{i},badtrials,ppparms,faparms);
                % reset channel information
                baseline_data.timefreq.num_trials = bl_trials(1);
            case 'groupchannels'
                data_out(i) = remove_by_group  (data_in,parms.conditions{i},badtrials_bychan,data_out_temp,ppparms,faparms,in_ts,in_ft);
                baseline_data.timefreq.num_trials = bl_trials;
            case 'bychannel'
                data_out(i) = remove_by_channel(data_in,parms.conditions{i},data_out_temp,in_ts,in_ft,badtrials_bychan,faparms);
                baseline_data.timefreq.num_trials = bl_trials;
            otherwise
        end
        tmp = badtrials_bychan;

        badtrials_bychan = [];
        for c =1:length(in_ts)
            badtrials_bychan{in_ts(c)} = tmp{in_ft(c)};
        end
        bt_out{i} = badtrials_bychan;
        % reset the baseline_data
        baseline_data.timefreq.power      = bl_power;
        baseline_data.timefreq.data       = bl_data;

				clear badtrials_bychan;		
    else
        tic         
        epoch_template             = rmfield(data_in,'epochs');
        epoch_template.epochs      = rmfield(data_in.epochs(parms.conditions{i}),'data');
        epoch_template.epochs.data = [];
        epoch_template.sensor_info = [];	        
        if parms.save_flag && parms.trials_flag && ~(isfield(faparms,'output') && strcmp(faparms.output,'powandcsd'))
					faparms.keeptrials = 'yes';
          %% Modify the data_in set for channel by channel analysis			
          for ch = 1:length(ft_timelock.label)
            faparms.channel = ch; 			
            warning off MATLAB:NonIntegerInput						
            ft_timefreq = freqanalysis(faparms,ft_timelock,1);				
            epoch_template.sensor_info = data_in.sensor_info(in_ts(ch));
            epoch_template.num_sensors = 1;
            tfchan      = ts_fieldtrip2timefreq (ft_timefreq,epoch_template);	
            clear ft_timefreq;
            parms.filename{end+1} = save_trials(tfchan,parms,faparms,rejparms);
            clear tfchan;
          end
				end; 
				if parms.trials_flag && (isfield(faparms,'output') && strcmp(faparms.output,'powandcsd'))
          chcmb = unique({faparms.channelcmb{:}});
          [jnk1 faparms.channel jnk2] = intersect(ft_timelock.label,chcmb);
          faparms.keeptrials = 'yes';
          faparms.trials_flag = 1;
          warning off MATLAB:NonIntegerInput
          ft_timefreq = freqanalysis(faparms,ft_timelock,1);
          clear ft_timelock
          [sel1 sel2] = match_str({data_in.sensor_info.label},chcmb);
          epoch_template.sensor_info = data_in.sensor_info(sel1);
          epoch_template.num_sensors = length(sel1);
          data_out(i) = ts_fieldtrip2timefreq(ft_timefreq,epoch_template);
        else
          if parms.trials_flag && ~parms.save_flag
            faparms.channel = 'all';
            faparms.keeptrials  ='yes'; 
            faparms.trials_flag = 1;
            warning off MATLAB:NonIntegerInput
            ft_timefreq         = freqanalysis(faparms,ft_timelock,1);		
            clear ft_timelock
            data_out(i) = ts_fieldtrip2timefreq(ft_timefreq,data_in);            
          else
            faparms.channel = 'all';
            faparms.keeptrials  ='no'; 
            faparms.trials_flag = 0;
            warning off MATLAB:NonIntegerInput
            ft_timefreq         = freqanalysis(faparms,ft_timelock,0);		
            clear ft_timelock
            data_out(i) = ts_fieldtrip2timefreq(ft_timefreq,data_in);
          end
        end
        data_out(i).timefreq.event_code = data_in.epochs(parms.conditions{i}).event_code;
        data_out(i).timefreq.num_rejects= data_in.epochs(parms.conditions{i}).num_rejects;
        data_out(i).timefreq.num_trials = data_in.epochs(parms.conditions{i}).num_trials;
        clear ft_timefreq
        for c = 1:length(data_in.sensor_info)
            bt_out{i}{c} = [];
        end
        mmil_logstr(parms,'Elapsed time %.0f seconds.',toc);
    end
end
% if ~parms.savecomplex_flag
%   data_out.timefreq = rmfield(data_out.timefreq,'cmplx');
%   data_out.timefreq = rmfield(data_out.timefreq,'data');
% end

%% Consolidate Conditions

if length(data_out) > 1
    data_out = ts_combine_data(data_out);
end
errors = ts_checkdata(data_out);
if ~isempty(errors)
  mmil_error(parms, 'Errors in final timefreq data structure: %s.', sprintf('\t%s\n', errors{:}));
end;
try parms.prepro   = ppparms; end
try parms.faparms  = faparms; end
try data_out.parms = parms;   end
badtrials = bt_out;

%keyboard;

warning on all


%% SUB FUNCTIONS

function data_out = remove_across_all (data_in,condition,badtrials,ppparms,faparms)
% reruns the freqency analysis on the data set removing the badtrials
% across all the channels
data_in.epochs = data_in.epochs(condition);
mmil_logstr(parms,'Removing %s of %s trials across all channels.',...
    num2str(length(badtrials)),num2str(data_in.epochs(condition).num_trials));
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
data_in.epochs               = data_in.epochs(condition);
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
   faparms.channel    = 'all';
   warning off MATLAB:NonIntegerInput
   ft_timefreq  = freqanalysis(faparms,ft_timelock,faparms.trials_flag);
   clear ft_timelock
   tmp_timefreq = ts_fieldtrip2timefreq(ft_timefreq,data_in);
   clear ft_timefreq
   [a,includedchans] =intersect({data_out.sensor_info.label}, ft_timefreq.label);
   includedchans     = sort(includedchans);
   %% Sort it into the proper locations
   data_out.timefreq.power(includedchans,:,:)  = tmp_timefreq.timefreq.power(includedchans,:,:);
   data_out.timefreq.data (includedchans,:,:)  = tmp_timefreq.timefreq.data (includedchans,:,:);
   data_out.timefreq.num_trials(includedchans) = data_in.epochs.num_trials - length(badtrialsbygroup);
   clear tmp_data tmp_timefreq
end
[data_out.timefreq.num_trials(find([data_out.sensor_info.badchan]))] = deal(0);

function data_out = remove_by_channel(data_in,condition,data_out_temp,in_ts,in_ft,badtrials_bychan,faparms)
% remove trials by channel - just uses the data_out_tmp
% added conditional statement to save trial data - by Jason Sherfey on 05-Sep-2008
if faparms.trials_flag && length(data_out_temp)==1
	ch =  length(badtrials_bychan); %find(strcmp({data_in.sensor_info.label},data_out_temp.sensor_info.label));
	data_out = data_out_temp;
	data_out.num_sensors = size(data_out_temp.timefreq.power,1);
	data_out.sensor_info = data_in.sensor_info(in_ts(ch));
	if size(data_out_temp.timefreq.power,4) > length(badtrials_bychan{end})
		% not all trials are bad
		data_out.timefreq.power = data_out_temp.timefreq.power(:,:,:,setdiff(1:size(data_out_temp.timefreq.power,4),badtrials_bychan{end}));
		data_out.timefreq.data = data_out_temp.timefreq.data(:,:,:,setdiff(1:size(data_out_temp.timefreq.data,4),badtrials_bychan{end}));		
%		data_out.timefreq.data 			= nan(length(data_in.sensor_info),...
%																				size(data_out_temp(1).timefreq.data,2),...
%																				size(data_out_temp(1).timefreq.data,3),...
%																				size(data_out_temp(1).timefreq.data,4));	
%		data_out.timefreq.power 			= nan(length(data_in.sensor_info),...
%																				size(data_out_temp(1).timefreq.power,2),...
%																				size(data_out_temp(1).timefreq.power,3),...
%																				size(data_out_temp(1).timefreq.power,4));	
%		goodtrials 	 = setdiff(1:size(data_out_temp.timefreq.power,4),badtrials_bychan{end});
%		data_out.timefreq.power(:,:,:,goodtrials)  = data_out_temp.timefreq.power(:,:,:,goodtrials);
%		goodtrials 	 = setdiff(1:size(data_out_temp.timefreq.data,4),badtrials_bychan{end});
%		data_out.timefreq.data(:,:,:,goodtrials) 	 = data_out_temp.timefreq.data(:,:,:,goodtrials);				
		data_out.timefreq.num_trials = data_out.timefreq.num_trials - length(badtrials_bychan{end});
	else
		% all trials are bad
		data_out.timefreq.power(1,:,:,:) = nan(1,size(data_out_temp.timefreq.power,2),size(data_out_temp.timefreq.power,3),size(data_out_temp.timefreq.power,4));
		data_out.timefreq.data(1,:,:,:) = nan(1,size(data_out_temp.timefreq.data,2),size(data_out_temp.timefreq.data,3),size(data_out_temp.timefreq.data,4));		
		data_out.timefreq.num_trials = 0;
	end
else
	data_in.epochs               = data_in.epochs(condition);
	data_out                     = data_out_temp(1);
	data_out.sensor_info         = data_in.sensor_info;
	data_out.timefreq.num_trials = repmat(data_out.timefreq.num_trials,[1 length(data_out.sensor_info)]);
	data_out.timefreq.data       = nan(length(data_in.sensor_info),...
	                                   size(data_out_temp(1).timefreq.power,2),...
	                                   size(data_out_temp(1).timefreq.power,3));
	data_out.timefreq.power      = nan(length(data_in.sensor_info),...
	                                   size(data_out_temp(1).timefreq.data,2),...
	                                   size(data_out_temp(1).timefreq.power,3));
	for ch = 1:length(in_ft)
	    data_out.timefreq.data (in_ts(ch),:,:)   = data_out_temp(in_ft(ch)).timefreq.data (1,:,:);
	    data_out.timefreq.power(in_ts(ch),:,:)   = data_out_temp(in_ft(ch)).timefreq.power(1,:,:);
	    data_out.timefreq.num_trials(in_ts(ch))  = data_out.timefreq.num_trials(in_ts(ch)) - length(badtrials_bychan{in_ft(ch)});
	end
	[data_out.timefreq.num_trials(find([data_out.sensor_info.badchan]))] = deal(0);
end


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

% function added by Jason Sherfey on 05-Sep-2008
function tf_file_name = save_trials(tfchan,parms,faparms,rejparms)
% save trial data for one channel for one condition
if parms.savesingles_flag
    try tfchan.timefreq.power = single(tfchan.timefreq.power); end
    try tfchan.timefreq.cmplx = single(tfchan.timefreq.cmplx); end
    try tfchan.timefreq.data  = single(tfchan.timefreq.data); end
    try tfchan.timefreq.cross = single(tfchan.timefreq.cross); end
end 
if parms.savecomplex_flag
    try tfchan.timefreq = rmfield(tfchan.timefreq,'power'); end
    try tfchan.timefreq = rmfield(tfchan.timefreq,'data'); end
else
    try tfchan.timefreq = rmfield(tfchan.timefreq,'cmplx'); end
%     try tfchan.timefreq = rmfield(tfchan.timefreq,'data'); end
end
% only save 'data' if not all NaNs
if isfield(tfchan.timefreq,'data')
  tmpdat = tfchan.timefreq.data;
  if sum(isnan(tmpdat(:)))~=numel(tmpdat)
    tfchan.timefreq.data = tmpdat; 
  else
    tfchan.timefreq = rmfield(tfchan.timefreq,'data');
  end
end
timefreq_data = tfchan; 
clear tfchan tmpdat
if ~isempty(parms.freq_path)
  outdir = parms.freq_path;
else
  outdir = pwd;
end
if ~exist(outdir,'dir')
  mmil_logstr(parms,'Creating directory: %s.',outdir);
  unix(['mkdir -p ' outdir]);
end

% Set up filename (this should match the format in ts_iEEG_FreqAnalysis)
if isfield(parms,'filename') && ~isempty(parms.filename)
  [pathstr,name,ext] = fileparts(parms.filename{1});
  tf_file_name = name;
else
  [dum, file_tag, ext] = fileparts(parms.f_name{1});

  clear dum ext;
  file_tag = strrep(file_tag,'.epoch','');
  file_name='';
  if strfind(file_tag,'event')
      file_name = file_tag(1:strfind(file_tag,'event')-1);
  elseif strfind(file_tag,'cond')
      file_name = file_tag(1:strfind(file_tag,'cond')-1);
  elseif strfind(file_tag,'epoch_data')
      % filename generated by ts_process_fif_data
      file_name = [file_tag(1:strfind(file_tag,'epoch_data')+9) '.'];
  end
  %file_name = strcat(file_name,...
  %    'foi',strrep(num2str(parms.foi(1)),'.',''),'-',strrep(num2str(parms.foi(end)),'.',''),...
  %    '.toi',strrep(num2str(parms.toi(1)/1000),'.',''),'-',strrep(num2str(parms.toi(end)/1000),'.',''));
  str_foi = sprintf('foi%d-%d',parms.foi(1),parms.foi(end));
  str_toi = sprintf('toi%.3g-%.3g',parms.toi(1)/1000,parms.toi(end)/1000);
  file_name = sprintf('%s%s.%s',file_name,str_foi,str_toi);

  if rejparms.reject_flag
      switch rejparms.reject_exclude
          case 'inall'
              file_name = sprintf('%s.tfrejall%s',file_name,num2str(rejparms.trial_threshold));
          case 'bychannel'
              file_name = sprintf('%s.tfrejbychan%s',file_name,num2str(rejparms.trial_threshold));
          case 'groupchannels'
              file_name = sprintf('%s.tfrejgroups%s',file_name,num2str(rejparms.trial_threshold));
      end
      file_name = sprintf('%s.blthresh%s',file_name,num2str(rejparms.bl_threshold));
      file_name = sprintf('%s.chanthresh%s',file_name,num2str(rejparms.chan_threshold));		
  end
  tf_name = strcat(file_name,'.timefreq.');
  if strfind(file_tag,'event')
    tf_file_name = strcat(tf_name,'event',num2str(timefreq_data.timefreq.event_code));
  else
    tf_file_name = strcat(tf_name,'cond',num2str(timefreq_data.timefreq.event_code));
  end
end
tf_file_name = sprintf('%s_chan%03i.mat',tf_file_name,parms.channels(faparms.channel));%timefreq_data.sensor_info.lognum);
tf_file_name = fullfile(outdir,tf_file_name);

if exist(tf_file_name,'file') && ~parms.overwrite
	mmil_logstr(parms,'Not overwriting trial data file: %s\n',tf_file_name);
else
  mmil_logstr(parms,'Saving trial data for channel %d of %d: %s\n',faparms.channel,length(parms.channels),tf_file_name);
%   fprintf('Saving trial data for channel %d of %d: %s\n',faparms.channel,length(parms.channels),tf_file_name);
  if isfield(timefreq_data,'parms')
    timefreq_data.parms.filename{1} = tf_file_name;
  else
    parms.filename{1} = tf_file_name;
    timefreq_data.parms = parms;    
  end
	save(tf_file_name,'timefreq_data');
end

