function data = ts_process_data(data,varargin)
% Purpose: import, epoch/re-epoch, and/or preprocess continuous or epoched
% data. This function calls ts_load_data, ts_epoch_data, and ts_preproc.
%
% Usage: data = ts_process_data(data,'key1',val1, ... )
%           -OR-
%        data = ts_process_data(datafile,'key1',val1, ... )
%
% Example: cont_data = ts_process_data(datafile,...
%                             'bandpass_flag'   ,1 ,...
%                             'bandpass_low_cf' ,.5,...
%                             'bandpass_high_cf',40,...
%                             'bandpass_low_tb' ,.2,...
%                             'bandpass_high_tb',5 ,...
%                             'verbose'         ,1);
%                                               
% Example: epoch_data = ts_process_data(cont_data,...
%                               'trigchan'      ,'STI 101',...
%                               'baseline_flag' ,1,...
%                               'blcwindow'     ,[.3 0],...
%                               'dsfact'        ,2 ,...
%                               'prestim'       ,.3,...
%                               'poststim'      ,.8,...
%                               'write_events_file',1,...
%                               'rootoutdir'    ,'/home/user/rootoutdir/',...
%                               'prefix'        ,'prefix',...
%                               'verbose'       ,1);
%
% Inputs: 
%     data: timesurfer data struct of continuous or epoch data 
%       -OR-
%     datafile: the full path to a datafile in the .fif, .edf, .eeg, .set,
%               or .mat format ex. '/home/user/subjectid_raw_data.fif' For subjects
%               with multiple runs a cell array of stings is an acceptable output
%               ex. {'/home/user/subjectid_run1_raw_data.fif' '/home/user/subjectid_run2_raw_data.fif'}
%
% Outputs:
%     data: timesurfer data struct of epoch data (if an ev2, evt, .txt
%           events file, sampling array, or trigger channel(s) is given as a parameter, see
%           ts_epoch_data)
%       -OR-
%     data: timesurfer data struct of continous data (if none of the above
%           is given as a parameter)
%
% Parameters: default : options :  description 
%     channels: [] :: vector of channel indices - overrides chantype
%     chantype: [] :: 'all','mag','grad1','grad2','eeg','other','grad','meg' : select which type of channels to analyse
%     chanlabels: [] :: cell array of strings listing channel labels
%     badchanfile: [] :: name of text file containing bad channel labels
%     evntfile: [] :: full filepath to ev2, evt (BESA), or .txt file
%                     containing events information to be used for epoching
%     samples: [] :: an array of sample points to be used for epoching
%     trigchan: [] :: either a amplitude-coded channel label eg. 'STI 101' or an array of
%                     binary trigger channel labels eg. {'STI001' 'STI002' 'STI 003'}
%     events: [] :: vector of event codes of conditions to process
%     lpfilter: 'no' : 'no','yes' : lowpass filter
%     bpfilter: 'no' : 'no','yes' : bandpass filter
%     lnfilter: 'no' : 'no','yes' : linepass filter
%     prestim: .2 :: duration of prestimulus period (sec)
%     poststim: .6 :: duration of poststimulus period (sec)
%     padding: [] : 1.1*(postim-prestim) : time before & after each epoch to 
%                                          include when filtering (sec) 
%                                          Note: by default padding is set to 110% of epoch length 
%     toilim: [-inf inf] :: time limits ([begin end]) in seconds
%     old_event_codes: [] :: vector of old event codes
%     new_event_codes: [] :: vector of new event codes
%                            Note: old_event_codes and new_event_codes
%                            parameters are used in conjection to change
%                            event codes.
%     old_chan_labels: [] :: vector of old channel labels
%     new_chan_labels: [] :: vector of new channel labels
%                            Note: old_chan_labels and new_chan_labels
%                            parameters are used in conjection to change
%                            channel labels
%     stim_delay: [] :: duration of stimulus onset delay after trigger (sec)
%     max_num_trials: Inf :: maximum number of trials per condition
%     blcwindow: [] :: baseline correction window, [begin end] in seconds,
%                      the default is the complete trial ([-inf inf])
%     calc_ncov_flag: 1 : 0,1 : whether to calculate/recalculate noise
%                               covariance matrix after epoching/re-epoching
%     noisewindow: [] :: noise covariance matrix calculation window, [begin end] in seconds,
%                        the default is the blcwindow
%     ncov_ex_evnts: [] :: vector of event codes that should not be used in calculating the noise covariance matrix
%                          e.g. to exclude events with pre-stimulus activity that is not noise -e.g. a button press
%     concat_ncov_flag: 0 : 0,1 : whether to concatenate noise before calculating noise covarience (1) or 
%                                 calculate noise covarience for each trial and then average (0)
%     max_noise_samples: 5000 :: max # of samples used for noise estimate
%                                Note: especially significant if concat_ncov_flag = 1
%     verbose: 0 : 0,1 : whether to produce full output 
%
% For additional data input parameters, see ts_load_data
% For additional epoching parameters, see ts_epoch_data
% For additional prepocessing parameters, see ts_preproc
%
% Created:  05-Jan-2011, Jason Sherfey
% Modified: 14-Jun-2011, Jason Sherfey
% Last Mod: 20-Mar-2012  Steve Deiss

parms = mmil_args2parms( varargin, ...
                         {  'channels',[],[],...
                            'chantype',[],[],...
                            'chanlabels',[],[],...
                            'badchanfile',[],[],...
                            'evntfile',[],[],...
                            'samples',[],[],...
                            'trigchan',[],[],...
                            'events',[],[],...
                            'lpfilter','no',{'yes','no'},...
                            'bpfilter','no',{'yes','no'},...
                            'lnfilter','no',{'yes','no'},...
                            'prestim' ,.2,[],...
                            'poststim',.6,[],...
                            'padding',[],[],...
                            'toilim',[-inf inf],[],...
                            'old_event_codes',[],[],...
                            'new_event_codes',[],[],...
                            'old_chan_labels',[],[],...
                            'new_chan_labels',[],[],...
                            'stim_delay',[],[],...
                            'max_num_trials',Inf,[1 Inf],...
                            'verbose',0,{0,1},...
                            'blcwindow',[],[],...
                            'calc_ncov_flag',true,[false true],...
                            'noisewindow',[],[],...
                            'ncov_ex_evnts',[],[],...
                            'max_noise_samples',5000,[],...
                            'concat_ncov_flag',false,[false true],...
                         }, ...
                         false );

parms    = backcompatible(parms,varargin{:});
if isempty(parms.blcwindow)  , parms.blcwindow   = [-inf inf];      end
if isempty(parms.noisewindow), parms.noisewindow = parms.blcwindow; end
val      = struct2cell(parms);
fld      = fieldnames(parms);
addfld   = setdiff(varargin(1:2:end),fld);
if ~isempty(addfld)
  idx      = 2 * find(ismember(varargin(1:2:end),addfld)) -  1;
  allparms = cell2struct({varargin{idx+1},val{:}},{varargin{idx},fld{:}},2);
else
  allparms = parms;
end

if ~iscell(parms.new_event_codes), parms.new_event_codes = {parms.new_event_codes}; end
if ~iscell(parms.old_event_codes), parms.old_event_codes = {parms.old_event_codes}; end

% Check whether TimeSurfer data or a data file was input
if ischar(data) || (iscell(data) && exist(data{1},'file'))
  % data is a file that needs to be loaded
  datafile = data; clear data
  if ~iscell(datafile), datafile = {datafile}; end
  % we don't want to remove channels based on indices at this step
  tmpparms = rmfield(allparms,{'channels','chantype','chanlabels','badchanfile'});
  tmpargs  = mmil_parms2args(tmpparms);
  % load data
  for i = 1:length(datafile)
    data(i) = ts_load_data(datafile{i},tmpargs{:});
    if length(parms.old_event_codes)==length(datafile) && length(parms.new_event_codes)==length(datafile)
      % rename event codes
      [c,ia,ib] = intersect([data(i).epochs.event_code],parms.old_event_codes{i});
      if ~isempty(c)
        newcodes  = num2cell(parms.new_event_codes{i}(ib));
        [data(i).epochs(ia).event_code] = deal(newcodes{:});
      end
    end
  end
  clear tmpparms tmpargs
  % NOTE: do not concatenate yet; we do not know if data is continuous or 
  % epoched already at this point.
elseif isnumeric(data)
  data     = ts_matrix2data(data);
  datafile = {[]};
else
  datafile = {[]};
% ---------------
  if isstruct(data)
    for i = 1:length(data)
      if length(parms.old_event_codes)==length(data) && length(parms.new_event_codes)==length(data)
        % rename event codes
        [c,ia,ib] = intersect([data(i).epochs.event_code],parms.old_event_codes{i});
        if ~isempty(c)
          newcodes  = num2cell(parms.new_event_codes{i}(ib));
          [data(i).epochs(ia).event_code] = deal(newcodes{:});
        end
      end
    end
  end
% ---------------
end

% Check that all data structures are the same type
for i = 1:length(data)
  if i == 1, type = ts_object_info(data(i),'verbose',0); end
  if ~isequal(type,ts_object_info(data(i),'verbose',0))
    error('All input data must be of the same type (epochs, avgs, TFR).\n');
    return;
  end
end
[datatype,datafield,dataparam] = ts_object_info(data(1));

%%  SRD 03202012, This section is the old way to do stim_delay and was
%%  replaced by an adjustment to ts_epoch_data.  
% % shift time vector (do this before any time-dependent processing)
% %  SIMPLY MODIFIES THE TIME VECTOR ACROSS THE CONDITIONS SO THAT THE 0 TIME
% %  POINT NOW CORRESPONDS TO THE STIMULUS RATHER THAN TRIGGER ONSET
% if ~isempty(parms.stim_delay)
%   mmil_logstr(parms,'Introducing stimulus delay of %s seconds.\n',num2str(parms.stim_delay))
%   fprintf('Introducing stimulus delay of %s seconds.\n',num2str(parms.stim_delay))
% %   sampling_rate = data(k).sfreq;  %  This is in the wrong place, see inside the loop. SRD
% %   time_steps    = 1/sampling_rate;
%   for k = 1:length(data)
%     sampling_rate = data(k).sfreq;
%     time_steps    = 1/sampling_rate;
%     for j = 1:length(data(k).epochs)
%       orig_time     = data(k).epochs(j).time;
%       prestim_samp  = find(orig_time <= parms.stim_delay,1,'last')-1;
%       poststim_samp = length(orig_time) - (prestim_samp + 1);
%       time_min      = -(prestim_samp*time_steps);
%       time_max      =  (poststim_samp*time_steps);
%       parms.prestim = time_min;
%       parms.poststim = time_max;
%       data(k).epochs(j).time = [];
%       data(k).epochs(j).time = [time_min:time_steps:time_max];
%       mmil_logstr(parms,'Condition %2s, prestimulus period: %s sec / poststimulus period: %s sec.',num2str(j),num2str(time_min),num2str(time_max));
% %       fprintf('Condition %2s, prestimulus period: %s sec / poststimulus period: %s sec.\n',num2str(j),num2str(time_min),num2str(time_max));
%     end
%   end
%   clear orig_time prestim_samp poststim_samp time_min time_max sampling_rate time_steps
% end
%%

% set flags
filter_flag     = strcmpi(parms.bpfilter,'yes') || strcmpi(parms.lnfilter,'yes') || strcmpi(parms.lpfilter,'yes');
% assume data are continuous if all num_trials = 1
input_cont_flag = all(arrayfun(@(x)all([x.(datafield).num_trials]==1),data));
% set flags (cont_flag - whether to OUTPUT continuous data, filter_flag)
trigchan_flag   = ~isempty(parms.trigchan) && input_cont_flag && any(ismember(parms.trigchan,{data(1).sensor_info.label}));
  % whether a trigger channel should be used to epoch continuous data
evntfile_flag   = ~isempty(parms.evntfile);% && exist(parms.evntfile,'file');
  % whether an event file should be used to epoch continuous or epoched data
samples_flag    = (isnumeric(parms.samples) && ~isempty(parms.samples));
  % whether an array of sample points should be used for epoching
% epoch_flag      = ~(input_cont_flag && (trigchan_flag || evntfile_flag || samples_flag));
epoch_flag      = trigchan_flag || evntfile_flag || samples_flag;
  % whether to return epochs or not
dsfact_ix       = find(cellfun(@(x)isequal(x,'dsfact'),varargin));
downsample_flag = ~isempty(dsfact_ix) && ~isempty(varargin{dsfact_ix+1}) && varargin{dsfact_ix+1}>1;
if downsample_flag
  downsample_fact = varargin{dsfact_ix+1}; 
end

if epoch_flag
  evntfile = parms.evntfile;
  if ~iscell(evntfile)
    evntfile = {evntfile};
  end
  if length(evntfile) ~= length(datafile)
    evntfile(1:length(datafile)) = evntfile;
  end
end

% filter padding
if isempty(parms.padding)
  parms.padding = (parms.prestim + parms.poststim)*1.1;
end

% Process data
if epoch_flag   % epoch/re-epoch the data
  % prepare options for epoching
  tmpcfg  = rmfield(allparms,{'prestim','poststim','toilim','evntfile','events',...
                    'old_event_codes','new_event_codes','calc_ncov_flag',...
                    'channels','chantype','chanlabels','badchanfile'});
  tmpargs = mmil_parms2args(tmpcfg);
  old_event_codes = [];
  new_event_codes = [];
  for i = 1:length(data)
    % filter padding and time limits
    T            = data(i).(datafield)(1).time;
    parms.toilim = [max(T(1),parms.toilim(1)) min(T(end),parms.toilim(2))];
    if filter_flag
      filtertime  = [-parms.prestim-parms.padding parms.poststim+parms.padding];
      timelimits  = [-parms.prestim parms.poststim];
    else
      filtertime = [-parms.prestim parms.poststim];
      timelimits = filtertime;
    end
    % NOTE: do not select channels before epoching b/c we might need
    % trigger channel for epoching
%     if ~isempty(parms.channels) || ~isempty(parms.chantype) || ~isempty(parms.chanlabels) || ~isempty(parms.badchanfile)
%       % select channels before epoching
%       data(i) = ts_data_selection(data(i),'channels',parms.channels,'chantype',parms.chantype,'badchanfile',parms.badchanfile,'chanlabels',parms.chanlabels);
%     end
    if length(parms.old_event_codes)==length(data) && length(parms.new_event_codes)==length(data)
      old_event_codes = parms.old_event_codes{i};
      new_event_codes = parms.new_event_codes{i};
    end    
    % epoch the data (wrt events_fnames, ev2, evt, trigchans, samples, ...)
    data(i) = ts_epoch_data(data(i),'prestim',-filtertime(1),'poststim',filtertime(2),'calc_ncov_flag',false,... 
      'datafile',datafile{i},'evntfile',evntfile{i},'old_event_codes',old_event_codes,'new_event_codes',new_event_codes,'events',parms.events,tmpargs{:});
    % select subset of channels before preprocessing
    if ~isempty(parms.channels) || ~isempty(parms.chantype) || ~isempty(parms.chanlabels) || ~isempty(parms.badchanfile)
      % select channels before epoching
      data(i) = ts_data_selection(data(i),'channels',parms.channels,'chantype',parms.chantype,'badchanfile',parms.badchanfile,'chanlabels',parms.chanlabels);
    end
    % preprocess (filter, downsample, detrend, baseline correct)
    data(i) = ts_preproc(data(i),tmpargs{:});
    % remove padding and apply selection criteria
    data(i) = ts_data_selection(data(i),'toilim',timelimits,tmpargs{:});
  end
  clear tmpcfg tmpargs
  if i > 1
    data = ts_combine_data(data); % concatenate trials & trial_info fields
  end
  if parms.calc_ncov_flag
    data = ts_calc_ncov(data,'noisewindow',parms.noisewindow,'ncov_ex_evnts',parms.ncov_ex_evnts,...
                        'concat_ncov_flag',parms.concat_ncov_flag,'max_noise_samples',parms.max_noise_samples); 
  end
else            % do not epoch or re-epoch the data
  tmpcfg  = rmfield(allparms,{'toilim','channels','chantype','chanlabels','badchanfile'});
  tmpargs = mmil_parms2args(tmpcfg);
  for i = 1:length(data)
    % filter padding and time limits
    T            = data(i).(datafield)(1).time;
    parms.toilim = allparms.toilim;
    parms.toilim = [max(T(1),parms.toilim(1)) min(T(end),parms.toilim(2))];
    if filter_flag
      filtertime  = [max(T(1),parms.toilim(1)-parms.padding) min(T(end),parms.toilim(2)+parms.padding)];
      timelimits  = parms.toilim;
    else
      filtertime = parms.toilim;
      timelimits = parms.toilim;
    end
    % select channels before preprocessing continuous data
    if ~isempty(parms.channels) || ~isempty(parms.chantype) || ~isempty(parms.chanlabels) || ~isempty(parms.badchanfile)
      data(i) = orderfields(ts_data_selection(data(i),'channels',parms.channels,'chantype',parms.chantype,'badchanfile',parms.badchanfile,'chanlabels',parms.chanlabels),data(i));
    end
    % select time limits + padding for filtering
    data(i) = orderfields(ts_data_selection(data(i),'toilim',filtertime,tmpargs{:}),data(i));
    % hold unfiltered trigger data
    sti_i   = strncmp('STI',{data(i).sensor_info.label},3);
    % have epoch data that is not going to be re-epoched:
    for c = 1:length(data(i).epochs)
      if isempty(find(sti_i)), break; end
      sti_x{c} = data(i).epochs(c).data(sti_i,:);
      if downsample_flag
        sti_x{c} = downsample(sti_x{c}',downsample_fact)';
      end
    end   
    % filter all channels
    data(i) = ts_preproc(data(i),tmpargs{:});
    % replace unfiltered trigger data
    for c = 1:length(data(i).epochs)
      if isempty(find(sti_i)), break; end
      data(i).epochs(c).data(sti_i,:) = sti_x{c};
    end
    clear sti_x
    % ---------------------------------------
    % ---------------------------------------
    % select time you want to return
    data(i) = ts_data_selection(data(i),'toilim',timelimits);
    % concatenate data structures in time
    if i == 1
      alldata = data(i);
    elseif input_cont_flag
      alltime = alldata.epochs(1).time;
      curtime = data(i).epochs(1).time - min(data(i).epochs(1).time);
        % NOTE: this shifts all time vectors after the first to zero so
        % that there are no gaps between sequential files even when the raw
        % data time vector starts at t > 0.
      alldata.epochs.time = [alltime curtime+alltime(end)+(1/alldata.sfreq)];
      alldata.epochs.data = cat(2,alldata.epochs.data,data(i).epochs.data);
    else
      alldata(i) = data(i);
    end
    [data(i).epochs.data] = deal([]);
    clear alltime curtime
  end
  clear data
  if ~input_cont_flag && i>1
    alldata = ts_combine_data(alldata);
  end
  data = alldata;
  clear alldata alldatatmp
end

% rename sensor labels
if iscellstr(parms.new_chan_labels)
  [sel1,sel2] = match_str({data.sensor_info.label},parms.old_chan_labels);
  [data.sensor_info(sel1).label] = deal(parms.new_chan_labels{:});
  clear sel1 sel2
end

data.parms = allparms;

%% SUBFUNCTIONS
function parms = backcompatible(parms,varargin)
  opt = mmil_args2parms(varargin,{...
          'event',[],[],...
          'event_codes',[],[],...
          'valid_event_codes',[],[]...
          'eventvals',[],[],...
          'events_fnames',[],[],...
          'event_fnames',[],[],...
          'trigchans',[],[],...
          'bandpass_flag',0,[],...
          'notch_flag',0,[],...
          'lowpass_flag',0,[],...
          'noise_start',[],[],...
          'noise_end',[],[],...
          'baseline_start',[],[],...
          'baseline_end',[],[],...
          },false);      
        
  if isempty(parms.events)
    if ~isempty(opt.event)
      parms.events = opt.event;
    elseif ~isempty(opt.event_codes)
      parms.events = opt.event_codes;
    elseif ~isempty(opt.eventvals)
      parms.events = opt.eventvals;
    elseif ~isempty(opt.valid_event_codes)
      parms.events = opt.valid_event_codes;
    end
  end
  if isempty(parms.evntfile) && ~isempty(opt.events_fnames)
    parms.evntfile = opt.events_fnames;
  elseif isempty(parms.evntfile) && ~isempty(opt.event_fnames)
    parms.evntfile = opt.event_fnames; 
  end
  if isempty(parms.trigchan) && ~isempty(opt.trigchans)
    parms.trigchan = opt.trigchans;
  end
  if strcmp(parms.bpfilter,'no') && opt.bandpass_flag
    parms.bpfilter = 'yes';
  end
  if strcmp(parms.lnfilter,'no') && opt.notch_flag
    parms.bpfilter = 'yes';
  end
  if strcmp(parms.lpfilter,'no') && opt.lowpass_flag
    parms.bpfilter = 'yes';
  end
  if isempty(parms.noisewindow) && ~isempty(opt.noise_start) && ~isempty(opt.noise_end)
    parms.noisewindow = [opt.noise_start opt.noise_end];
  end
  if isempty(parms.blcwindow) && ~isempty(opt.baseline_start) && ~isempty(opt.baseline_end)
    parms.blcwindow = [opt.baseline_start opt.baseline_end];
  end
  