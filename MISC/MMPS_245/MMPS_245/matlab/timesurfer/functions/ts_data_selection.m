function [outdata,opt,err] = ts_data_selection(data,varargin)
% Purpose: select a subset of data, one chan choose ranges to keep by
%          selecting  conditions, channels, trials, toi, or foi; or one can
%          use rejection information to remove rejects 
%
% Usage: data = ts_process_data(data,'key1',val1, ... )
%
% Example: epoch_data = ts_data_selection(epoch_data,...
%                                 'chantype'          ,'grad',...
%                                 'events'            ,[1 2 3],...
%                                 'keepbadtrials_flag',0,...
%                                 'keepbadchans_flag' ,0,...
%                                 'verbose'           ,1);
%
% Inputs:
%     data: timesurfer data struct of continuous, epoch, or average data
%
% Outputs:
%     data: timesurfer data struct of continuous, epoch, or average data
%           containing only the selected data
%
% Parameters: default : options :  description 
%     condition: [] :: condition number (not event code) used to index data inside structure
%     events: [] :: array of eventcodes to keep
%     evntfile: [] :: full filepath to a events .txt file 
%     samples: [] :: array of samples to keep
%     channels: [] :: array of channel indices to keep - overrides chantype
%     chantype: 'all' : 'all','mag','grad1','grad2','eeg','other','grad','meg' : select which type of channels keep (overrides 'channels' parameter)
%     chanlabels: [] :: cell array of strings listing channel labels
%     badchans: [] :: array of bad channel indices
%     badlabels: [] :: arral of bad channel labels
%     keepbadchans_flag: 0 : 0,1 : whether to keep channels previously marked as bad
%     keepbadtrials_flag: 0 : 0,1 : whether to keep trials previously marked as bad
%     markbadchans_flag: 0 : 0,1 : if the keepbadchans_flag = 1, wether to keep bad channels marked as bad, 
%                                  if markbadchans_flag = 0 & keepbadchans_flag = 1 all channels will be set to good,
%                                  if keepbadchans_flag = 0,
%                                  markbadchans_flag is ignored
%     trials:  [] :: array of trials to keep
%     toilim: [] :: limits on the times of interest, [begin end] in seconds
%     toi: [] :: array of times to keep
%     foilim: [] :: limits on the frequencies of interest, [begin end] in seconds
%     foi: [] :: array of frequencies to keep
%     badchanfile: [] :: full file path to a .txt file containing a line
%                        delimited line of bad channel labels
%     rejectfile: [] : full filepath to a timesurfer .mat file containing reject data
%     reject_data: [] : timesurfer reject data structure
%     reject_data_latflag: 0 : 0,1 : whether to use the latency rather than
%                                    the trial index when applying reject data
%     sensor_info: [] : sensor info, by default this is read from the
%     timesurfer data struct
%     calc_ncov_flag: 0 : 0,1 : whether to calculate/recalculate noise
%                               covariance matrix after
%                               epoching/re-epoching
% last modified by BQR on 08/13/2012



parms = mmil_args2parms(varargin,{...
        'condition',[],[],...
        'event',[],[],...
        'evntfile',[],[],...
        'samples',[],[],...        
        'channel',[],[],...
        'channels',[],[],...
        'chantype','all', [],...
        'chanlabel',[],[],...
        'badchans',[],[],...
        'badlabels',[],[],...
        'keepbadchans',[],[],...
        'keepbadtrials',0,{0,1},...
        'markbadchans_flag',0,{0,1},...
        'reject_data_latflag',0,{0,1},...
        'removebadchans',[],[],...data.sensor_info
        'trial',[],[],...
        'rejects',[],[],...
        'toilim',[],[],...
        'toi',[],[],...
        'foilim',[],[],...
        'foi',[],[],...
        'badchanfile',[],[],...
        'rejectfile',[],[],...
        'reject_data',[],[],...
        'sensor_info',[],[],...
        'opt',[],[],...
        'dataparam',[],[],...
        'verbose',1,{0,1},...
        'logfile',      [],[],...
        'logfid',       [1],[], ...        
        'calc_ncov_flag',false,[false true],...
        },false);
      
parms = backcompatible(parms,varargin{:});     
data  = ts_checkdata_header(data,'events',parms.event);
if isempty(parms.sensor_info), parms.sensor_info = data.sensor_info; end

% Validate inputs
% note: only keepbadchans_flag would exist in an ideal world.  poor
% development is the only reason keepbadchans, keepbadchans_flag, &
% removebadchans all exist in this function.
if ~isempty(parms.keepbadchans) && ~isempty(parms.removebadchans) && parms.keepbadchans==parms.removebadchans
  error('you cannot set removebadchans equal to keepbadchans or keepbadchans_flag');
end
if isempty(parms.keepbadchans) && ~isempty(parms.removebadchans), parms.keepbadchans = ~parms.removebadchans; end
if ~isempty(parms.keepbadchans) && isempty(parms.removebadchans), parms.removebadchans = ~parms.keepbadchans; end
if isempty(parms.keepbadchans) && isempty(parms.removebadchans)
  parms.keepbadchans = 0;
  parms.removebadchans = 1;
end

err = 0; opt = [];
[object,datafield,dataparam] = ts_object_info(data,varargin{:});

if ~isempty(parms.dataparam)
  if all(~ismember(parms.dataparam,dataparam))
    parms.dataparam = dataparam;
  else
    dataparam = parms.dataparam; 
  end
end  
if ~iscell(dataparam), dataparam = {dataparam}; end
if ~isempty(parms.chanlabel) && ~iscell(parms.chanlabel), parms.chanlabel = {parms.chanlabel}; end

% check trial_info and try to fix organization if necessary
if isfield(data.(datafield),'trial_info'), check_trial_info; end

%% conditions
if isempty(parms.event) && isempty(parms.condition)
  cond = 1:length(data.(datafield));
elseif ~isempty(parms.condition)
  cond = parms.condition;
elseif ~isempty(parms.event)
  if iscell(parms.event)
    parms.event = unique([parms.event{:}]);
  end
  cond = find(ismember([data.(datafield).event_code],parms.event));
end
events = [data.(datafield)(cond).event_code];

%% rejects
badchans                      = parms.badchans;
[badtrials{1:length(cond)}]   = deal([]);

% reject matfile with reject_data
if ~isempty(parms.rejectfile) && exist(parms.rejectfile,'file')
  mmil_logstr(parms,'loading reject file: %s',parms.rejectfile);
  reject_data = []; %added 8/22/11 BQR
  load(parms.rejectfile);
  parms.reject_data = reject_data;
  if isfield(parms.reject_data,'badchanlabels')
    [sel1 sel2] = match_str({parms.sensor_info.label},parms.reject_data.badchanlabels);
    parms.reject_data.badchans = sel1;
  elseif ~isempty(parms.reject_data.badchans)
    mmil_logstr(parms,'warning: the following chans will be rejected based on indices, not labels:');
    {parms.sensor_info(parms.reject_data.badchans).label}
  end
  badchans = union(badchans,parms.reject_data.badchans);
  if isfield(parms.reject_data,'event_code')
    for k = 1:length(parms.reject_data.badtrials)
      if ismember(parms.reject_data.event_code(k),events)
        badtrials{parms.reject_data.event_code(k)==events} = parms.reject_data.badtrials{k};
      end
    end
  elseif length(parms.reject_data.badtrials) == length(events)
    [badtrials{:}]  = deal(parms.reject_data.badtrials{cond});
  end
elseif isstruct(parms.reject_data)
  if isfield(parms.reject_data,'badchanlabels')
    [sel1 sel2] = match_str({parms.sensor_info.label},parms.reject_data.badchanlabels);
    parms.reject_data.badchans = sel1;
  elseif ~isempty(parms.reject_data.badchans)
    mmil_logstr(parms,'warning: the following chans will be rejected based on indices, not labels:');
    {parms.sensor_info(parms.reject_data.badchans).label}
  end
  badchans        = union(badchans,parms.reject_data.badchans);
  if isfield(parms.reject_data,'event_code')
    for k = 1:length(parms.reject_data.badtrials)
      if ismember(parms.reject_data.event_code(k),events)
          if parms.reject_data_latflag
                badtrials{parms.reject_data.event_code(k)==events} = ...
                    data.(datafield)(k).trial_info.number(ismember(...
                            data.(datafield)(k).trial_info.latency,parms.reject_data.badtrial_info(k).latency));
          else
                badtrials{parms.reject_data.event_code(k)==events} = parms.reject_data.badtrials{k};
          end
      end
    end
  elseif length(parms.reject_data.badtrials) == length(events)
    [badtrials{:}]  = deal(parms.reject_data.badtrials{cond});
  end
end

% ascii file listing bad channel labels
if ~isempty(parms.badchanfile) && exist(parms.badchanfile,'file')
  badchans = union(badchans,ts_read_txt_badchans(parms.badchanfile,{parms.sensor_info.label}));
end

% bad channels flagged in sensor_info
badchans = union(badchans,find([parms.sensor_info.badchan]));

% bad channel labels
if ~isempty(parms.badlabels) 
  if iscellstr(parms.badlabels)
    [sel,jnk] = match_str({parms.sensor_info.label},parms.badlabels);
  elseif ischar(parms.badlabels)
    sel = strmatch(parms.badlabels,{parms.sensor_info.label});
  else
    sel = [];
  end
  badchans = union(badchans,sel);
end
badchans = sort(unique(badchans));

if parms.keepbadchans
	[data.sensor_info.badchan]  = deal(0);
    [parms.sensor_info.badchan] = deal(0); 
    if parms.markbadchans_flag %added 2/02/12 BQR
        [data.sensor_info(badchans).badchan]  = deal(1);
        [parms.sensor_info(badchans).badchan] = deal(1);
    end
    badchans = [];
end


%% channels
chans = [];
if ~isempty(parms.chanlabel) && iscell(parms.chanlabel)
  [sel1 chans] = match_str(parms.chanlabel,{parms.sensor_info.label});
elseif ~isempty(parms.channel)
  chans = parms.channel;  
elseif ~isempty(parms.channels)
  chans = parms.channels;
elseif ~isempty(parms.chantype)
  if ~isempty(parms.sensor_info)
    if ~iscell(parms.chantype), parms.chantype = {parms.chantype}; end
    for k = 1:length(parms.chantype)      
      switch lower(parms.chantype{k})
        case {'mag' 'grad1' 'grad2' 'eeg', 'other','sti','eog','ekg'}
          chans     = [chans find(strcmp(parms.chantype{k},{parms.sensor_info.typestring}))];
        case 'ieeg'
          chans     = [chans find(strcmp('eeg',{parms.sensor_info.typestring}))];
          if isempty(chans)
            chans   = [chans find(strcmp('ieeg',{parms.sensor_info.typestring}))];
          end
        case {'grad'}
          chans     = [chans find(strncmp(parms.chantype{k},{parms.sensor_info.typestring},length(parms.chantype{k})))];
        case 'meg'
          [a,chans] = [chans find(ismember({parms.sensor_info.typestring},{'mag', 'grad1', 'grad2'}))];
        case 'all'
          try chans = [chans setdiff(1:data.num_sensors,find(strcmp('other',{parms.sensor_info.typestring})))]; end;
        otherwise
          chans     = [chans find(strcmp(parms.chantype{k},{parms.sensor_info.typestring}))];
      end
    end
  end
else
  chans = 1:length(parms.sensor_info);
end
% [badchans sel2] = match_str({parms.sensor_info(chans).label},{parms.sensor_info(badchans).label});
badchans = sort(intersect(badchans,chans));
chans    = sort(setdiff(chans,badchans));
badchans = sort(setdiff(1:data.num_sensors,chans));
if any(ismember(badchans,parms.channels))
  warning('You requested "channels" that are marked as bad.');
end
if isempty(chans)
  mmil_logstr(parms,'no good channels selected');
  err = 1; outdata = [];
  return;
end;


for i = 1:length(cond)
  %% trials
  if ~isempty(parms.trial)
    if iscell(parms.trial) && length(parms.trial)==length(cond)
      trials{i} = parms.trial{i};
    elseif isnumeric(parms.trial)
      trials{i} = parms.trial;
    end
  elseif isfield(data.(datafield),'num_trials')
    trials{i} = 1:data.(datafield)(cond(i)).num_trials;
  else
    trials{i} = 1;
  end
  if isfield(data.(datafield),'trial_info') && ~parms.keepbadtrials % added 25-Jan-2011 by JSS
    badtrials{i} = [badtrials{i} find([data.(datafield)(cond(i)).trial_info.badtrial]==1)];
  end
  trials{i} = sort(setdiff(trials{i},badtrials{i}));
  if isempty(trials{i})
    mmil_logstr(parms,'no good trials selected for event %g',data.(datafield)(cond(i)).event_code);
    trials{i} = [];
%     err = 1;  outdata = [];
%     return;
  end
  trialvals(i) = length(trials{i});
  if strcmp(datafield,'averages') || ~isfield(data.(datafield),'num_trials') || ~ismember(data.(datafield)(cond(i)).num_trials,size(data.(datafield)(cond(i)).(dataparam{1}))) ...
      || (isfield(data.(datafield)(cond(i)),'frequencies') && ndims(data.(datafield)(cond(i)).(dataparam{1})) < 4)
    trials{i} = 1;  % this is an average
    isavg_flag = 1;
    badtrials{i} = [];
  else
    isavg_flag = 0;
  end
  %% time
  if ~isempty(parms.toi)
    [c tidx{i} ib] = intersect(data.(datafield)(cond(1)).time,parms.toi);
  else
    if isempty(parms.toilim) || length(parms.toilim)~=2
      parms.toilim = [data.(datafield)(cond(i)).time(1) data.(datafield)(cond(i)).time(end)];
    end
    tidx{i} = find(data.(datafield)(cond(i)).time>=parms.toilim(1) & ...
                   data.(datafield)(cond(i)).time<=parms.toilim(end));
%     tidx{i} = nearest(data.(datafield)(cond(i)).time,parms.toilim(1)):...
%               nearest(data.(datafield)(cond(i)).time,parms.toilim(end));  
  end
  %% frequency     
  freq = 0;
  if isfield(data.(datafield),'frequencies')
    freq = 1;
    if ~isempty(parms.foi)
      [c fidx{i} ib] = intersect(data.(datafield)(cond(i)).frequencies,parms.foi);
    else
      if isempty(parms.foilim) || length(parms.foilim)~=2
        parms.foilim = [data.(datafield)(cond(i)).frequencies(1) data.(datafield)(cond(i)).frequencies(end)];
      end
      fidx{i} = find(data.(datafield)(cond(i)).frequencies>=parms.foilim(1) & ...
                     data.(datafield)(cond(i)).frequencies<=parms.foilim(end));
%       fidx{i} = nearest(data.(datafield)(cond(i)).frequencies,parms.foilim(1)):...
%                 nearest(data.(datafield)(cond(i)).frequencies,parms.foilim(end));          
    end
  end
end;
for k = 1:length(badtrials)
  badtrials{k} = unique(badtrials{k});
end

%% correct channels for PLV and other 4D data sets
% chans will be for indices into the data matrix
% sens will be for indices into the sensor_info array
senschans    = chans;
sensbadchans = badchans;
if issubfield(data,'timefreq.labelcmb')
  sen = {data.sensor_info(senschans).label};
  cmb = data.timefreq(1).labelcmb;
  keep=[];
  for i = 1:length(sen)
    keep = [keep find(ismember(cmb(:,1),sen{i}) | ismember(cmb(:,2),sen{i}))];
  end
  chans = unique(keep);
  sen = {data.sensor_info(sensbadchans).label};
  cmb = data.timefreq(1).labelcmb;
  keep=[];
  for i = 1:length(sen)
    keep = [keep find(ismember(cmb(:,1),sen{i}) | ismember(cmb(:,2),sen{i}))];
  end
  badchans = unique(keep);
end

%% data selection

outdata = rmfield(data,datafield);
outdata.(datafield) = data.(datafield)(cond);
for i = 1:length(cond)
  for j = 1:length(dataparam)
    dat = dataparam{j};
    if ~isfield(data.(datafield)(cond(i)),dat) || isempty(data.(datafield)(cond(i)).(dat)), continue; end
    if freq == 0
      if parms.removebadchans
        outdata.(datafield)(i).(dat)        = data.(datafield)(cond(i)).(dat)(chans,tidx{i},trials{i});
      else
        outdata.(datafield)(i).(dat)        = data.(datafield)(cond(i)).(dat)(:,tidx{i},trials{i});
      end
    else
      if parms.removebadchans
        outdata.(datafield)(i).(dat)        = data.(datafield)(cond(i)).(dat)(chans,tidx{i},fidx{i},trials{i});
      else
        outdata.(datafield)(i).(dat)        = data.(datafield)(cond(i)).(dat)(:,tidx{i},fidx{i},trials{i});
      end
      outdata.(datafield)(i).frequencies  = data.(datafield)(cond(i)).frequencies(fidx{i});
    end
    if ~parms.keepbadtrials && isfield(outdata.(datafield),'trial_info') && j==1 && ~isavg_flag
      if isfield(outdata.(datafield)(i).trial_info,'number')       , outdata.(datafield)(i).trial_info.number        = outdata.(datafield)(i).trial_info.number(trials{i});     end
      if isfield(outdata.(datafield)(i).trial_info,'latency')      , outdata.(datafield)(i).trial_info.latency       = outdata.(datafield)(i).trial_info.latency(trials{i});    end
      if isfield(outdata.(datafield)(i).trial_info,'badtrial')     , outdata.(datafield)(i).trial_info.badtrial      = outdata.(datafield)(i).trial_info.badtrial(trials{i});   end
      if isfield(outdata.(datafield)(i).trial_info,'event_code')   , outdata.(datafield)(i).trial_info.event_code    = outdata.(datafield)(i).trial_info.event_code(trials{i}); end
      if isfield(outdata.(datafield)(i).trial_info,'duration')     , outdata.(datafield)(i).trial_info.duration      = outdata.(datafield)(i).trial_info.duration(trials{i});   end
      if isfield(outdata.(datafield)(i).trial_info,'datafile')     , outdata.(datafield)(i).trial_info.datafile      = outdata.(datafield)(i).trial_info.datafile(trials{i});   end
      if isfield(outdata.(datafield)(i).trial_info,'events_fnames'), outdata.(datafield)(i).trial_info.events_fnames = outdata.(datafield)(i).trial_info.events_fnames(trials{i}); end
    end
    if islogical(outdata.(datafield)(i).(dat)) && ~parms.removebadchans      % mask
      outdata.(datafield)(i).(dat)(badchans,:,:,:) = 0;
    elseif ~parms.removebadchans
      if freq == 1 % TF data
        outdata.(datafield)(i).(dat)(badchans,:,:,:) = nan;
      else
        outdata.(datafield)(i).(dat)(badchans,:,:,:) = 0;
      end
    end
  end
  outdata.(datafield)(i).num_trials       = trialvals(i);
  outdata.(datafield)(i).time             = data.(datafield)(cond(i)).time(tidx{i});
  if parms.removebadchans
    outdata.sensor_info = parms.sensor_info(senschans);
    outdata.num_sensors = length(senschans);
  else
    [outdata.sensor_info(sensbadchans).badchan] = deal([1]);
  end
  if ~isempty(badtrials{i})
    outdata.(datafield)(i).num_rejects.manual = outdata.(datafield)(i).num_rejects.manual + length(badtrials{i});
    if parms.verbose
      mmil_logstr(parms,'event %g: %g trials removed (%s)',data.(datafield)(cond(i)).event_code,length(badtrials{i}),num2str(badtrials{i}));
    end
  end
  if issubfield(data,'timefreq.labelcmb')
    outdata.timefreq(i).labelcmb = data.timefreq(i).labelcmb(chans,:);
  end
end
if ~isempty(sensbadchans)
  labels = {parms.sensor_info(sensbadchans).label};
  str    = labels{1}; 
  for k=2:length(labels), str=[str ' ' labels{k}]; end
  if parms.removebadchans
    if parms.verbose
      mmil_logstr(parms,'%g channels removed (%s)',length(sensbadchans),str);
    end
  else
    if parms.verbose
      mmil_logstr(parms,'%g channels set to zero and marked bad (%s)',length(sensbadchans),str);
    end
  end
end

% epoch or re-epoch data based on evntfile or vector of sample points
if ~isempty(parms.evntfile) || ~isempty(parms.samples)
  outdata = ts_epoch_data(outdata,varargin{:});
elseif parms.calc_ncov_flag
  % calculate noise covariance matrix
   outdata = ts_calc_ncov(outdata,varargin{:}); 
end

% update options structure
if ~isempty(parms.opt)
  opt = parms.opt;
%   opt.datafile = opt.datafile(cond);
end

outdata = orderfields(outdata,data);

  % nested function
  function check_trial_info
    for cnum = 1:length(data.(datafield))
      trl = data.(datafield)(cnum).trial_info;
      if length(trl) > 1
        data.(datafield)(cnum).trial_info = [];
        if isfield(trl,'number')        , data.(datafield)(cnum).trial_info.number     = [trl.number]; end
        if isfield(trl,'latency')       , data.(datafield)(cnum).trial_info.latency    = [trl.latency]; end
        if isfield(trl,'badtrial')      , data.(datafield)(cnum).trial_info.badtrial   = [trl.badtrial]; end
        if isfield(trl,'event_code')    , data.(datafield)(cnum).trial_info.event_code = [trl.event_code]; end
        if isfield(trl,'duration')      , data.(datafield)(cnum).trial_info.duration   = [trl.duration]; end
        if isfield(trl,'datafile')      , data.(datafield)(cnum).trial_info.datafile   = {trl.datafile}; end
        if isfield(trl,'events_fnames') , data.(datafield)(cnum).trial_info.events_fnames = {trl.events_fnames}; end
        if isfield(trl,'epochnum')      , data.(datafield)(cnum).trial_info.epochnum   = [trl.epochnum]; end
%         data.(datafield)(cnum).trial_info = trl;
      end
    end
  end
end
% subfunction
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
        'conditions',[],[],...
        'cond',[],[],...
        'events',[],[],...
        'event_codes',[],[],...
        'eventvals',[],[],...
        'channels',[],[],...
        'chanlabels',[],[],...
        'chantype','all',[],...
        'badchans',[],[],...
        'c',[],[],...
        'keepbadtrials_flag',[],[],...
        'keepbadchans_flag',[],[],...
        'badlabel',[],[],...
        'trial',[],[],...
        'trials',[],[],...
        'toilim',[],[],...
        'toi',[],[],...
        'foilim',[],[],...
        'foi',[],[],...
        'badchanfile',[],[],...
        'reject_file',[],[],...
        'opt',[],[],...
        'dataparam',[],[],...
        },false);

if ~isempty(opt.cond) && isempty(opt.conditions)
  opt.conditions = opt.cond;
end
if isempty(parms.condition) && ~isempty(opt.conditions)
  parms.condition = opt.conditions; 
end

if isempty(parms.event)
  if ~isempty(opt.events)
    parms.event = opt.events;
  elseif ~isempty(opt.event_codes)
    parms.event = opt.event_codes;
  elseif ~isempty(opt.eventvals)
    parms.event = opt.eventvals;
  end
end

if ~isempty(opt.keepbadtrials_flag)
  parms.keepbadtrials = opt.keepbadtrials_flag;
end
if ~isempty(opt.keepbadchans_flag)
  parms.keepbadchans   = opt.keepbadchans_flag;
  parms.removebadchans = ~opt.keepbadchans_flag;
end

if isempty(parms.chanlabel) && ~isempty(opt.chanlabels)
  parms.chanlabel = opt.chanlabels;
end

if isempty(parms.badlabels) && ~isempty(opt.badlabel)
  parms.badlabels = opt.badlabel;
end

if isempty(parms.channel) && ~isempty(opt.channels)
  parms.channel = opt.channels;
end

if isempty(parms.rejectfile) && ~isempty(opt.reject_file)
  parms.rejectfile = opt.reject_file;
end

if isempty(parms.trial) && ~isempty(opt.trials)
  parms.trial = opt.trials;
end
end
