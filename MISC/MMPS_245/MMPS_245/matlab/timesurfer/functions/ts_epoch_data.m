function data = ts_epoch_data(data,varargin)
% Purpose: input continuous or epoched data in TimeSurfer format along with
% event info and return epoched (or re-epoched) data.
%
% Usage: data = ts_epoch_data(data,'key1',val1, ... )
%
% Example: epoch_data = ts_epoch_data(cont_data,...
%                             'trigchan'         ,'STI 101'...
%                             'write_events_file',1,...
%                             'rootoutdir'       ,'/home/user/rootoutdir',...
%                             'prefix'           ,'prefix',...
%                             'prestim'          ,.3,...
%                             'poststim'         ,.8);
%
% Example: epoch_data = ts_epoch_data(epoch_data,...
%                             'evntfile','/home/user/events.txt',...
%                             'events'  ,[16 32]);
%
% Inputs: 
%     data: timesurfer data struct of continuous or epoch data 
%     evntfile: full filepath to ev2, evt (BESA), nev, or .txt file
%       -OR-         containing events information to be used for epoching  
%     trigchan: either a amplitude-coded channel label eg. 'STI 101' or an array of
%       -OR-         binary trigger channel labels eg. {'STI001' 'STI002' 'STI 003'}    
%     samples: an array of sample points to be used for epoching
%
% Outputs:
%     data: timesurfer data struct of epoch data 
%
% Parameters: default : options :  description 
%     trigthresh: [] :: trigger threshold, by default: the mean plus one
%     standard deviation of trigger amplitudes
%     trig_minduration: 5 ::
%     conds: [] ::
%     events: [] :: vector of event codes of conditions to process
%     old_event_codes: [] :: vector of old event codes
%     new_event_codes: [] :: vector of new event codes
%                            Note: old_event_codes and new_event_codes
%                            parameters are used in conjection to change
%                            event codes.
%     datafile: '' :: full path of raw datafile to be included in trial_info
%                     section of created data structure
%     prestim: .2 :: duration of prestimulus period (sec)
%     poststim: .6 :: duration of poststimulus period (sec)
%     channels: [] :: vector of channel indices - overrides chantype
%     write_events_file: 0 : 0,1 : whether to generate a .txt file containing events information
%     prefix: 'proc' :: string to be prepended to any generated files
%     rootoutdir: [] :: path to directory where any generated files will be
%                       saved, do not include trailing '/'
%     calc_ncov_flag: 1 : 0,1 : whether to calculate/recalculate noise
%                               covariance matrix after epoching/re-epoching
%     noisewindow: [] :: noise covariance matrix calculation window, [begin end] in seconds,
%                        the default is the entire epoch ([-inf inf])
%     ncov_ex_evnts: [] :: vector of event codes that should not be used in calculating the noise covariance matrix
%                          e.g. to exclude events with pre-stimulus activity that is not noise -e.g. a button press
%     concat_ncov_flag: 0 : 0,1 : whether to concatenate noise before calculating noise covarience (1) or 
%                                 calculate noise covarience for each trial and then average (0)
%     max_noise_samples: 5000 :: max # of samples used for noise estimate
%                                Note: especially significant if concat_ncov_flag = 1
%     max_num_trials: Inf :: maximum number of trials per condition
%     zparam : [] :: name of the data field to epoch (string; default
%                    depends on the data type)
%     stim_delay: 0 :: time delay to add to offset the data accounting for delay from trigger to actual stimulus.  SRD 032112
%     pulsetrain: 0 : 0,1 : whether the tigger channel contains pulsetrains (1) or not (0)
%     pt_trigthresh: 13000 : [] : trigger threshold 
%     pt_digthresh: 16000 : [] : digital marker theshold 
%     npulses: 4 : [] : # of marker pulses (not including trigger pulse)
%
% Methods for specifying event information.
% evntfile:
%   txt (events_fnames), ev2, evt (BESA), nev
% trigchan:
%   sensor label (amplitude code) or labels (binary code)
% vector:
%   absolute sample points from start of file (latencies) + corresponding
%   condition codes (conds)
% evnts structure (not yet implemented)
%
% Re-epoch = select subset of trials from epoch_data based on info in an event file
%
% trig_minduration (in samples)
%
% Created:  05-Jan-2011, Jason Sherfey
% Modified: 23-Feb-2011, Jason Sherfey
% Modified: 21-Mar-2012, Steve Deiss (fixing stim delay)
% Modified: 17-May-2012, Steve Deiss (fix double underscore in event file name written.)
% Modified: 15-Mar-2013, Burke Rosen (fix rounding errors and add .nev eventsfile type
% Modified: 15-Mar-2013, Burke Rosen added support for pulsetrain triggers

%  TO ADD CHECK FOR EPOCH DATA W/O TRIAL_INFO IN TS_EPOCH_DATA & ADD IT

% set default trig_minduration to 5?   
parms = mmil_args2parms( varargin, ...
    { 'datafile','',[],...
    'evntfile',[],[],...
    'exp_file',[],[],...
    'samples',[],[],...
    'trigchan',[],[],...
    'trigthresh',[],[],...
    'trig_minduration',5,[],...
    'trig_minseparation',[],[],...
    'events',[],[],...
    'old_event_codes',[],[],...
    'new_event_codes',[],[],...
    'conds',[],[],...
    'prestim' ,.2,[],...
    'poststim',.6,[],...
    'channels',[],[],...
    'write_events_file',0,[],...
    'rootoutdir',[],[],...
    'prefix','proc',[],...
    'calc_ncov_flag',true,[false true],...
    'noisewindow',[-inf inf],[],...
    'ncov_ex_evnts',[],[],...
    'max_noise_samples',5000,[],...
    'concat_ncov_flag',false,[false true],...
    'max_num_trials',Inf,[1 Inf],...
    'zparam',[],[],...
    'stim_delay',0,[],...   %%% SRD, 032012: Added stim_delay parameter with 0 as default
    'pulsetrain',0,[],...
    'pt_trigthresh',13000,[],...
    'pt_digthresh',16000,[],...
    'npulses',4,[],...
    }, ...
    false ); 

if ischar(data)
  parms.datafile = data;
  data = ts_load_data(data);
elseif iscellstr(data)
  error('The input must be a data structure or filename. Use ts_process_data to process multiple files at once.');
elseif ~isstruct(data)
  error('The input must be TimeSurfer data or a single filename.');
end

[datatype,datafield,dataparam] = ts_object_info(data);
if isempty(parms.zparam) || ~ismember(parms.zparam,dataparam)
  parms.zparam = dataparam{1}; % name of data field to epoch
end

parms = backcompatible(parms,varargin{:});   
nsamp = length(data.(datafield)(1).time);
Fs    = data.sfreq;

% assume data are continuous if num_trials = 1 & only one condition in data
continuous_flag = all([data.(datafield).num_trials]==1) && length(data.(datafield))==1;
  % whether the INPUT is continuous or epoched
trigchan_flag   = ~isempty(parms.trigchan) && continuous_flag && any(ismember(parms.trigchan,{data.sensor_info.label}));
  % whether a trigger channel should be used to epoch continuous data
evntfile_flag   = ~isempty(parms.evntfile);% && exist(parms.evntfile,'file');
  % whether an event file should be used to epoch continuous or epoched data
samples_flag    = (isnumeric(parms.samples) && ~isempty(parms.samples));
  % whether an array of sample points should be used for epoching
if ~(trigchan_flag || evntfile_flag || samples_flag)
  fprintf('Improper event specification. Must supply existing event file, trigchan with continuous data, or array of samples for epoching.\n');
  return;
end

if continuous_flag
  parms.trialdim = ndims(data.(datafield)(1).(parms.zparam)) + 1; % trial dimension
else
  parms.trialdim = ndims(data.(datafield)(1).(parms.zparam));
end

%% get event info
%   => samp: absolute sample number relative to start of datafile
%   => cond: condition corresponding to each sample
samp = [];
cond = [];

if trigchan_flag
  % if an evntfile or samples are also supplied, epoch with trigchan 1st
  if  evntfile_flag || samples_flag
    % *RECURSIVE* call to ts_epoch_data()
    % - to epoch data based on all triggers before applying constraints in
    % evntfile or samples
    % - ensures that absolute trial numbers and base conditions in trial_info
    % are determined by the experimental setup & not the evntfile or otherwise.
    tmpcfg = parms;
    tmpcfg.samples  = [];
    tmpcfg.evntfile = [];
    args = mmil_parms2args(tmpcfg);
    data = ts_epoch_data(data,args{:}); % cont_data => epoch_data
    continuous_flag = 0;
    % NOTE: cond & samp are still empty at this point
  else  % epoch continuous data using trigchan(s)
    % extract event info from trigger channels
    if ischar(parms.trigchan)
      % single trigger channel (amplitude, decimal code)
      trigdata    = ts_data_selection(data,'chanlabel',parms.trigchan,'verbose',0);
      if ~parms.pulsetrain
        [samp,cond] = extract_amplitude_code(trigdata,parms);
      else
        [samp,cond] = extract_pulsetrain_code(trigdata,parms);
      end
    elseif iscell(parms.trigchan)
      % multiple trigger channels (distributed, binary code)
      % pass triggers to subfunction to extract binary code & convert to decimal
      trigdata    = ts_data_selection(data,'chanlabel',parms.trigchan,'verbose',0);
      [samp,cond] = extract_binary_code(trigdata,parms);
      clear trigdata
    end
    name = repmat({'unknown'},1,length(samp));  % TODO: can this info be obtained?
    dura = ones(1,length(samp));              % TODO: can this info be obtained?
    % NOTE: cond & samp are no longer empty at this point
  end
end

% get event info from evntfile or sample_vec if not specified
if isempty(samp)
  if evntfile_flag
    % Event info specified in an ascii text file (.txt, .ev2, .evt)
    if ~exist(parms.evntfile,'file')
      fprintf('Event file does not exist: %s\n',parms.evntfile);
      return;
    end
    [path,name,ext] = fileparts(parms.evntfile);
    if length(name)>4 && strcmp(name(end-3:end),'.dio')
      ext = ['.dio' ext];
    elseif ~strcmp(ext,'.nev')
      [content,res]   = mmil_readtext(parms.evntfile,'\t');
    end
    switch lower(ext)
      case '.dio.txt'
        % *.dio file used to epoch NSpike data
        if ~exist(parms.exp_file,'file')
          error('An experiment file must be given to epoch NSpike data.');
        end
        load(parms.exp_file,'experiment');
        [time_codes, A, B, C, D] = ntools_dio_parse(parms.evntfile);
        NSPIKE_SAMPLING_RATE 				= 30000;  %  NSpike timestamps are in Nspike system sampling rate
        sfreq_to_nspike_sfreq_ratio = data.sfreq / NSPIKE_SAMPLING_RATE;
        cond = nan(1,length(C));
        samp = nan(1,length(C));
        name = repmat({''},1,length(C));
        dura = nan(1,length(C));
        for i = 1:length(experiment.event)
          experiment.event(i).timestamps = time_codes(ismember(C, experiment.event(i).source_code));
          code = experiment.event(i).source_code;
          ind  = ismember(C,code);
          cond(ind) = code;
          samp(ind) = round(experiment.event(i).timestamps * sfreq_to_nspike_sfreq_ratio);
          name(ind) = {experiment.event(i).name};
          dura(ind) = 1;
        end
%         if isfield(data,'experiment') && ~isempty(data.experiment) && length(fieldnames(data.experiment))>2
%           data.experiment.dio_file = parms.evntfile;
%         else
          data.experiment = experiment;
          data.experiment.dat_file = parms.datafile;
          data.experiment.dio_file = parms.evntfile;
%         end
      case '.txt'
        % events_fnames (TimeSurfer MEG stream at MMIL)
        evnt = content(2:end,:);
        name = evnt(:,1)';    % type: Trigger, Response
        samp = [evnt{:,2}]';
        cond = [evnt{:,3}]';
        dura = [evnt{:,4}]';
        rejects_ix = find(cellfun(@(x)isequal(x,'Reject'),name));
        if ~isempty(rejects_ix)
          name(rejects_ix) = [];
          samp(rejects_ix) = [];
          cond(rejects_ix) = [];
          dura(rejects_ix) = [];
        end
      case '.ev2'
        % *.ev2 file (NYU, MGH)
        cond = [content{:,2}]';
        samp = [content{:,6}]';
        name = repmat({'unknown'},1,length(samp));
%         name = content(:,3)';           % check this
        dura = ones(1,length(samp));    % TODO: can this info be obtained?
      case '.evt'
        % BESA event file
        Tunits = content{1,1};
        if strcmp(Tunits,'Tms')
          c  = 1E3;
        elseif strcmp(Tunits,'Tmu')
          c  = 1E6;
        else
          c  = 1;
        end
        evnt = content(2:end,:);
        samp = round(Fs*([evnt{:,1}]')/c);
        cond = [evnt{:,2}]';
        if size(evnt,2) > 3
          name = evnt(:,4)';            % check this
        else
          name = repmat({'unknown'},1,length(samp));
        end
        dura = ones(1,length(samp));    % TODO: can this info be obtained?
        clear Tunits c
      case '.edfevt'
        % not supported at this time (old format used in ts_loadedf)
        fprintf('Legacy evntfile format created for ts_loadedf is no longer supported.\n');
        return;
        case '.nev' %BQR 13.03.15
            % Blackrock binary trigger file
            NEV = openNEV(parms.evntfile,'nosave','nomat','8bits');
            BLKRK_SAMPLING_RATE = double(NEV.MetaTags.SampleRes);  %  Blackrock timestamps are in Blackrock system sampling rate (30kHz for ns5)
            sfreq_to_blkrk_sfreq_ratio = data.sfreq / BLKRK_SAMPLING_RATE;
            cond = NEV.Data.SerialDigitalIO.UnparsedData;
            samp = round(NEV.Data.SerialDigitalIO.TimeStamp .* sfreq_to_blkrk_sfreq_ratio);
            dura = ones(1,length(samp));
            name = repmat({'unknown'},1,length(samp));
    end
  elseif samples_flag
    % array of sample points
    samp = parms.samples;
    if isempty(parms.conds), parms.conds = data.(datafield)(1).event_code; end
    if length(parms.conds) == 1
      parms.conds = parms.conds * ones(1,length(samp));
    elseif length(parms.conds) ~= length(samp)
      fprintf('The number of conditions given must equal the number of triggers.\n');
      return;
    end
    cond = parms.conds;
    name = repmat({'unknown'},1,length(samp));
    dura = ones(1,length(samp));
  else
    fprintf('Failed to find event information for epoching.\n');
    return;
  end
end

% keep valid event codes
if ~isempty(parms.events)
  % note: new event codes could be specified in the evntfile that do not
  % correspond to the original base conditions. parms.events indicates which
  % event codes of those specified in the evntfile should be kept (not the
  % base conditions). Original base conditions are recorded in trial_info.
  name = name(ismember(cond,parms.events));
  dura = dura(ismember(cond,parms.events));
  samp = samp(ismember(cond,parms.events));
  cond = cond(ismember(cond,parms.events));
end

% rename event codes
if ~isempty(parms.old_event_codes) && length(parms.old_event_codes)==length(parms.new_event_codes)
  for k  = 1:length(parms.old_event_codes)
    cond(ismember(cond,parms.old_event_codes(k))) = parms.new_event_codes(k);
  end
end

if size(samp,1) == size(cond,2)
  samp = samp';
end

% prepare event info
if continuous_flag
  keep        = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
  keep        = keep & ~isnan(cond);
  samp        = samp(keep);
  cond        = cond(keep);
  name        = name(keep);
  dura        = dura(keep);
end
[conds,ii]  = unique(cond);
nevnt       = length(samp);
ncond       = length(conds);
samp        = samp + (parms.stim_delay * Fs) ; %%% SRD 032112: Making stim_delay work for data epoched without triggers (e.g. Nspike event file) 
begsample   = round(samp - parms.prestim * Fs ) + 1;
endsample   = round(samp + parms.poststim * Fs);

if nevnt == 0
  % the data does not contain any of the desired samples
  fprintf('Error: No triggers were found in the given data.\n');
  return;
end

% sort by latency
[jnk,idx] = sort(samp);
samp      = samp(idx); 
cond      = cond(idx);
name      = name(idx);
dura      = dura(idx);

% % switched order of begsample and endsample def from above, Andrew Il Yang,
% % 08/08/11
% begsample   = round(samp - parms.prestim*Fs ) + 1;
% endsample   = round(samp + parms.poststim*Fs);

%% Epoch continuous data or re-epoch epoched data
if continuous_flag
  % continuous data
  contdata         = data.(datafield).(parms.zparam);
  data.(datafield).time = [];
  data.(datafield).(parms.zparam) = [];
  data.(datafield)(1:ncond) = data.(datafield);
  if isempty(parms.channels)
    chans = 1:size(contdata,1);
  else
    chans = parms.channels;
  end
%   T     = [-parms.prestim:1/Fs:parms.poststim];
  T = linspace(-parms.prestim,parms.poststim,endsample(1)-begsample(1)+1);
  for c = 1:ncond
    ii  = cond==conds(c);
    s0  = num2cell(begsample(ii));
    sf  = num2cell(endsample(ii));
    % fix rounding errors, added bqr 13.02.22
    while numel(unique(cellfun(@(x,y) x-y,sf,s0)))>1
        [~,maxidx] = max(cellfun(@(x,y) x-y,sf,s0));
        sf{maxidx} = sf{maxidx}-1;
    end
    %
    tmp = cellfun(@(x,y)(contdata(chans,x:y,:)),s0,sf,'UniformOutput',false);
    tmpname = unique(name(ii));
    if iscell(tmpname) && length(tmpname)==1, tmpname=tmpname{1}; end
    data.(datafield)(c).name       = tmpname;
    data.(datafield)(c).event_code = conds(c);
    data.(datafield)(c).num_trials = length(s0);
    data.(datafield)(c).trial_info = [];
    data.(datafield)(c).trial_info = add_trial_info;
    data.(datafield)(c).time       = T;%(1:size(tmp{1},2))/Fs - parms.prestim;%T;
    data.(datafield)(c).(parms.zparam)       = cat(parms.trialdim,tmp{:}); % channel x time x trial
    if iscell(data.(datafield)(c).name) && length(data.(datafield)(c).name)==1
      data.(datafield)(c).name     = data.(datafield)(c).name{1};
    end
    clear tmp
  end
  clear contdata
else
  % epoched data (re-epoch)
  re_epoch;
end

% select subset of trials if maximum number was specified
if parms.max_num_trials < max([data.(datafield).num_trials])
  for k = 1:length(data.(datafield))
    num_trials{k} = 1:min(data.(datafield)(k).num_trials,parms.max_num_trials);
  end
  data = ts_data_selection(data,'trials',num_trials);
end

% calculate noise covariance matrix
if parms.calc_ncov_flag && ~strcmp(datafield,'timefreq')
  calc_ncov;
else
  data.noise.num_trials  = 0;
  data.noise.num_samples = 0;
  data.noise.covar       = [];    
end

% Write log
if parms.write_events_file && exist('samp','var')
  write_log;
end

  %% NESTED FUNCTIONS

  function trial_info = add_trial_info
    % Build the trial info field.
    trial_info=[];
    trial_info.number(1,:)                           = find(ii);
    trial_info.latency(1,:)                          = samp(ii);
    trial_info.badtrial(1,:)                         = zeros(1,length(s0));
    trial_info.event_code(1,:)                       = cond(ii);
    trial_info.duration(1,:)                         = dura(ii);
    [trial_info.datafile{1,1:length(s0)}]       = deal(parms.datafile);
    [trial_info.events_fnames{1,1:length(s0)}]  = deal(parms.evntfile);
  end

  function re_epoch
    % shuffle epochs in epoch_data based on event info
    epoch_data  = data;
    clear data
    data        = rmfield(epoch_data,datafield);
    tmp         = rmfield(epoch_data.(datafield)(1),dataparam);%parms.zparam);
    [data.(datafield)(1:ncond)] = deal(tmp);
    clear tmp
    latency = arrayfun(@(x)x.trial_info.latency,epoch_data.(datafield),'uniformoutput',false);
    N       = cellfun(@length,latency);                         % if N   = [40 50 60];
    cum     = cellfun(@(x)sum(N(1:x)),num2cell(1:length(N)));   %    cum = [40 90 150]
    latvec  = [latency{:}];
    
    if all(cellfun(@(x)ismember(x-1,latvec),num2cell(samp))), samp = samp-1; end
    if all(cellfun(@(x)ismember(x+1,latvec),num2cell(samp))), samp = samp+1; end
    
    keep = cellfun(@(x)ismember(x,latvec),num2cell(samp));
    samp = samp(keep);
    cond=cond(keep);
    for k = 1:ncond
      s = cond == conds(k);
      tmpsamp = samp(s);  
      allind  = cellfun(@(x)find(x==latvec),num2cell(tmpsamp));
      c       = cellfun(@(x)find(cum >= x,1,'first'),num2cell(allind));                     % c = index into epochs for each sample of this cond
      i       = cellfun(@(x,y)find(ismember(latency{x},y)),num2cell(c),num2cell(tmpsamp));  % i = index into latency in epochs(c).trial_info
      data.(datafield)(k).trial_info.number        = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.number(y)       ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.latency       = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.latency(y)      ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.badtrial      = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.badtrial(y)     ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.event_code    = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.event_code(y)   ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.duration      = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.duration(y)     ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.datafile      = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.datafile(y)     ,num2cell(c),num2cell(i))';
      data.(datafield)(k).trial_info.events_fnames = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.events_fnames(y),num2cell(c),num2cell(i))';
      if evntfile_flag
        [data.(datafield)(k).trial_info.events_fnames{1:length(c)}] = deal(parms.evntfile);
      end
      if isfield(data.(datafield)(k).trial_info,'epoch_num')
        data.(datafield)(k).trial_info.epoch_num   = cellfun(@(x,y)epoch_data.(datafield)(x).trial_info.epoch_num(y)   ,num2cell(c),num2cell(i))';
      end
      if parms.trialdim==3
        dat                                        = cellfun(@(x,y)epoch_data.(datafield)(x).(parms.zparam)(:,:,y),num2cell(c),num2cell(i),'uniformoutput',false);
      elseif parms.trialdim==4
        dat                                        = cellfun(@(x,y)epoch_data.(datafield)(x).(parms.zparam)(:,:,:,y),num2cell(c),num2cell(i),'uniformoutput',false);  
      else
        error('Trials must be dimension 3 or 4.');
      end
      data.(datafield)(k).(parms.zparam)           = cat(parms.trialdim,dat{:});  
      data.(datafield)(k).num_rejects              = cell2struct({0,0,0,0,0,0},{'mag','grad','eeg','eog','manual','skip'},2);
      clear dat
    end
    ntrl  = num2cell(arrayfun(@(x)size(x.(parms.zparam),parms.trialdim),data.(datafield)));
    conds = num2cell(conds);
    [data.(datafield).num_trials] = deal(ntrl{:});
    [data.(datafield).event_code] = deal(conds{:});
    clear epoch_data
  end

  function write_log
    % where to save?
    if ~isempty(parms.rootoutdir)
      outdir = parms.rootoutdir;
    elseif ~isempty(parms.evntfile)
      outdir = fileparts(parms.evntfile);
    else
      outdir = pwd;
    end
    name = sprintf('%s',parms.prefix);
    if ischar(parms.datafile)   
      if ~strcmp(name(end),'_'); 
          name = [name '_']; %%%%%  SRD 12.05.17 Without this if qualif, we get two '_' in the evnt file name. (Annoying) -updated 13.03.15 BQR
      end
    end
    name = [name 'events.txt'];
    logfile = fullfile(outdir,name);
    % write the log
    fid     = fopen(logfile,'wt');
%     [~,b]   = unix('echo $HOST');
%     fprintf(fid,'--------------------------------------------------------------------\n');
%     fprintf(fid,'Time: %s\n',datestr(now));
%     fprintf(fid,'Host: %s',b);
%     fprintf(fid,'Epoching %s data:\n%s\n',ext,parms.datafile);
%     fprintf(fid,'Number channels: %i \n',data.num_sensors);
%     fprintf(fid,'Sampling frequency: %5.6g Hz\n',Fs);
%     fprintf(fid,'--------------------------------------------------------------------\n');
    fprintf(fid, '%s\t%s\t%s\t%s\n','type','latency','condition','duration'); %% bqr 12.03.13
    for k = 1:nevnt
        fprintf(fid,'%s\t%i\t%i\t%i\n','Trigger',samp(k),cond(k),dura(k));
    end
%     fprintf(fid,'--------------------------------------------------------------------\n');
    fclose(fid);
  end

  function calc_ncov
    % add noise estimate
    noise_events          = setdiff([data.(datafield).event_code],parms.ncov_ex_evnts);
    noise                 = ts_data_selection(data,'events',noise_events,'toilim',parms.noisewindow);
    noise                 = {noise.(datafield).(parms.zparam)};
    data.noise            = [];
    data.noise.num_trials = sum([data.(datafield).num_trials]);
    if parms.concat_ncov_flag   % concatenate noise from all trials then calc cov
      noise = cellfun(@(x)reshape(x,[size(x,1) size(x,2)*size(x,3)]),noise,'uniformoutput',false);
      noise = cat(2,noise{:});  % chan x samp, contains noise samps from all trials and conditions
      if size(noise,2) > parms.max_noise_samples
        noise = noise(:,1:parms.max_noise_samples);
      end
      data.noise.num_samples = size(noise,2);
      data.noise.covar       = cov(noise');
    else                        % average over single trial noise cov matrices
      noise = cat(3,noise{:});
      if size(noise,2) > parms.max_noise_samples
        noise = noise(:,1:parms.max_noise_samples,:);
      end
      noise                  = squeeze(mat2cell(noise,size(noise,1),size(noise,2),ones(1,size(noise,3))));
      data.noise.num_samples = sum(cellfun(@(x)size(x,2),noise));
      data.noise.covar       = cellfun(@(x)cov(x'),noise,'uniformoutput',false);
      data.noise.covar       = cat(3,data.noise.covar{:});
      data.noise.covar       = mean(data.noise.covar,3);
    end       
  end

end

%% SUBFUNCTIONS
function [samp,cond] = extract_pulsetrain_code(trigdata,parms)
    [~,datafield] = ts_object_info(trigdata);
    data = trigdata.(datafield).(parms.zparam);
    ind = ts_crossing(data,[],parms.pt_trigthresh);
    ind = ind(1:2:end);% leading edge
    ind = ind+2;%  plateau offset
    rnddata = round(data(ind)/1000)*1000;
    digkey = zeros(parms.npulses,1);
    digtrig = zeros(parms.npulses,size(rnddata,2)/(parms.npulses+1));
    for ipul = 1:parms.npulses
        digkey(ipul) = 2^(ipul-1);
        digtrig(ipul,:) = rnddata(ipul+1:parms.npulses+1:end);
    end
    digkey = repmat(digkey,1,size(digtrig,2));
    digtrig(digtrig<parms.pt_digthresh) = 0;
    digtrig(digtrig>parms.pt_digthresh) = 1;
    digtrig = digtrig.*digkey;
    cond = sum(digtrig);
    ind  = ind-2;% correct for offset
    samp = ind(1:parms.npulses+1:end);
end

function [samp,cond] = extract_amplitude_code(trigdata,parms)
  % based on ts_read_fif_events()
  [datatype,datafield] = ts_object_info(trigdata);
  X   = trigdata.(datafield).(parms.zparam); % 1-D array
  dX  = diff(X);
  change_idx = find(dX~=0);
  change_idx = change_idx(change_idx+parms.trig_minduration<length(dX));
  change_len = diff(change_idx);
  onset_idx  = find(dX >0);
  onset_idx  = onset_idx(onset_idx+parms.trig_minduration<length(dX));
  noise_idx  = change_idx(change_len < parms.trig_minduration);
  noise_onset= [];
  for k = 1:length(noise_idx)
    trig_onset_idx = noise_idx(k);
    trig_code      = dX(trig_onset_idx);
    next_samples   = dX(trig_onset_idx+1:trig_onset_idx+parms.trig_minduration);
    trig_offset_idx= trig_onset_idx + find(next_samples == -trig_code) - 1;
    % noise found (remove from 'differences')
    if ~isempty(trig_offset_idx)
      dX(trig_onset_idx)  = 0;
      dX(trig_offset_idx) = 0;
      noise_onset(end+1)  = trig_onset_idx;
    end
  end
  onset_idx  = find(dX >0);   % refresh onset w/o noise
  trig_onset = onset_idx;% + 1; % adjust to account for the element lost by diff
  cond = dX(onset_idx);
  samp = trig_onset;      
  % try to correct for non-integer event code multiples
  if ~isequal(cond,round(cond))
    tmpcond = cond / min(cond);
    if isequal(tmpcond,round(tmpcond))
      cond = tmpcond;
    end
  end
end

function [samp,cond] = extract_binary_code(trigdata,parms)
  % assume these channels carry a binary code for epoching
  Fs = trigdata.sfreq;
  [datatype,datafield] = ts_object_info(trigdata);
  % set threshold for trigger detection
  if isempty(parms.trigthresh) || ~isnumeric(parms.trigthresh)
    tmp              = trigdata.(datafield).(parms.zparam)(:);
    parms.trigthresh = mean(tmp) + std(tmp);
    clear tmp
  end
  % find threshold crossings
  t     = trigdata.(datafield).time;
  nsamp = length(t);
  ntrigchan = trigdata.num_sensors;
  [crsind{1:ntrigchan}] = deal([]);
  for k = 1:ntrigchan
    tmpdat  = trigdata.(datafield).(parms.zparam)(k,:);
    ind     = ts_crossing(tmpdat,[],parms.trigthresh); %%%% changed call of crossing to ts_crossing, SRD 032612
    % only keep left edge of pulse
    sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
    ind     = ind(sel);
%     % only keep crossings that occur > minduration since last crossing, if specified
%     if ~isempty(parms.trig_minduration)
%       sel = [1 find((t(ind(2:end))-t(ind(1:end-1))) >= parms.trig_minduration/Fs)];
%       ind = ind(sel);
%     end
    crsind{k} = ind;
    clear ind sel tmpdat
  end
  bincode = zeros(ntrigchan,nsamp);
  for k   = 1:length(crsind), bincode(k,crsind{k}) = 1; end
  detind  = any(bincode,1); % whether there's a detection in any trigchan
  % make sure detections for the same pulse are grouped before conversion
    if isempty(parms.trig_minseparation)
      parms.trig_minseparation = .01;
    end
    delta   = round(parms.trig_minseparation*Fs/2);
    % get indices at which a detection occurs in any trigchan
    ind     = find(detind);
    % find detections with a second detection > delta steps before it
    tmp     = cellfun(@(x)any((ind>x-delta & (ind<x))),num2cell(ind));
    R_ind   = find(tmp);
    % find detections with a second detection < delta steps after it
    tmp     = cellfun(@(x)any((ind<x+delta & (ind>x))),num2cell(ind));
    L_ind   = find(tmp);
    % remove detections that are between two other detections in the same pulse
    Rsel    = ~ismember(R_ind,L_ind);
    Lsel    = ~ismember(L_ind,R_ind);
    R_ind   = R_ind(Rsel); clear Rsel
    L_ind   = L_ind(Lsel); clear Lsel
    % for each pair (L,R), set [ind(L):ind(R)] to 1 in each bincode row w/ a detection
    for k = 1:ntrigchan
      sel = cellfun(@(x,y)any(bincode(k,ind(x):ind(y))),num2cell(L_ind),num2cell(R_ind));
      sel = cellfun(@(x,y)ind(x):ind(y),num2cell(L_ind(sel)),num2cell(R_ind(sel)),'uniformoutput',false);
      sel = unique([sel{:}]);
      bincode(k,sel) = 1;
    end
  detind  = any(bincode,1); % whether there's a detection in any trigchan    
  % only keep the first pulse in contiguous detections
  samp    = find(detind(1:end-1)==0 & detind(2:end)==1) + 1;
  bincode = flipud(bincode(:,samp));
  bincode = mat2cell(bincode,ntrigchan,ones(1,length(samp)));
  % convert binary to decimal
  cond    = cellfun(@(x)polyval(x,2),bincode);
end

function parms = backcompatible(parms,varargin)
  opt = mmil_args2parms(varargin,{...
          'event',[],[],...
          'event_codes',[],[],...
          'valid_event_codes',[],[]...
          'eventvals',[],[],...
          'events_fnames',[],[],...
          'event_fnames',[],[],...
          'trigchans',[],[],...
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
end

