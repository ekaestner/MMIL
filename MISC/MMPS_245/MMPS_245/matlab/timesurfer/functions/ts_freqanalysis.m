function [timefreq_data,outfiles] = ts_freqanalysis(data,varargin)
%   cfg.output     = 'pow'       return the power-spectra
%                    'powandcsd' return the power and the cross-spectra

if nargin < 1, help(mfilename); end       % check number of inputs
error(nargoutchk(1,2,nargout,'string'));  % check number of outputs

% set up parms structure
parms = mmil_args2parms(varargin,...
                        {...
                         'complex_flag',1,[0 1],...
                         'method','wltconvol',[],...
                         'foi','all',[],...
                         'sf',[],[],...
                         'toi',[],[],...
                         'toilim',[],[],...
                         'events',[],[],...
                         'conditions',[],[],...
                         'channels',      [],   [],...
                         'chantype',      [], {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                         'channelcmb',    {'all' 'all'},[],...
                         'output','pow',[],...
                         'verbose',1,{0,1},...
                         'logfile',[],[],...
                         'logfid',1,[],...
												 'trials_flag', 0, sort([false true]),...
                         'save_flag', 1, sort([false true]),...
                         'overwrite',0,{0,1},...
                         'feedback','no',[],...
                         'removemean','no',[],...
                         'prefix','proc_tf',[],...
                         'rootoutdir',pwd,[],...
                         'blc','no',[],...
                         'blcwindow',[],[],...
                       },...
                       false);
  
% Backwards compatibility
parms = backcompatible(parms,varargin{:});     
if iscell(parms.events), parms.events = unique([parms.events{:}]); end
parms.keeptrials = 'yes'; % must start as yes for timelockanalysis before TF

% Select data to process
data  = ts_checkdata_header(data,'precision','double','events',parms.events);
data  = ts_data_selection(data,'channels',parms.channels,'chantype',parms.chantype,'events',parms.events,'conditions',parms.conditions); 
ncond = length(data.epochs);
parms.channels = 1:data.num_sensors;

% Check input data and parameters
if ~isfield(data,'epochs')
  mmil_error(parms,'Input data must be epoch_data. Set num_trials to 1 for continuous data.');
end

% check frequencies of interest
foi = [2:12 14:2:24 25:5:55 70:10:200];   % default TimeSurfer values
sf  = [1 1 2*ones(1,9) 3*ones(1,6) 5*ones(1,7) 10*ones(1,14)];
if ischar(parms.foi)                      % use pre-defined frequencies
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
      || (~isempty(parms.chantype) && ismember('eeg',{data.sensor_info.typestring})) ...
      || (~isempty(parms.chantype) && ismember('ieeg',{data.sensor_info.typestring}))
    parms.foi = foi;
    parms.sf  = sf;
  else
    parms.foi = foi;
    parms.sf  = sf;
  end
elseif isempty(parms.sf) || ~all(ismember(parms.sf,sf)) || ~all(ismember(parms.foi,foi)) || length(parms.sf)~=length(parms.foi)    
  % warn that we are not using the MMIL standard spec
  mmil_logstr(parms,'Not using the TimeSurfer-recommended frequencies of interest and spectral resolution!');
end
% otherwise, analyze 1-60Hz
if isempty(parms.foi) || ischar(parms.foi)
  parms.foi = 1:60;
end
clear foi sf

% check events
parms.conditions = num2cell(1:ncond);
parms.events     = num2cell([data.epochs.event_code]);

% check times of interest
if ~isempty(parms.toilim)
  parms.toi = parms.toilim(1):1/data.sfreq:parms.toilim(2);
elseif isempty(parms.toi)
  % set toi = intersection of all conditions' time vectors
  for c = 1:ncond
    if c == 1
      t0 = data.epochs(c).time(1);
      tf = data.epochs(c).time(end);
    else
      t0 = max(t0,data.epochs(c).time(1));
      tf = min(tf,data.epochs(c).time(end));
    end
  end
  parms.toi = t0:1/data.sfreq:tf;
end
clear t0 tf c

% Enforce Nyquist criterion
if data.sfreq < 2*max(parms.foi)
  mmil_error(parms,'The sampling frequency (%d) of the data is not high enough to compute frequencies as high as %d Hz.',data.sfreq,max(parms.foi));
end

% Create cfg structure with all input parameters
val      = struct2cell(parms);
fld      = fieldnames(parms);
addfld   = setdiff(varargin(1:2:end),fld);
if ~isempty(addfld)
  idx      = 2 * find(ismember(varargin(1:2:end),addfld)) -  1;
  allparms = cell2struct({varargin{idx+1},val{:}},{varargin{idx},fld{:}},2);
else
  allparms = parms;
end
clear val fld addfld idx

% define important flags: trials_flag, save_flag, complex_flag
trials_flag = parms.trials_flag || all([data.epochs.num_trials]==1);
save_flag   = parms.save_flag;

% store initial sensor info
orig_sensor_info = getfield(data,'sensor_info');

%% Time-Frequency Analysis
% Loop over conditions
warning off
outfiles = {};
for c = 1:ncond
  tic
  mmil_logstr(parms,'Converting to fieldtrip.');
  mmil_logstr(parms, 'Running t/f analysis, method=%s, on %d frequencies and %d times.',...
        parms.method, length(parms.foi), length(parms.toi));
  allparms.keeptrials = 'yes';
  allparms.channel    = 1:data.num_sensors;
	% convert to double precisiona
  if 0 % ~isa(data.epochs(c).data,'double')
    data.epochs(c).data = double(data.epochs(c).data);
  end
  
  % convert epoch_data to FieldTrip format
  FT_data = ts_data2fieldtrip(data,'condition',parms.conditions{c},'dimord','chan_time');
  
  % preprocessing with FieldTrip
  FT_data = timelockanalysis(allparms,FT_data);
  
	% keep track of indices
  [in_ts,in_ft] = match_str({data.sensor_info.label},FT_data.label);
	nchan = length(in_ts);
  
  % initialize some variables
  epoch_template             = rmfield(data,'epochs');
  epoch_template.epochs      = rmfield(data.epochs(c),'data');
  epoch_template.epochs.data = [];
  epoch_template.sensor_info = [];	        
  
  % clear input data for this condition
  data.epochs(c).data = [];
  data.epochs(c).time = [];
  
  % call FieldTrip function for TF analysis
  if trials_flag && ~save_flag  % proc trials & return structure
    allparms.channel           = 'all';
    allparms.keeptrials        = 'yes';
    FT_TFdata                  = freqanalysis(allparms,FT_data,1); clear FT_data
    epoch_template.sensor_info = data.sensor_info(in_ts);
    epoch_template.num_sensors = length(epoch_template.sensor_info);
    timefreq_data(c)           = ts_fieldtrip2timefreq(FT_TFdata,epoch_template);
    clear FT_TFdata
  elseif trials_flag            % proc trials & save files; return average
    allparms.keeptrials = 'yes';
    for ch = 1:nchan
      allparms.channel  = ch;
      FT_TFdata         = freqanalysis(allparms,FT_data,1);
      epoch_template.sensor_info = data.sensor_info(in_ts(ch));
      epoch_template.num_sensors = 1;
      tfdata                     = ts_fieldtrip2timefreq(FT_TFdata,epoch_template);
      clear FT_TFdata
      outfiles{end+1}            = save_trials(tfdata,allparms);      
      if ch==1
        % initialize average data structure for TF power
        timefreq_data(c)                = tfdata;
        timefreq_data(c).sensor_info    = data.sensor_info(in_ts);
        timefreq_data(c).num_sensors    = length(timefreq_data(c).sensor_info);
        timefreq_data(c).timefreq.power = zeros(nchan,length(tfdata.timefreq.time),length(parms.foi),'single');
        timefreq_data(c).timefreq.cmplx = single(complex(zeros(nchan,length(tfdata.timefreq.time),length(parms.foi)),...
                                                         zeros(nchan,length(tfdata.timefreq.time),length(parms.foi))));
      end
      timefreq_data(c).timefreq.power(ch,:,:) = single(mean(tfdata.timefreq.power,4));
      timefreq_data(c).timefreq.cmplx(ch,:,:) = single(mean(tfdata.timefreq.cmplx,4));
%         if isfield(timefreq_data(c).timefreq,'cmplx')
%           timefreq_data(c).timefreq     = rmfield(timefreq_data(c).timefreq,'cmplx');
%         end
      clear tfdata
    end
  else                          % proc average & return structure
    allparms.channel           = 'all';
    allparms.keeptrials        = 'no';
    FT_TFdata                  = freqanalysis(allparms,FT_data,0); clear FT_data
    epoch_template.sensor_info = data.sensor_info(in_ts);
    epoch_template.num_sensors = length(epoch_template.sensor_info);
    timefreq_data(c)           = ts_fieldtrip2timefreq(FT_TFdata,epoch_template);
    clear FT_TFdata
  end
  % copy meta-info for this condition
  timefreq_data(c).timefreq.event_code   = data.epochs(c).event_code;
  timefreq_data(c).timefreq.num_rejects  = data.epochs(c).num_rejects;
  timefreq_data(c).timefreq.num_trials   = data.epochs(c).num_trials;
  if isfield(data.epochs,'trial_info')
    timefreq_data(c).timefreq.trial_info = data.epochs(c).trial_info;
  end
%   % remove complex spectra if not specified
%   if ~parms.complex_flag && isfield(timefreq_data(c).timefreq,'cmplx')
%     timefreq_data(c).timefreq = rmfield(timefreq_data(c).timefreq,'cmplx');
%   end
  mmil_logstr(parms,'Elapsed time %.0f seconds.',toc);
end
% combine results
if length(timefreq_data) > 1
  timefreq_data = ts_combine_data(timefreq_data);
end
% store parameters
if isfield(data,'parms')
  allparms.procparms = data.parms;
end
allparms.filename   = outfiles;
timefreq_data.parms = allparms;

% restore all warning messages
warning on all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outfile = save_trials(timefreq_data,parms)
% save_timefreq_data => saves data split into files based on save_format

% remove complex spectra if not specified
if ~parms.complex_flag && isfield(timefreq_data.timefreq,'cmplx')
  timefreq_data.timefreq = rmfield(timefreq_data.timefreq,'cmplx');
else % keep only the complex spectra
  timefreq_data.timefreq = rmfield(timefreq_data.timefreq,'power');
end

% convert to single precision
if isfield(timefreq_data.timefreq,'power'), timefreq_data.timefreq.power = single(timefreq_data.timefreq.power); end
if isfield(timefreq_data.timefreq,'cmplx'), timefreq_data.timefreq.cmplx = single(timefreq_data.timefreq.cmplx); end
if isfield(timefreq_data.timefreq,'data') , timefreq_data.timefreq.data  = single(timefreq_data.timefreq.data);  end
if isfield(timefreq_data.timefreq,'cross'), timefreq_data.timefreq.cross = single(timefreq_data.timefreq.cross); end

str = parms.prefix;
if isequal(parms.blc,'yes'), str = [str '_blc']; end
str = sprintf('%s_timefreq_data_%d-%dHz',str,parms.foi(1),parms.foi(end));

% construct output file name: outfile
outdir  = [parms.rootoutdir '/matfiles'];
if ~exist(outdir,'dir'), mkdir(outdir); end
outfile = sprintf('%s/%s_event%02i_chan%03i.mat',outdir,str,timefreq_data.timefreq.event_code,parms.channel);

if exist(outfile,'file') && ~parms.overwrite
	mmil_logstr(parms,'Not overwriting trial data file: %s',outfile);
  error('Remove files or set "overwrite" to 1.'); % consider treating this differently...
else
  mmil_logstr(parms,'Saving trial data for channel %d of %d: %s',parms.channel,length(parms.channels),outfile);
  timefreq_data.parms = parms;
  timefreq_data.parms.filename = {outfile};
	save(outfile,'timefreq_data','-v7.3');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
        'event_codes',[],[],...
        'savecomplex_flag',[],[],...
        'blc',[],[],...
        'overwrite',[],[],...
        },false);
  
if isempty(parms.events) && ~isempty(opt.event_codes)
  parms.events = opt.event_codes;
end
if ~isempty(opt.savecomplex_flag)
  parms.complex_flag = opt.savecomplex_flag;
end
if ~isempty(opt.blc)
  if isequal(opt.blc,false)
    parms.blc = 'no';
  elseif isequal(opt.blc,true)
    parms.blc = 'yes';
  end
end
if ischar(opt.overwrite)
  if strcmpi(opt.overwrite,'yes')
    parms.overwrite = 1;
  elseif strcmpi(opt.overwrite,'no')
    parms.overwrite = 0;
  end
end