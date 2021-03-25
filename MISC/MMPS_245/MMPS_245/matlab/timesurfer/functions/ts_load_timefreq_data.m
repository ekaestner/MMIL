function [timefreq_data,files] = ts_load_timefreq_data(varargin)
% Purpose: load & combine TF trial spectra from several files.
%
% Different ways to specify timefreq_data files to load:
% Method 1: simple cell array of strings listing the files
% Method 2: pass structure with parms.filename listing the files
% Method 3: search directory using info from key/value pairs
%
% Assumes: each file contains one channel for one condition; all times,
% frequencies, and trials.
%
% Examples #1 & #2:
% foi           = [8 10 12];
% data          = ts_matrix2data(rand(10,100,5));
% tf_avg_data   = ts_freqanalysis(data,'foi',foi,'trials_flag',1,'save_flag',1);
% tf_epoch_data = ts_load_timefreq_data(tf_avg_data.parms.filename);
%   or
% tf_epoch_data = ts_load_timefreq_data(tf_avg_data);
%
% Example #3:
% datapath = '/space/mdkm1/2/kmdev/projects/NV/NY25/matfiles/timefreq';
% events   = [2011 2021];
% channels = 1:64; 
% data     = ts_load_timefreq_data('datapath',datapath,'events',events,'channels',channels);
%
% Created by JSS on 29-Apr-2011
% Modified last by JSS on 29-Apr-2011

% check inputs
data = [];
if mod(nargin,2)      % odd number of inputs
  data = varargin{1};
  if nargin > 1
    varargin = varargin(2:end);
  else
    varargin = {};
  end
end

% set up parms structure using key/value pairs
parms = mmil_args2parms(varargin,...
                        {...
                         'datapath',pwd,[],...
                         'events',[],[],...
                         'channels',[],[],...
                         'foi',[],[],...
                         'foilim',[],[],...
                         'toi',[],[],...
                         'toilim',[],[],...
                         'num_sensors',[],[],...
                         'prefix',[],[],...
                         'verbose',1,{0,1},...
                        },false);

% Get list of files
if iscellstr(data)            % Method 1 (list of files)
  files = data;
elseif isstruct(data)         % Method 2 (structure w/ parms.filename)
  if issubfield(data,'parms.filename')
    files = data.parms.filename;
    if isfield(data,'num_sensors')
      parms.num_sensors = data.num_sensors;
    end
  else
    error('Structure must contain list of files in parms.filename');
  end
else                          % Method 3 (search directory)
  contents = dir(parms.datapath);
  files    = {contents.name};
  files    = files(~[contents.isdir]);
end
clear data contents

% check for existence & create full filenames if necessary
rej = [];
for k = 1:length(files)
  if exist(fullfile(parms.datapath,files{k}),'file')
    % add path to file name
    files{k} = fullfile(parms.datapath,files{k});
  elseif ~exist(files{k},'file')
    rej = [rej k];
  end    
%   if ~exist(files{k},'file')
%     % add path to file name
%     files{k} = fullfile(parms.datapath,files{k});
%     if ~exist(files{k},'file'), rej = [rej k]; end
%   end
end
files(rej) = [];
clear rej

% select all files for the given events and channels
pat = [parms.prefix '(\S*)'];
if ~isempty(parms.events)
  tmp = sprintf('(%02i)',parms.events);
  tmp = [tmp sprintf('(%g)',parms.events)];
  tmp = strrep(tmp,')(',')|(');
  pat = [pat '(event)(' tmp ')+'];
end
if ~isempty(parms.channels)
  tmp = sprintf('(%03i)',parms.channels);
  tmp = [tmp sprintf('(%g)',parms.channels)];
  tmp = strrep(tmp,')(',')|(');
  pat = [pat '_(chan)(' tmp ')'];
%   pat = [pat '_(chan)(' tmp ').mat'];
  parms.num_sensors = length(parms.channels);
  pat = [pat '.mat'];
else
  pat = [pat '(\S*)chan(\S*).mat'];
end

if ~isempty(pat)
  ix   = regexp(files,pat);
  ix   = ~cellfun('isempty',ix);
  files = files(ix);
end
clear tmp pat ix

% did we find any files?
if isempty(files)
  error('No files found');
end

% Determine the number of channels and number of trials per condition and
% preallocate TF data matrices
rej     = [];
ntrials = [];
nchan   = 0;
nfreq   = 0;
F       = [];
T       = [];
for k = 1:length(parms.events)
  pat = [parms.prefix '(\S*)'];
  tmp = [sprintf('(%02i)',parms.events(k)) sprintf('(%g)',parms.events(k))];
  tmp = strrep(tmp,')(',')|(');
  pat = [pat '(event)(' tmp ')+'];
  sel = find(~cellfun('isempty',regexp(files,pat)));
  cnt = length(sel);
  if cnt==0
    rej = [rej k]; 
  else
    nchan = max(nchan,cnt);
    if isempty(parms.num_sensors), parms.num_sensors = nchan; end
    res = load(files{sel(1)},'timefreq_data');
    res.timefreq_data = ts_data_selection(res.timefreq_data,'foi',parms.foi,'foilim',parms.foilim,'toi',parms.toi,'toilim',parms.toilim);
    if k == 1
      timefreq_data = res.timefreq_data;
      complex_flag  = isfield(timefreq_data.timefreq,'cmplx');
      power_flag    = isfield(timefreq_data.timefreq,'power');
      F     = timefreq_data.timefreq.frequencies;
      T     = timefreq_data.timefreq.time;
      ntime = length(T);
      nfreq = length(F);
    else
      timefreq_data.timefreq(k) = res.timefreq_data.timefreq;
    end
    ntrials = [ntrials timefreq_data.timefreq(k).num_trials];
    if complex_flag
      tmp1 = zeros(nchan,ntime,nfreq,ntrials(k),'single');
      tmp2 = zeros(nchan,ntime,nfreq,ntrials(k),'single');
      timefreq_data.timefreq(k).cmplx = complex(tmp1,tmp2); 
      clear tmp1 tmp2
    end
    if power_flag, timefreq_data.timefreq(k).power   = zeros(nchan,ntime,nfreq,ntrials(k),'single'); end    
    if parms.verbose, mmil_logstr(parms,'Preallocating event %g: %g channels, %g frequencies, %g trials, %g samples',timefreq_data.timefreq(k).event_code,nchan,nfreq,ntrials(end),ntime); end
  end
end
parms.events(rej) = [];

% % preallocate data matrix if necessary info is provided
% if nchan ~=0 && ~isempty(ntrials)
%   fprintf('Processing %g channels, %g frequencies, and %g conditions (%s trials, %g samples each)\n',nchan,nfreq,length(ntrials),num2str(ntrials),ntime);
%   for k = 1:length(ntrials)
%     if complex_flag
%       tmp1 = zeros(nchan,ntime,nfreq,ntrials(k),'single');
%       tmp2 = zeros(nchan,ntime,nfreq,ntrials(k),'single');
%       timefreq_data.timefreq(k).cmplx = complex(tmp1,tmp2); 
%       clear tmp1 tmp2
%     end
%     if power_flag, timefreq_data.timefreq(k).power   = zeros(nchan,ntime,nfreq,ntrials(k),'single'); end
%   end
% end
  
% Load all TF data from selected files
event_codes = [];
chan_labels = {};
ncond       = 1;
for k = 1:length(files)
%   if parms.verbose, mmil_logstr(parms,'Loading file %g of %g: %s',k,length(files),files{k}); end
  res = load(files{k},'timefreq_data');
  % check that file contains data satisfying the necessary criteria
  if isempty(fieldnames(res)) % file is missing timefreq_data
    continue;
  elseif res.timefreq_data.num_sensors > 1
    error('File must contain exactly 1 sensor: %s',files{k});
  elseif length(res.timefreq_data.timefreq) > 1
    error('File must contain exactly 1 condition: %s',files{k});
  end
  complex_flag = isfield(res.timefreq_data.timefreq,'cmplx');
  power_flag   = isfield(res.timefreq_data.timefreq,'power');
  % initialize output structure
  if k == 1             
    chan_labels   = {res.timefreq_data.sensor_info.label};
    event_codes   = [res.timefreq_data.timefreq.event_code];
    if isempty(F)
      timefreq_data = res.timefreq_data;
      F       = timefreq_data.timefreq.frequencies;
      T       = timefreq_data.timefreq.time;
      ntrials = [ntrials timefreq_data.timefreq.num_trials];
      ntime   = length(T);
      nfreq   = length(F);
      Fix     = 1:nfreq;
      Tix     = 1:ntime;
    else
      [jnk1,Fix,jnk2] = intersect(res.timefreq_data.timefreq.frequencies,F);
      [jnk1,Tix,jnk2] = intersect(res.timefreq_data.timefreq.time,T);
    end
  end
  % what channel & condition does this file contain?
  label = res.timefreq_data.sensor_info.label;
  code  = res.timefreq_data.timefreq.event_code;
  % is this a new channel?
  if ~ismember(label,chan_labels)
    chan_labels = {chan_labels{:} label};
    timefreq_data.sensor_info(end+1) = res.timefreq_data.sensor_info;
%     timefreq_data.sensor_info(k) = res.timefreq_data.sensor_info;
%     timefreq_data.num_sensors        = timefreq_data.num_sensors + 1;
  end
%   timefreq_data.sensor_info(k) = res.timefreq_data.sensor_info;
  % is this a new condition?
  if ~ismember(code,event_codes) % yes
    ntrials     = [ntrials res.timefreq_data.timefreq.num_trials];
    event_codes = [event_codes code];
    ncond       = length(event_codes);
    timefreq_data.timefreq(ncond) = res.timefreq_data.timefreq;
    timefreq_data.timefreq(ncond).time        = T;
    timefreq_data.timefreq(ncond).frequencies = F;
    % preallocate memory if we know how many channels there are
    if ~isempty(parms.num_sensors)
      if complex_flag
        timefreq_data.timefreq(ncond).cmplx = single(complex(zeros(parms.num_sensors,ntime,nfreq,ntrials(ncond)),...
                                                             zeros(parms.num_sensors,ntime,nfreq,ntrials(ncond))));
      end
      if power_flag
        timefreq_data.timefreq(ncond).power = zeros(parms.num_sensors,ntime,nfreq,ntrials(ncond),'single');
      end
    end
  end
  % get indices to this channel & condition
  chan = strmatch(label,chan_labels,'exact');
  cond = find(code == event_codes);
  % add this data to timefreq_data
  if complex_flag, timefreq_data.timefreq(cond).cmplx(chan,:,:,:) = res.timefreq_data.timefreq.cmplx(1,Tix,Fix,:); end
  if power_flag  , timefreq_data.timefreq(cond).power(chan,:,:,:) = res.timefreq_data.timefreq.power(1,Tix,Fix,:); end  
end
timefreq_data.num_sensors = length(timefreq_data.sensor_info);