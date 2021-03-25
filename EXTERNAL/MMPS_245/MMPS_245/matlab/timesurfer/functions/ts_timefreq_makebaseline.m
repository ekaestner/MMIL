function [baseline_data] = ts_timefreq_makebaseline(varargin)
% created by Jason Sherfey

if mod(nargin,2) && isstruct(varargin{1})
  % odd # of parameters AND 1st parm is data struct
  data = varargin{1};
  if nargin > 1
    varargin = varargin(2:end);
  else
    varargin = {};
  end
  data = ts_checkdata_header(data,varargin{:});
elseif mod(nargin,2)
  % odd # of parameters AND 1st parm is not data struct
  error('invalid specification of input parameters.');
else
  % even # of parameters (assume key/value pairs & search for data)
  scanflag = 1;
end

parms = mmil_args2parms(varargin,...
    {'datafile'     , [],[],...
     'conditions'   , [],[],...
     'events'       , [],[],...
     'channel'    , [],[],...
     'chantype',      'all', {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
     'frequency', 'all',[],...
     'findex',[],[],...
     'cfg',[],[],...
     'baseline','no',[],...
     'freqcorr',0,{0,1},...
     'toilim',[],[],...
     'ai_struct',0,{0,1},...
     'filename',[],[],...
     'overwrite',0,{0,1},...
     'rejectfile',[],[],...
     'reject_data',[],[],...
     'prefix',[],[],...
     'verbose',1,{0,1},...
     'logfile',      [],[],...
     'logfid',       [1],[], ... 
     'blcwindow',[],[],...
     'keeptrials','no',[],...
    },...
    false);
if isempty(parms.toilim) && ~isempty(parms.blcwindow), parms.toilim = parms.blcwindow; end
if isempty(parms.filename)
  outfile = parms.prefix;
else
  outfile = parms.filename;
end
if ~parms.verbose, parms.feedback = 'no'; end
if ~exist('data','var')
  % parms.datafile should contain a list of all files with TF trials
  % timefreq_data in first datafile should contain sensor_info w/ all sensors
  if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
  hdr.parms.filename = parms.datafile{1};
  data = ts_checkdata_header(hdr);
%   if iscell(parms.events)
%     data = ts_checkdata_header(hdr,'events',[parms.events{:}]);
%   else
%     data = ts_checkdata_header(hdr,'events',parms.events);
%   end
  loadflag = 1;
else
  loadflag = 0;
end
[datatype datafield dataparam] = ts_object_info(data);
hdr  = rmfield(data,datafield);
if loadflag, clear data; end

% channels (index to channels)
chans = []; 
if ~isempty(parms.channel)
  chans = parms.channel;
elseif ~isempty(parms.chantype)
  switch parms.chantype
    case {'mag' 'grad1' 'grad2' 'eeg', 'other'}
      chans     = find(strcmp(parms.chantype,{hdr.sensor_info.typestring}));
    case {'grad'}
      chans     = find(strncmp(parms.chantype,{hdr.sensor_info.typestring},length(parms.chantype)));
    case 'meg'
      [a,chans] = find(ismember({hdr.sensor_info.typestring},{'mag', 'grad1', 'grad2'}));
    case 'all'
      try chans = setdiff(1:hdr.num_sensors,find(strcmp('other',{hdr.sensor_info.typestring}))); end;
  end
else
  chans = 1:length(hdr.sensor_info);
end
chans = chans(~[hdr.sensor_info(chans).badchan]);
if isempty(chans)
  error('%s: no channels selected',mfilename);
end;

if loadflag
  if isempty(parms.conditions), parms.conditions = 1:length(parms.events); end
  if isempty(parms.events),     parms.events     = parms.conditions;       end
else
  if isempty(parms.conditions), parms.conditions = 1:length(data.(datafield));                      end
  if isempty(parms.events),     parms.events     = [data.(datafield)(parms.conditions).event_code]; end
end

if ischar(parms.rejectfile) && exist(parms.rejectfile,'file')
  mmil_logstr(parms,'%s: loading reject data: %s\n',mfilename,parms.rejectfile);
%   fprintf('%s: loading reject data: %s\n',mfilename,parms.rejectfile);
  load(parms.rejectfile);
  parms.reject_data = reject_data;
  clear reject_data;
end

if loadflag
  % sortfiles is a dirty little function that could be replaced by a few
  % lines using cellfun and cell arrays of filename substrings
  [datafile chans] = sortfiles(parms.datafile,chans,parms.events,parms);
  nchan    = length(chans);  
  ncdfiles = size(datafile,1);
  ncond    = size(datafile,2);
  mmil_logstr(parms,'%s: selecting %i files for %i channels (for each of %i conditions)\n',mfilename,ncdfiles,nchan,ncond);
  if parms.verbose, fprintf('%s: selecting %i files for %i channels (for each of %i conditions)\n',...
    mfilename,ncdfiles,nchan,ncond); end
  nfiles   = ncdfiles*ncond;
else
  nchan    = data.num_sensors;
  ncdfiles = 0;
  ncond    = length(data.(datafield));
  nfiles   = 0;
end
% if isequal(parms.baseline_data,1) || isequal(parms.baseline_data,'yes')
  fprintf('combining %g conditions for baseline_data\n',ncond);
  cnt   = 1;
  for c = 1:length(parms.conditions)
    for ch  = 1:length(chans)
      if parms.verbose, fprintf('%g of %g\n',cnt,ncond*nchan); end
      fname = datafile{ch,c};
      S     = load(fname);
      s     = fieldnames(S);
      dat   = S.(s{1});
      clear S s
      dat   = ts_data_selection(dat,'toilim',parms.toilim);
      % calculate power if necessary and convert to double precision
      if ~isfield(dat.timefreq,'power') && isfield(dat.timefreq,'cmplx')
        dat.timefreq.power = abs(double(dat.timefreq.cmplx)).^2;
      elseif isa(dat.timefreq.power,'single')
        dat.timefreq.power = double(dat.timefreq.power);
      end
      % remove complex spectra
      if isfield(dat.timefreq,'cmplx')
        dat.timefreq = rmfield(dat.timefreq,'cmplx');
      end      
      if isstruct(parms.reject_data)
        dat = ts_data_selection(dat,'reject_data',parms.reject_data);
      end
      ntrl  = dat.timefreq.num_trials;
      if c==1 && ch==1
        trials = 1:ntrl;
        parms.baseline_data = dat;
        parms.baseline_data.timefreq.power      = zeros(nchan,length(dat.timefreq(1).time),length(dat.timefreq(1).frequencies),ntrl);
      elseif ch==1
        trials = parms.baseline_data.timefreq.num_trials + [1:ntrl];
        parms.baseline_data.timefreq.num_trials = parms.baseline_data.timefreq.num_trials + ntrl;
        parms.baseline_data.timefreq.power(:,:,:,trials) = zeros(nchan,length(dat.timefreq.time),length(dat.timefreq(1).frequencies),ntrl);
        parms.baseline_data.timefreq.event_code = [parms.baseline_data.timefreq.event_code dat.timefreq.event_code];
      end
      if c==1 && ch>1
        parms.baseline_data.sensor_info(ch) = dat.sensor_info;
        parms.baseline_data.num_sensors     = parms.baseline_data.num_sensors + 1;
      end
      parms.baseline_data.timefreq.power(ch,:,:,trials) = dat.timefreq.power;
      cnt = cnt + 1;
    end
  end
% end
if isequal(parms.keeptrials,'no') || isequal(parms.keeptrials,0)
elseif isequal(parms.keeptrials,'cat')
  tmp = parms.baseline_data.timefreq.power;
  tmp = permute(tmp,[2 4 1 3]);
  tmp = [tmp(:)];
  ntm = length(parms.baseline_data.timefreq.time);
  nfq = length(parms.baseline_data.timefreq.frequencies);
  ntr = parms.baseline_data.timefreq.num_trials;
  nch = parms.baseline_data.num_sensors;
  tmp = reshape(tmp,[ntm*ntr,nch,nfq]);
  tmp = permute(tmp,[2 1 3]);
  parms.baseline_data.timefreq.power = tmp;
  clear tmp
  t = parms.baseline_data.timefreq.time;
  parms.baseline_data.timefreq.time = repmat(t,[1 ntr]);
  clear t
else
  parms.baseline_data.timefreq.power = nanmean(parms.baseline_data.timefreq.power,4);
end
baseline_data       = parms.baseline_data;
parms               = rmfield(parms,'baseline_data');
baseline_data.parms = parms;

function [datafile outchans] = sortfiles(matfiles,chans,conds,parms)
outchans = chans;
k = 1;
for ch = 1:length(chans)
  x  = regexp(matfiles,sprintf('chan%03i',chans(ch)));
  if isempty(x)
    x  = regexp(matfiles,sprintf('channel%03i',chans(ch)));
  end
  if isempty(x)
    x  = regexp(matfiles,sprintf('channel_%03i',chans(ch)));
  end  
  ix = []; 
  for i = 1:length(x)
    if ~isempty(x{i}),ix = [ix i]; end; 
  end
  tmpfiles = matfiles(ix);
  if isempty(tmpfiles)
    mmil_logstr(parms,'%s: Warning: could not find a timefreq_data filename containing ''chan%03i''\n',mfilename,chans(ch));
%     fprintf('%s: Warning: could not find a timefreq_data filename containing ''chan%03i''\n',mfilename,chans(ch));
    outchans = setdiff(outchans,chans(ch));
    continue;
  end
  tmp = {};
  for c = 1:length(conds)    
    x = regexp(tmpfiles,sprintf('cond%i_',conds(c)));
    ix = [];
    for i = 1:length(x)
      if ~isempty(x{i}),ix = [ix i]; end; 
    end
    if length(ix) > 1, error('found more than one timefreq_data filename containing ''chan%03i'' and ''cond%i''',chans(ch),conds(c)); end
    if ~isempty(ix), tmp(end+1) = tmpfiles(ix); end
  end
  if isempty(tmp)
    for c = 1:length(conds)    
      x = regexp(tmpfiles,sprintf('conds%i_',conds(c)));
      ix = [];
      for i = 1:length(x)
        if ~isempty(x{i}),ix = [ix i]; end; 
      end
      if length(ix) > 1, error('found more than one timefreq_data filename containing ''chan%03i'' and ''cond%i''',chans(ch),conds(c)); end
      if ~isempty(ix), tmp(end+1) = tmpfiles(ix); end
    end
  end
  if isempty(tmp)
    for c = 1:length(conds)
      x = regexp(tmpfiles,sprintf('event%i_',conds(c)));
      ix = [];
      for i = 1:length(x)
        if ~isempty(x{i}),ix = [ix i]; end; 
      end
      if length(ix) > 1, error('found more than one timefreq_data filename containing ''chan%03i'' and ''event%i''',chans(ch),conds(c)); end      
      tmp(end+1) = tmpfiles(ix);
    end
  end
  if isempty(tmp), error('could not find a timefreq_data filename containing ''cond%i'' or ''event%i''',conds(c),conds(c)); end
  if length(tmp) ~= length(conds), error('number of files not equal to number of conditions.'); end
  datafile(k,1:length(conds)) = tmp;
  k = k + 1;
end  

    


  
