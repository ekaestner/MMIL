function [outdata] = ts_timefreq_average(varargin)

% freq_bands: cell array of numeric array of band limits
%   (ex. {[2 8] [8 12] [12 25]})
%
% mods:
%   09-Sep-2009: added rejectfile option
% todo:
%   - move outfile to beginning and abort if exists and ~overwrite
%   - finish mods to work with data struct inputs (instead of file lists)
%
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
     'blc',[],[],...
     'baselinetype','zscore',[],...
     'baselinefile',[],[],...
     'blcwindow',[-inf 0],[],...
     'baseline_data',[],[],...
     'freqband','all',[],...
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
     'avgovertrials',0,{0,1},...
    },...
    false);
if (isnumeric(parms.blc) && ~isempty(parms.blc) && parms.blc==1) || (ischar(parms.blc) && ~strcmpi(parms.blc,'no')), parms.baseline = 'yes'; end
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
if ischar(parms.freqband) && strcmp(parms.freqband,'all')
  parms.freqband = [data.(datafield)(1).frequencies(1) data.(datafield)(1).frequencies(end)];
end
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
  if ~iscell(parms.freqband),   parms.freqband   = {parms.freqband};       end
else
  if isempty(parms.conditions), parms.conditions = 1:length(data.(datafield));                      end
  if isempty(parms.events),     parms.events     = [data.(datafield)(parms.conditions).event_code]; end
  if ~iscell(parms.freqband),   parms.freqband   = {parms.freqband};                                end  
end

if ischar(parms.baselinefile) && exist(parms.baselinefile,'file')
  mmil_logstr(parms,'%s: loading baseline data: %s\n',mfilename,parms.baselinefile);
%   fprintf('%s: loading baseline data: %s\n',mfilename,parms.baselinefile);
  load(parms.baselinefile);
  parms.baseline_data = baseline_data;
  clear baseline_data;
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
if isequal(parms.baseline_data,1) || isequal(parms.baseline_data,'yes')
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
      dat   = ts_data_selection(dat,'toilim',parms.blcwindow);
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
      ntrl  = dat.timefreq.num_trials;
      if c==1 && ch==1
        trials = 1:ntrl;
        parms.baseline_data = dat;
        parms.baseline_data.timefreq.power      = zeros(nchan,length(dat.timefreq(1).time),length(dat.timefreq(1).frequencies),ntrl);
      elseif ch==1
        trials = parms.baseline_data.timefreq.num_trials + [1:ntrl];
        parms.baseline_data.timefreq.num_trials = parms.baseline_data.timefreq.num_trials + ntrl;
        parms.baseline_data.timefreq.power(:,:,:,trials) = zeros(nchan,length(dat.timefreq(1).time),length(dat.timefreq(1).frequencies),ntrl);
      end
      if c==1 && ch>1
        parms.baseline_data.sensor_info(ch) = dat.sensor_info;
        parms.baseline_data.num_sensors     = parms.baseline_data.num_sensors + 1;
      end
      parms.baseline_data.timefreq.power(ch,:,:,trials) = dat.timefreq.power;
      cnt = cnt + 1;
    end
  end
end
%% calculate frequency band averages
fcount = 0;
for b = 1:length(parms.freqband)
  band = parms.freqband{b};
  if length(band) > 2
    band = [min(band) max(band)];
  end
  fcount = 0;
  % load/select and concatenate timefreq_data structures
  for c = 1:length(parms.conditions)
    for ch = 1:length(chans)
      % select file for this condition and channel
      if loadflag
        fcount = fcount + 1;
        fname  = datafile{ch,c};
        if parms.verbose
          mmil_logstr(parms,'%s: processing data file %g of %g: %s\n',mfilename,fcount,nfiles,fname);
%           fprintf('%s: processing data file %g of %g: %s\n',mfilename,fcount,nfiles,fname);
        end
        % load timefreq_data
        S   = load(fname);
        s   = fieldnames(S);
        dat = S.(s{1});
        clear S s
      else
        dat = ts_data_selection(data,'channels',ch,'events',data.(datafield)(c).event_code,'removebadchans',1);
      end
      % remove unecessary data to save memory
      if isfield(dat.timefreq,'data')
        dat.timefreq = rmfield(dat.timefreq,'data');
      end
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
      % correct for alternative timefreq_data format
      if parms.ai_struct
        dat.num_sensors = 1;
        dat.sensor_info = dat.sensor_info(chans(ch));
      end
      if ~isempty(parms.toilim)
        dat = ts_data_selection(dat,'toilim',parms.toilim);
      end
      % SK: remove bad trials but keep all channels, we will get rid off them later
      if ~isempty(parms.reject_data) && isstruct(parms.reject_data)
        dat = ts_data_selection(dat,'reject_data',parms.reject_data,'keepbadchans',1);
      end

      if ~(ischar(parms.baseline) && strcmp(parms.baseline,'no')) && isempty(parms.baseline_data)
        % create baseline data
        mmil_logstr(parms,'%s: averaging baseline over trials...\n',mfilename);
%         fprintf('%s: averaging baseline over trials...\n',mfilename);
        parms.baseline_data = rmfield(dat,{'timefreq','sensor_info'});
        parms.baseline_data.timefreq = rmfield(dat.timefreq,'power');
        parms.baseline_data.sensor_info = dat.sensor_info;
        bdat = ts_data_selection(dat,'toilim',parms.blcwindow);
        [sz1 sz2 sz3 sz4] = size(bdat.timefreq.power);
        parms.baseline_data.timefreq.power = zeros(1,1,sz3,sz4);
        parms.baseline_data.timefreq.power(1,:,:,:) = mean(bdat.timefreq.power,2);
        parms.baseline_data.timefreq.power = repmat(parms.baseline_data.timefreq.power,[1 length(dat.timefreq.time) 1 1]);
        parms.baseline_data.num_sensors = 1;
        parms.baseline_data.type = 'cond';
        tdim = 4;
      else
        tdim = 2;
      end;
      
      % whitening (correction for 1/f^2 scaling of power spectrum)
      if parms.freqcorr ~= 0
        mmil_logstr(parms,'%s: correcting for 1/f^2 scaling of power spectrum\n',mfilename);
%         fprintf('%s: correcting for 1/f^2 scaling of power spectrum\n',mfilename);
        freqcorr = permute(dat.timefreq.frequencies.^2,[1 3 2]); % 1 x 1 x freqs
        freqcorr = repmat(freqcorr,[1 length(dat.timefreq.time) 1 dat.timefreq.num_trials]);
        dat.timefreq.power = dat.timefreq.power .* freqcorr;
        clear freqcorr
      end
      % average over trials
      if parms.avgovertrials
        dat.timefreq.power = nanmean(dat.timefreq.power,4);
      end
      % baseline correction
      if exist('baseline_data','var')
        dat = ts_zscore(dat,'baseline_data',baseline_data,'blcwindow',parms.blcwindow,'zparam','power','baselinetype',parms.baselinetype,'tdim',tdim);
      elseif ~isempty(parms.baseline_data)
        dat = ts_zscore(dat,'baseline_data',ts_data_selection(parms.baseline_data,'chanlabel',{dat.sensor_info.label}),'blcwindow',parms.blcwindow,'zparam','power','baselinetype',parms.baselinetype,'tdim',tdim);
        if isfield(parms.baseline_data,'type') && ~isempty(parms.baseline_data.type) && strcmp(parms.baseline_data.type,'cond')
          parms.baseline_data = [];
        end
      elseif ~(ischar(parms.baseline) && strcmp(parms.baseline,'no'))
        dat = ts_zscore(dat,'blcwindow',parms.blcwindow,'zparam','power','baselinetype',parms.baselinetype,'tdim',tdim);
      end

%       % whitening (correction for 1/f^2 scaling of power spectrum)
%       if parms.freqcorr ~= 0
%         fprintf('%s: correcting for 1/f^2 scaling of power spectrum\n',mfilename);
%         freqcorr = permute(dat.timefreq.frequencies.^2,[1 3 2]); % 1 x 1 x freqs
%         freqcorr = repmat(freqcorr,[1 length(dat.timefreq.time) 1 dat.timefreq.num_trials]);
%         dat.timefreq.power = dat.timefreq.power .* freqcorr;
%         clear freqcorr
%       end

      % add spectral power to outdata structure
      if c == 1 && ch == 1
        outdata = rmfield(dat,'timefreq');
      end
      if ch == 1
        outdata.timefreq(c)       = dat.timefreq;
        outdata.timefreq(c).power = nan([length(chans) length(dat.timefreq.time) length(dat.timefreq.frequencies)]);
      elseif c == 1
        outdata.num_sensors = outdata.num_sensors + 1;
        outdata.sensor_info = cat(2,outdata.sensor_info,dat.sensor_info);
      end
      % calculate average over trials
      if ndims(dat.timefreq.power)==4
        outdata.timefreq(c).power(ch,:,:) = squeeze(nanmean(dat.timefreq.power,4));
      else
        outdata.timefreq(c).power(ch,:,:) = dat.timefreq.power;
      end
      clear dat
    end
  end
end
% save epoch_data
if isfield(outdata,'parms')
  parms.previous = outdata.parms;
end
outdata.parms = parms;

clear foi args
%SK get rid of bad channels but keep all trials
if ~isempty(parms.reject_data) && isstruct(parms.reject_data)
  parms.reject_data.badtrials = [];
  outdata = ts_data_selection(outdata,'reject_data',parms.reject_data,'removebadchans',1);
end

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

    


  
