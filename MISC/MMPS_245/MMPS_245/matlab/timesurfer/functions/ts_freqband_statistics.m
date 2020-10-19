function [filenames] = ts_timefreq_statistics(varargin)

% freq_bands: cell array of numeric array of band limits
%   (ex. {[2 8] [8 12] [12 25]})

parms = mmil_args2parms(varargin,...
    {'datafile'     , [],[],...
     'conditions'   , [],[],...
     'events'       , [],[],...
     'method'       , 'montecarlo',{'montecarlo','analytic','stats','glm'},...
     'statistic'    , 'indepsamplesT',{'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                        'indepsamplesZcoh','depsamplesT','depsamplesF' ,...
                        'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
                        'anova1','kruskalwallis'},...
     'correctm'   , 'cluster',{'no','max','cluster','bonferoni','fdr'},...
     'alpha'      , 0.05, [],...
     'tail'       ,    0, {-1, 1, 0},...
     'ivar'       ,   1,[],...
     'uvar'       ,   [],[],...
     'wvar'       ,   [],[],...
     'feedback'   , 'text',{'gui', 'text', 'textbar', 'no'},...
     'clusterstatistic', 'maxsum', {'maxsum','maxsize','wcm'},...
     'clusterthreshold', 'parametric', {'parametric','nonparametric'},...
     'clusteralpha',        0.05 , [],...
     'clustercrtival',        [],[],...
     'clustertail',0,{-1, 1, 0},...
     'channel'    , [],[],...
     'chantype',      'all', {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
     'latency'    , 'all',[],...
     'frequency', 'all',[],...
     'findex',[],[],...
     'avgoverchan', 'no',[],...
     'avgovertime', 'no',[],...
     'numrandomization',500,[],...
     'minnbchan',[],[],...
     'neighbours',[],[],...
     'design',[],[],...
     'cfg',[],[],...
     'baseline','no',[],...
     'baselinetype','zscore',{'absolute','relchange','relative','no'},...
     'baselinefile',[],[],...
     'blcwindow',[-inf 0],[],...
     'verbose',1,{1,0},...
     'freqband','all',[],...
     'skipstats_flag',1,{0,1},...
     'trials_flag',0,{0,1},...
     'freqcorr',0,{0,1},...
     'toilim',[],[],...
     'ai_struct',0,{0,1},...
     'filename',[],[],...
     'overwrite',0,{0,1},...
    },...
    false);
outfile = parms.filename;
if ~parms.verbose, parms.feedback = 'no'; end

% parms.datafile should contain a list of all files with TF trials
% timefreq_data in first datafile should contain sensor_info w/ all sensors
if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
hdr.parms.filename = parms.datafile{1};
if iscell(parms.events)
  data = ts_checkdata_header(hdr,'events',[parms.events{:}]);
else
  data = ts_checkdata_header(hdr,'events',parms.events);
end
[datatype datafield dataparam] = ts_object_info(data);
hdr  = rmfield(data,datafield);
if ischar(parms.freqband) && strcmp(parms.freqband,'all')
  parms.freqband = [data.(datafield)(1).frequencies(1) data.(datafield)(1).frequencies(end)];
end
clear data;
% if issubfield(hdr,'parms.filename')
%   parms.datafile = hdr.parms.filename;
% end

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

if isempty(parms.conditions), parms.conditions = 1:length(parms.events); end
if isempty(parms.events),     parms.events     = parms.conditions;       end
if ~iscell(parms.freqband),   parms.freqband   = {parms.freqband};          end

if length(parms.events) > 2 && ~parms.skipstats_flag
  error('%s: freqband statistics can only compare two conditions at this time.\n',mfilename);
end
if ischar(parms.baselinefile) && exist(parms.baselinefile,'file')
  fprintf('%s: loading baseline data: %s\n',mfilename,parms.baselinefile);
  load(parms.baselinefile);
end

[datafile chans] = sortfiles(parms.datafile,chans,parms.events);
nchan    = length(chans);
ncdfiles = size(datafile,1);
ncond    = size(datafile,2);
if parms.verbose, fprintf('%s: selecting %i files for %i channels (for each of %i conditions)\n',...
  mfilename,ncdfiles,nchan,ncond); end
nfiles   = ncdfiles*ncond;

% stat_data = [];
%% calculate frequency band averages
fcount = 0;
for b = 1:length(parms.freqband)
  band = parms.freqband{b};
  if length(band) > 2
    band = [min(band) max(band)];
  end
  % load and concatenate timefreq_data structures
  for c = 1:length(parms.conditions)
    for ch = 1:length(chans)
      % select file for this condition and channel
      fcount = fcount + 1;
      fname  = datafile{ch,c};
      if parms.verbose
        fprintf('%s: processing data file %g of %g: %s\n',mfilename,fcount,nfiles,fname);
      end
      % load timefreq_data
      S   = load(fname);
      s   = fieldnames(S);
      dat = S.(s{1}); 
      clear S s
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
%       if ~parms.trials_flag && ndims(dat.timefreq.power)==4
%         dat.timefreq.power = nanmean(dat.timefreq.power,4);
%       end
      % baseline correction
      % note: this correction is done for each trial
      % => may not make sense to use baseline_data (avg over trials & conds)
      if exist('baseline_data','var')
        dat = ts_zscore(dat,'baseline_data',baseline_data,'blcwindow',parms.blcwindow,'zparam','power','baselinetype',parms.baselinetype);
      elseif ~strcmp(parms.baseline,'no')
        dat = ts_zscore(dat,'blcwindow',parms.blcwindow,'zparam','power','baselinetype',parms.baselinetype);
%       elseif ~strcmp(parms.baseline,'no') && ~strcmp(parms.baselinetype,'no')
%         freq              = ts_data2fieldtrip(dat,'condition',1,'dimord','rpt_chan_freq_time','channels',parms.channel,'chantype',parms.chantype);
%         cfg               = [];
%         cfg.baseline      = parms.baseline;
%         cfg.baselinetype  = parms.baselinetype;
%         freq              = freqbaseline(cfg,freq);
%         dat               = ts_fieldtrip2data(freq,dat);
%         clear freq cfg
      end      

      % whitening (correction for 1/f^2 scaling of power spectrum)
      if parms.freqcorr ~= 0
        fprintf('%s: correcting for 1/f^2 scaling of power spectrum\n',mfilename);
        freqcorr = permute(dat.timefreq.frequencies.^2,[1 3 2]); % 1 x 1 x freqs
        freqcorr = repmat(freqcorr,[1 length(dat.timefreq.time) 1 dat.timefreq.num_trials]);
        dat.timefreq.power = dat.timefreq.power .* freqcorr;        
        clear freqcorr
      end
      
      % add spectral power to artificial epoch_data structure
      if c == 1 && ch == 1
        epoch_data = rmfield(dat,'timefreq');
      end
      if ch == 1
        dat.timefreq.data    = nan([length(chans) length(dat.timefreq.time) dat.timefreq.num_trials]);
        epoch_data.epochs(c) = rmfield(dat.timefreq,{'power','frequencies'});
      elseif c == 1 % ch > 1
        epoch_data.num_sensors = epoch_data.num_sensors + 1;
        epoch_data.sensor_info = cat(2,epoch_data.sensor_info,dat.sensor_info);          
%         pow = cat(1,epoch_data.epochs.power,dat.timefreq.power);
%       else
%         epoch_data.epochs.data = cat(1,epoch_data.epochs.power,dat.timefreq.power);
      end
      % calculate frequency band average and add to epoch_data
      if isempty(parms.findex)
        fidx1 = nearest(dat.timefreq.frequencies,band(1));
        fidx2 = nearest(dat.timefreq.frequencies,band(2));                
        foi   = dat.timefreq.frequencies(fidx1:fidx2);
        epoch_data.epochs(c).data(ch,:,:) = squeeze(nanmean(dat.timefreq.power(:,:,fidx1:fidx2,:),3));        
      else
        foi   = dat.timefreq.frequencies(parms.findex);
        epoch_data.epochs(c).data(ch,:,:) = squeeze(nanmean(dat.timefreq.power(:,:,parms.findex,:),3));        
      end
%       epoch_data.epochs(c).data(ch,:,:) = shiftdim(squeeze(nanmean(dat.timefreq.power(:,:,fidx1:fidx2,:),3)),-1);
      clear dat fidx1 fidx2
    end
  end
  % save epoch_data
  if isfield(epoch_data,'parms')
    parms.previous = epoch_data.parms;
  end
  epoch_data.parms = parms;
  tmpname = save_epochs(epoch_data,foi,outfile,parms.overwrite);
  if ~isempty(tmpname)
    if b==1
      parms.filename{1}     = tmpname;
    else
      parms.filename{end+1} = tmpname;
    end
    clear tmpname
  end
  if ~parms.skipstats_flag
    % call ts_statistics_wrapper
    args         = rmfield(parms,{'conditions','events','freqband'});
    args.events  = [epoch_data.epochs.event_code];
    args         = mmil_parms2args(args);
    stat_data(b) = ts_statistics_wrapper(epoch_data,args{:});
    [stat_data(b).stats(1:length(stat_data.stats)).frequencies] = deal(foi);
    stat_data(b).parms = parms;
  elseif b==length(parms.freqband)
    stat_data.parms = parms;
  end
  clear epoch_data foi args   
end
filenames = parms.filename;

function outfile = save_epochs(epoch_data,thisfoi,filename,overwrite)
if isempty(filename)
  return;
end
if iscell(filename)
  filename   = filename{1};
end
[fpath fname fext] = fileparts(filename);
foilim             = [thisfoi(1) thisfoi(end)];
if ~any(findstr(fname,'freq'))
  fname            = sprintf('%s_freq%g-%gHz_',foilim(1),foilim(2));
end
evcodes            = [epoch_data.epochs.event_code];
fname              = sprintf('%s_conds_%s',fname,strrep(num2str(evcodes),'  ','_'));
outfile            = fullfile(fpath,[fname '_avg' fext]);

epoch_data.parms.filename = {outfile};
if ~exist(outfile,'file') || overwrite
  avg_data = ts_trials2avg(epoch_data);
  fprintf('%s: saving file: %s\n',mfilename,outfile);
  save(outfile,'epoch_data','avg_data');
else
  outfile = [];
  fprintf('%s: not overwriting %s\n',mfilename,outfile);
end

function [datafile outchans] = sortfiles(matfiles,chans,conds)
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
    fprintf('%s: Warning: could not find a timefreq_data filename containing ''chan%03i''\n',mfilename,chans(ch));
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

    


  