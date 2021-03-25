function stat_data = ts_timefreq_statistics(varargin)

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
     'avgoverchan', 'no',[],...
     'avgovertime', 'no',[],...
     'numrandomization',500,[],...
     'minnbchan',[],[],...
     'neighbours',[],[],...
     'design',[],[],...
     'cfg',[],[],...
     'baseline','no',[],...
     'baselinetype','no',{'absolute','relchange','relative','no'},...
     'verbose',1,{1,0},...
    },...
    false);

if ~parms.verbose, parms.feedback = 'no'; end
if strcmp(parms.method,'montecarlo'), needtrials = 1; else needtrials = 0; end

if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
hdr.parms.filename = parms.datafile{1};
if iscell(parms.events)
  data = ts_checkdata_header(hdr,'events',[parms.events{:}]);
else
  data = ts_checkdata_header(hdr,'events',parms.events);
end
[datatype datafield dataparam] = ts_object_info(data);
hdr  = rmfield(data,datafield);
clear data;
if issubfield(hdr,'parms.filename')
  parms.datafile = hdr.parms.filename;
end

% if ~exist(parms.datafile{1},'file')  
%   error('%s: datafile not found: %s\n',mfilename,parms.datafile{1});
% else
%   data = getfield(load(parms.datafile{1},'timefreq_data'),'timefreq_data');
%   hdr = rmfield(data,'timefreq');
%   clear data;    
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
if isempty(chans)
  error('%s: no channels selected',mfilename);
end;

if ~needtrials
%% compute stats from averages  
  % load data
  for f = 1:length(parms.datafile)
    if ~exist(parms.datafile{f},'file')
      error('%s: datafile not found: %s\n',mfilename,parms.datafile{f});
    end
    data(f) = getfield(load(parms.datafile{f},'timefreq_data'),'timefreq_data');
  end
  if f > 1
    data = ts_combine_data(data);
  end
  
  % conditions
  if isempty(parms.events) && isempty(parms.conditions)
    conds = 1:length(data.timefreq);
  elseif ~isempty(parms.conditions)
    conds = parms.conditions;
  elseif ~isempty(parms.events)
    conds = find(ismember([data.timefreq.event_code],parms.events));
  end
  data.timefreq = data.timefreq(conds);
  parms.conditions = 1:length(conds);
  
  args = mmil_parms2args(parms);
  stat_data = ts_freqstatistics(data,args{:});
  clear data;
else
%% compute stats from trials  
  parms.channel   = 1;
  parms.chantype  = [];
  if isempty(parms.conditions), parms.conditions = 1:length(parms.events); end
  if isempty(parms.events),     parms.events     = parms.conditions;       end
  datafile = sortfiles(parms.datafile,chans,parms.events);
  if parms.verbose, fprintf('%s: selecting %i files for %i channels (for each of %i conditions)\n',...
    mfilename,size(datafile,1),length(chans),size(datafile,2)); end
  n = 1;
  for ch = 1:length(chans)
    try if hdr.sensor_info(chans(ch)).badchan, continue; end; end
    for c = 1:length(parms.conditions)
      fname = datafile{ch,c};
      if parms.verbose, fprintf('%s: loading data file: %s\n',mfilename,fname); end
      tmpdat  = getfield(load(fname,'timefreq_data'),'timefreq_data');
      try tmpdat.timefreq = rmfield(tmpdat.timefreq,'data');  end
      if ~isfield(tmpdat.timefreq,'power') && isfield(tmpdat.timefreq,'cmplx')
          tmpdat.timefreq.power = 2*abs(double(tmpdat.timefreq.cmplx)).^2;
      elseif isa(tmpdat.timefreq.power,'single')
        tmpdat.timefreq.power = double(tmpdat.timefreq.power);
      end
      try tmpdat.timefreq = rmfield(tmpdat.timefreq,'cmplx'); end
      data(c) = tmpdat;
      clear tmpdat;
    end
    if c > 1
      data = ts_combine_data(data);
    end
    if length(data.sensor_info) > 1
      data.sensor_info = data.sensor_info(chans(ch));
      data.num_sensors = 1;
    end
    if data.sensor_info.badchan
      fprintf('%s: skipping bad channel: %s\n',mfilename,data.sensor_info.label); 
      continue; 
    end
    if parms.verbose, fprintf('%s: comparing events %g and %g, channel %g of %g - %s (%s) spectral power\n',...
                      mfilename,data.timefreq(1).event_code,data.timefreq(end).event_code,...
                      ch,length(chans),data.sensor_info.label,data.sensor_info.typestring); end
    args = mmil_parms2args(parms);
    chan_stat(n)   = ts_freqstatistics(data,args{:});
    clear data;
    n = n + 1;
  end
  % data re-organization
  if parms.verbose, fprintf('%s: re-organizing stat_data...',mfilename); end
  % concatenate stat structures
  for ch = 1:length(chan_stat)
    stat = chan_stat(ch);
    if ch==1
        stat_data=stat;
    else
      for i=1:length(stat.stats)
        stat_data.stats(i).prob=cat(1,stat_data.stats(i).prob,stat.stats(i).prob);
        stat_data.stats(i).mask=cat(1,stat_data.stats(i).mask,stat.stats(i).mask);
        stat_data.stats(i).stat=cat(1,stat_data.stats(i).stat,stat.stats(i).stat);
        try
          stat_data.stats(i).posclusters	=	cat(2,stat_data.stats(i).posclusters,stat.stats(i).posclusters);
          for k=1:length(stat.stats(i).posclusters)
            stat.stats(i).posclusterslabelmat(stat.stats(i).posclusterslabelmat==k) = k+length(stat_data.stats(i).posclusters);
          end
          stat_data.stats(i).posclusterslabelmat=cat(1,stat_data.stats(i).posclusterslabelmat,stat.stats(i).posclusterslabelmat);
          stat_data.stats(i).posdistribution=cat(1,stat_data.stats(i).posdistribution,stat.stats(i).posdistribution);
        end			
        try
          stat_data.stats(i).negclusters=cat(2,stat_data.stats(i).posclusters,stat.stats(i).negclusters);
          for k=1:length(stat.stats(i).negclusters)
            stat.stats(i).negclusterslabelmat(stat.stats(i).negclusterslabelmat==k)=k+length(stat_data.stats(i).negclusters);
          end
          stat_data.stats(i).negclusterslabelmat=cat(1,stat_data.stats(i).negclusterslabelmat,stat.stats(i).negclusterslabelmat);
          stat_data.stats(i).negdistribution=cat(1,stat_data.stats(i).negdistribution,stat.stats(i).negdistribution);
        end
%         try stat_data.stats(i).label = {stat_data.stats(i).label{:} stat.stats(i).label{:}}; end
        try stat_data.stats(i).label=cat(1,stat_data.stats(i).label,stat.stats(i).label); end
      end
      stat_data.sensor_info   = cat(1,stat_data.sensor_info,stat.sensor_info);
    end		
  end
  stat_data.num_sensors = length(stat_data.sensor_info);

  % sort clusters
  for i=1:length(stat_data.stats)
    % sort positive clusters by p-value
   try
    [sortp,sort_ind]=sort([stat_data.stats(i).posclusters.prob],'ascend');
    stat_data.stats(i).posclusters=stat_data.stats(i).posclusters(sort_ind);
    labelmattemp=stat_data.stats(i).posclusterslabelmat;
    for i=1:length(stat_data.stats(i).posclusters)
      labelmattemp(stat_data.stats(i).posclusterslabelmat==sort_ind(i))=i;
    end
    stat_data.stats(i).posclusterslabelmat=labelmattemp;
   end
   try
    % sort negative clusters by p-value
    [sortp,sort_ind]=sort([stat_data.stats(i).negclusters.prob],'ascend');
    stat_data.stats(i).negclusters=stat_data.stats(i).negclusters(sort_ind);
    labelmattemp=stat_data.stats(i).negclusterslabelmat;
    for i=1:length(stat_data.stats(i).negclusters)
      labelmattemp(stat_data.stats(i).negclusterslabelmat==sort_ind(i))=i;
    end
    stat_data.stats(i).negclusterslabelmat=labelmattemp;
   end
  end
  if parms.verbose, fprintf('done.\n'); end
end
return;
    

function datafile = sortfiles(matfiles,chans,conds)
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
  if isempty(tmpfiles), error('could not find a timefreq_data filename containing ''chan%03i''',chans(ch)); end
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
  datafile(ch,1:length(conds)) = tmp;
end  

    


  