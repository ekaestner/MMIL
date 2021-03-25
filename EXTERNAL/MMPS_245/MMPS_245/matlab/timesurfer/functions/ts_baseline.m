function baseline_data = ts_baseline(data,varargin)

% toi & toilim units must match what is in the data structure

parms = mmil_args2parms(varargin,...
						{'events',[],[],...
						 'conditions',[],[],...
						 'toilim',[],[],...
             'toi',[],[],...
             'foi',[],[],...
             'keep_event_codes',1,[],...
             'out','yes',{'yes','avg_data','epoch_data','timefreq_data','timefreq'},...
						},false);
          
[datatype,datafield,dataparam]   = ts_object_info(data,varargin{:});

if strcmp(datatype,'avg_data') || (strcmp(datatype,'timefreq_data') && ndims(data.(datafield)(1).(dataparam{1})) == 3)
  trials_flag = 0;
elseif strcmp(datatype,'epoch_data') || (strcmp(datatype,'timefreq_data') && ndims(data.(datafield)(1).(dataparam{1})) == 4)
  trials_flag = 1;
else
  error('data must be avg_data, epoch_data, or timefreq_data'); 
end

% extract events
if isempty(parms.events) && isempty(parms.conditions)
  conds = 1:length(data.(datafield));
elseif ~isempty(parms.conditions)
  conds = parms.conditions;
elseif ~isempty(parms.events)
  conds = find(ismember([data.(datafield).event_code],parms.events));
end
if isempty(conds)
  error('could not find event codes');
end
data.(datafield) = data.(datafield)(conds);

if trials_flag && strcmpi(datatype,'epoch_data') && strcmpi(parms.out,'timefreq')
  % convert epoch_data to timefreq_data
  data = ts_freqanalysis_fieldtrip(data,varargin{:});
  [datatype,datafield,dataparam]   = ts_object_info(data,varargin{:});
  trials_flag = 0;
end
if trials_flag
  avg_data = ts_trials2avg(data); % average over trials
else
  avg_data = data;                % already an average over trials
  if isfield(avg_data,datafield)
    avg_data.averages = avg_data.(datafield);
    avg_data = rmfield(avg_data,datafield);
  end
end
baseline_data = rmfield(data,datafield);
baseline_data.(datafield) = data.(datafield)(1);
if parms.keep_event_codes
  baseline_data.(datafield).event_code = [data.(datafield).event_code];
else
  baseline_data.(datafield).event_code = 1;
end
clear data;

% toilim
if ~isempty(parms.toi)
  [c tidx ib] = intersect(avg_data.averages(1).time,parms.toi);
else
  if isempty(parms.toilim) || length(parms.toilim)~=2
    parms.toilim = [avg_data.averages(1).time(1) avg_data.averages(1).time(end)];
  end
  tidx = nearest(avg_data.averages(1).time,parms.toilim(1)):...
         nearest(avg_data.averages(1).time,parms.toilim(end));  
end
baseline_data.(datafield).time = avg_data.averages(1).time(tidx);

for p = 1:length(dataparam)
  try
    baseline_data.(datafield).(dataparam{p}) = [];
    ndim = 2;%ndims(avg_data.averages(1).(dataparam{p}));
    dat  = {avg_data.averages.(dataparam{p})};
    N    = [];
    S    = 0;
    for c = 1:length(dat)
      N = cat(1,N,avg_data.averages(c).num_trials);
      S = S + N(end)*dat{c}(:,tidx,:);
    end
    % weighted average
    baseline_data.(datafield).(dataparam{p}) = S ./ sum(N);
  end
end
baseline_data.(datafield).num_trials = sum(N);
