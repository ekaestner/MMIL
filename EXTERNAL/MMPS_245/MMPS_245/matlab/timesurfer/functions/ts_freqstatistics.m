function stat_data = ts_freqstatistics(data,varargin)

parms = mmil_args2parms(varargin,...
    {'conditions'   , [],[],...
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
    },...
    false);

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

% single-channel case
if issubfield(data,'timefreq.power') && size(data.timefreq(1).power,1) == 1
  parms.channel = 1;
end
  
%% between-conditions
for c = 1:length(parms.conditions)
  % convert to fieldtrip
  freq{c} = ts_data2fieldtrip(data,'condition',parms.conditions(c),'dimord','rpt_chan_freq_time','channels',parms.channel,'chantype',parms.chantype);
  % baseline correction
  if ~strcmp(parms.baseline,'no') && ~strcmp(parms.baselinetype,'no')
    cfg = [];
    cfg.baseline = parms.baseline;
    cfg.baselinetype = parms.baselinetype;
    freq{c} = freqbaseline(cfg,freq{c});
  end
end

% set up cfg for stats
cfg = [];
cfg = rmfield(parms,{'events','conditions','baseline','baselinetype'});

% create design matrix
ntrials = 0;
for c=1:length(parms.conditions)
  ntrials = ntrials+size(freq{c}.powspctrm,1);
end
design = zeros(1,ntrials);
trl1 = 1;
for c=1:length(parms.conditions)
  trl2 = trl1 + size(freq{c}.powspctrm,1) - 1;
  design(trl1:trl2) = c;
  trl1 = trl2 + 1;
end		
cfg.design = design;		
cfg.channel = freq{1}.label;

% compute stats
ftstat1 = freqstatistics(cfg,freq{:}); 
clear freq;
  
%% within-condition  

for c = 1:length(parms.conditions)
  % convert to fieldtrip
  freq = ts_data2fieldtrip(data,'condition',parms.conditions(c),'dimord','rpt_chan_freq_time','channels',parms.channel,'chantype',parms.chantype);
  
  % baseline correction
  if ~strcmp(parms.baseline,'no') && ~strcmp(parms.baselinetype,'no')
    cfg               = [];
    cfg.baseline      = parms.baseline;
    cfg.baselinetype  = parms.baselinetype;
    freq              = freqbaseline(cfg,freq);
  end  
  
  % set up cfg for stats
  cfg           = [];
  cfg.method    = 'stats';
  cfg.statistic = 'ttest';
  cfg.alpha     = parms.alpha;
  cfg.tail      = parms.tail;
  cfg.latency   = parms.latency;
  cfg.frequency = parms.frequency;
%   ntrials=size(freq.powspctrm,1);
%   design = zeros(2,2*ntrials);
%   design(1,1:ntrials)=1;
%   design(1,ntrials+1:2*ntrials)=2;
%   design(2,1:ntrials)=[1:ntrials];
%   design(2,ntrials+1:2*ntrials)=[1:ntrials];
%   cfg.design  = design;
  design = [];
  for tr=1:data.timefreq(c).num_trials
    design = cat(2,design,1);
  end		
  cfg.design = design;
  cfg.ivar    = 1;                         
  cfg.uvar    = 2;   
  cfg.design  = design;		
  cfg.channel = freq.label;
  
  ftstat2(c)    = freqstatistics(cfg,freq);
  clear freq;
end

%% data-reorganization (convert to stat_data)
evcodes = [-1 [data.timefreq(parms.conditions).event_code]];
ntrials = [sum([data.timefreq(parms.conditions).num_trials]) ...
          [data.timefreq(parms.conditions).num_trials]];
        
stat_data.num_sensors  = length(data.sensor_info);
stat_data.sensor_info  = data.sensor_info;
stat_data.sfreq        = data.sfreq;
clear data;

% quick fix (may not work for multiple channels)
for i = 1:length(ftstat1)+length(ftstat2)
  if i==1
    stat = ftstat1;
  else
    stat = ftstat2(i-1);
  end	
  [a,gchs]  = intersect({stat_data.sensor_info.label}, stat.label);                     
  [a,bchs]  = setdiff  ({stat_data.sensor_info.label}, stat.label);
  gch       = sort(gchs);
  bchs      = sort(bchs); 
  stat_data.stats(i).event_code    	= evcodes(i);
  stat_data.stats(i).prob(gch,:,:)  = permute(stat.prob(gch,:,:),[1 3 2]);
  stat_data.stats(i).prob(bchs,:,:)	= inf;            
  try stat_data.stats(i).stat(gch,:,:)  = permute(stat.stat(gch,:,:),[1 3 2]); end
  try stat_data.stats(i).stat(bchs,:,:)	= inf;     end
  try stat_data.stats(i).mask(gch,:,:)	= permute(stat.mask(gch,:,:),[1 3 2]); end
  try stat_data.stats(i).mask(bchs,:,:)	= 0;       end
  try stat_data.stats(i).time         = stat.time; end
	try stat_data.stats(i).frequencies 	= stat.freq; end
%   stat_data.stats(i).parms.method     = stat.cfg.method;
%   stat_data.stats(i).parms.statistic  = stat.cfg.statistic;
%   stat_data.stats(i).parms.alpha      = stat.cfg.alpha;
  stat_data.stats(i).num_trials       = ntrials(i);
	stat_data.stats(i).dimord						= stat.dimord;
	stat_data.stats(i).label						= stat.label;
	stat_data.stats(i).cfg							= stat.cfg;
  if isfield(stat,'posclusters'),         stat_data.stats(i).posclusters         = stat.posclusters; end
  if isfield(stat,'negclusters'),         stat_data.stats(i).negclusters         = stat.negclusters; end
  if isfield(stat,'posclusterslabelmat'), stat_data.stats(i).posclusterslabelmat = permute(stat.posclusterslabelmat,[1 3 2]); end
  if isfield(stat,'negclusterslabelmat'), stat_data.stats(i).negclusterslabelmat = permute(stat.negclusterslabelmat,[1 3 2]); end  
  if isfield(stat,'posdistribution'),     stat_data.stats(i).posdistribution     = stat.posdistribution; end
  if isfield(stat,'negdistribution'),     stat_data.stats(i).negdistribution     = stat.negdistribution; end
end
clear ftstat1 ftstat2 stat  
return;
  
