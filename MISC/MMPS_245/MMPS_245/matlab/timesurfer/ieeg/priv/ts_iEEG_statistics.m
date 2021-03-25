function [stat_data] = ts_iEEG_statistics (data_in,varargin)

% Use stats = ts_iEEG_statistics (epoch_data, 'option', value,...
%
% Function to perform statistical analysis using FieldTrip's
% timelockstatistics.
%
% Required input:
%
%   data_in    - epoch_data TimeSurfer data that contains all the conditions
%                necessary for stats requested
%   conditions - list of condition numbers (not event codes) to include in
%      or        analysis
%    events    - list of event_codes for conditions to include in analysis
%
%   method     - FieldTrip method for statistical analysis
%   statistic  - The statistics to produce for a given method
%
% Optional Inputs: (see FieldTrip documentation as each may vary based on
%                   method and statistic specified)
%
%   correctm, alpha, tail, ivar, uvar, wvar, feedback, clusterstatistic,
%   clusterthreshold, clusteralpha, clustercrtival, clustertail, channel,
%   latency,avgoverchan,avgovertime
%
% Ouput: 
%
%   stat_data  - a statistics data structure containing the
%                 following fields:
%      num_sensors - number of channels
%      sensor_info - information on the channels
%      sfreq       - sampling frequency
%      stats        - a structure array with following s
%        condition  - the condition number in relation to the original data (-1 for between conditions)
%        event_code - the event code (-1 for between conditions)
%        prob       - the result of the test p-value - a matrix [channels x time points]
%        mask       - mask (0 for bad channels)
%        time       - vector of time points 
%        parms      - method - the trip method used
%                     statistic - the trip statistic used
%                     alpha  - alpha value used
%        num_trials - number of trials
%        norm_bl    - result of Kolmogorov-Smirnov test for baseline data (0 for between conditions)
%        type       - 'Between-Conditions' or 'Within-Conditions'
%        posclusters         - a structure array with statistical results for positive clusters (sorted by p-value)
%                            - prob - p-values
%                            - clusterstat - cluster-level statistics
%        negclusters         - a structure array with statistical results for negative clusters (sorted by p-value)
%                            - prob - p-values
%                            - clusterstat - cluster-level statistics
%        posclusterslabelmat - matrix [channels x time points] indicating pos clusters to which (channel,time)-pairs belong
%        negclusterslabelmat - matrix [channels x time points] indicating neg clusters to which (channel,time)-pairs belong
%
%  See also: ts_iEEG_StatAnalysis, timelockstatistics, freqstatistics
%
% Created on ??? by Rajan Patel
% Last Modified on 07/17/08 by Jason Sherfey

%% Check Options

if ~mmil_check_nargs(nargin, 4), return; end;

opt = mmil_args2parms(varargin,...
                      {'conditions' , [],[],...
                       'events'     , [],[],...
                       'method'     , [],{'montecarlo','analytic','stats','glm'},...
                       'statistic'  , [],{'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                                          'indepsamplesZcoh','depsamplesT','depsamplesF' ,...
                                          'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
                                          'anova1','kruskalwallis'},...
                       'correctm'   , 'no',{'no','max','cluster','bonferoni','fdr'},...
                       'parameter', 'individual',[],...
                       'alpha'      , 0.05, [],...
                       'tail'       ,    0, {-1, 1, 0},...
                       'ivar'       ,   1,[],...
                       'uvar'       ,   [],[],...
                       'wvar'       ,   [],[],...
                       'feedback'   , 'textbar',{'gui', 'text', 'textbar', 'no'},...
                       'clusterstatistic', 'maxsum', {'maxsum','maxsize','wcm'},...
                       'clusterthreshold', 'parametric', {'parametric','nonparametric'},...
                       'clusteralpha'    ,        0.05 , [],...
                       'clustercrtival'  ,        [],[],...
                       'clustertail',0,{-1, 1, 0},...
                       'channel'    , 'all',[],...
											 'chantype',      'all', {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                       'latency'    , 'all',[],...
                       'avgoverchan', 'no',[],...
                       'avgovertime', 'no',[],...
                       'numrandomization',500,[],...
                       'minnbchan',0,[],...
                       'neighbours',[],[],...
                       'design',[],[],...
                       'cfg',[],[],...
                       'lpfilter','no',{'yes','no'},...
                       'hpfilter','no',{'yes','no'},...
                       'bpfilter','no',{'yes','no'},...
                       'bsfilter','no',{'yes','no'},...
                       'lnfilter','no',{'yes','no'},...
                       'dftfilter','no',{'yes','no'},...
                       'medianfilter','no',{'yes','no'},...
                       'lpfreq',[],[],...
                       'hpfreq',[],[],...
                       'bpfreq',[],[],...
                       'lnfreq',[],[],...
                       'dftfreq',[],[],... 
                       'medianfiltord',[],[],...
                       'blc','no',{'yes', 'no'},...
                       'blcwindow',[],[],...
                       'detrend','no',{'yes','no'},...                       
                      },...
                      false);
                    
%% Check conditions

if ismember('epochs',fieldnames(data_in))
    datafield = 'epochs';
elseif ismember('timefreq',fieldnames(data_in))
    datafield = 'timefreq';
else
    fprintf('ERROR: Invalid data supplied for stats.\n');
    return;
end

if ~isempty(opt.events)                                              % if event codes given convert to condition numbers
    if ~isempty(opt.conditions)
        fprintf('WARNING: Event values given will be used instead of condition numbers.\n');
    end
    conditions = [];
    for i = 1:length(data_in.(datafield))                            % convert events to conditions
        if ismember(data_in.(datafield)(i).event_code,opt.events)
            conditions = cat(1,conditions,i);
        end
    end
    if length(conditions) ~= length(opt.events)                      % make sure found all events
        fprintf('ERROR: Invalid event codes given.\n');
        return;
    end
elseif ~isempty(opt.conditions)
    conditions = opt.conditions;   
end

auto_design = 1;
if ~isempty(opt.design), auto_design = 0; end;

if isempty(opt.latency)
    opt.latency = 'all';
end


%% Initialize Some Data

events                 = [data_in.(datafield)(conditions).event_code];
stat_data.sensor_info  = data_in.sensor_info;
stat_data.num_sensors  = length(data_in.sensor_info);
stat_data.sfreq        = data_in.sfreq;
num_trials             = [];

%% Processing steps

tla = [];
tla.feedback = 'no';
tla.lpfilter      = opt.lpfilter; 
tla.hpfilter      = opt.hpfilter; 
tla.bpfilter      = opt.bpfilter; 
tla.lnfilter      = opt.lnfilter; 
tla.dftfilter     = opt.dftfilter;
tla.medianfilter  = opt.medianfilter; 
tla.lpfreq        = opt.lpfreq;                      
tla.hpfreq        = opt.hpfreq;   
tla.bpfreq        = opt.bpfreq; 
tla.lnfreq        = opt.lnfreq; 
tla.dftfreq       = opt.dftfreq; 
tla.medianfiltord = opt.medianfiltord; 
tla.detrend       = opt.detrend; 
tla.blc           = opt.blc; 
tla.blcwindow     = opt.blcwindow; 
% Currently Unavailable Options 
% tla.padding       = opt.padding;
% tla.lpfiltord     = opt.lpfiltord; 
% tla.hpfiltord     = opt.hpfiltord; 
% tla.bpfiltord     = opt.bpfiltord; 
% tla.lnfiltord     = opt.lnfiltord; 
% tla.lpfilttype    = opt.lpfilttype; 
% tla.hpfilttype    = opt.hpfilttype;
% tla.bpfilttype    = opt.bpfilttype; 
% tla.lpfiltdir     = opt.lpfiltdir; 
% tla.hpfiltdir     = opt.hpfiltdir; 
% tla.bpfiltdir     = opt.bpfiltdir; 
% tla.hilbert       = opt.hilbert; 
% tla.rectify       = opt.rectify; 
% tla.precision     = opt.precision; 
% tla.polyremoval   = opt.polyremoval; 
% tla.polyorder     = opt.polyorder; 

%% Run The Between Statistics
 %% Convert to Field Trip

 
fprintf('Converting to FieldTrip:\n');
if auto_design, opt.design       = []; end
for i=1:length(conditions)
  fprintf('Condition %s.\n',num2str(conditions(i)));
  switch datafield
      case 'epochs'
        ft_data        = ts_data2fieldtrip(data_in,'condition',conditions(i),'chantype',opt.chantype,'dimord','chan_time');
        tla.keeptrials = 'yes';
        warning off all
        %% process the data because epoch_data is does not have any filters
        %% applied
        tl_data {i}    = timelockanalysis(tla,ft_data);
        warning on all
        if isfield(tl_data{i},'dof'), tl_data {i} = rmfield(tl_data{i},'dof'); end
        clear ft_data
      case 'timefreq'
        %%% CHECK ME!
        tl_data {i}    = ts_data2fieldtrip(data_in,'condition',conditions(i),'chantype',opt.chantype,'dimord','chan_freq_time');
  end
  if auto_design
    switch opt.statistic
        case {'indepsamplesT','indepsamplesF','indepsamplesregrT',...
              'indepsamplesZcoh','ttest','ttest2','paired-ttest','anova1','kruskalwallis'}
         if strcmpi(datafield,'epochs')
            opt.design(1,end+1:end+data_in.(datafield)(i).num_trials) = i*ones(1,data_in.(datafield)(i).num_trials);
         elseif strcmpi (datafield,'timefreq')
            %% CHECK ME!
            opt.design(i) = i;
         end
        case {'depsamplesT','depsamplesF','depsamplesregrT'}
            error('You must supply your own design matrix for depenedent tests as a stat.design option in the setup file.');
    end
  end
end

 %% Run FieldTrip Analysis
fprintf('Testing between conditions.\n');
cfg = rmfield(opt,'conditions');
cfg = rmfield(cfg,'events');
if isstruct(opt.cfg)
   fields = fieldnames(opt.cfg);
   for f = fields
    cfg.(f) = opt.cfg.(f);
   end
end

switch datafield
    case 'epochs'
      ft_stats{1}   = timelockstatistics (cfg,tl_data{:});
    case 'timefreq'
      %% CHECK ME!
      cfg.parameter = 'powspctrm';
      ft_stats{1}   = freqstatistics (cfg,tl_data{:});
end
ft_stats{1}.type = 'Between-Conditions';
ft_stats{1}.norm_bl(1:length(data_in.sensor_info)) = 0;
num_trials(1) = 0;
event_codes   = events;

% save ft_stats ft_stats
%a=1:20;
%pd=pwd;
%display(pd);


%% Run the Individual Statistics vs Baseline

if length(conditions) > 1
  event_codes = [];
  event_codes(1) = -1;
  indiv.method     = 'stats';
  indiv.statistic  = 'ttest';
  indiv.feedback   = opt.feedback;
  indiv.channel    = opt.channel;
  indiv.latency    = opt.latency;
  indiv.design     = [];

  for i = 1:length(conditions)
    fprintf('Analyzing condition %s against baseline.\n',num2str(conditions(i)));
    switch datafield
        %%%%%%% AVG STATS
        case 'epochs'
            ft_data = ts_data2fieldtrip(data_in,'condition',conditions(i),'chantype',opt.chantype,'dimord','chan_time');
            for j=1:data_in.(datafield)(i).num_trials
                indiv.design = cat(2,indiv.design,1);
            end
            tla.keeptrials = 'yes';
            warning off all
            tl_data    = timelockanalysis(tla,ft_data);
            warning on all
            if isfield(tl_data,'dof'), tl_data = rmfield(tl_data,'dof'); end
            clear ft_data;

            ft_stats{i+1}     = timelockstatistics (indiv,tl_data);
            bl_samples = find(data_in.(datafield)(i).time >= opt.blcwindow(1) & data_in.(datafield)(i).time <= opt.blcwindow(2));
            for j=1:length(data_in.sensor_info)
                t = 0;
                if ~isnan(data_in.(datafield)(i).data(j,:,:))
                    t = kstest(data_in.(datafield)(i).data(j,bl_samples,:));
                end
                ft_stats{i+1}.norm_bl(j) = t;
            end
       
        %%%%%%%% TIME FREQ STATS
        case 'timefreq'
            %%% CHECK ME!
            tl_data        = ts_data2fieldtrip(data_in,'condition',conditions(i),'chantype',opt.chantype,'dimord','chan_freq_time');
            indiv.design   = 1;
            ft_stats{i+1}  = freqstatistics(indiv,tl_data);            
    end
    ft_stats{i+1}.type = 'Within-Conditions: baseline vs activation';    
    num_trials(i+1)   = data_in.(datafield)(i).num_trials;
    indiv.design = [];
    event_codes(i+1) = events(i);
    clear tl_data;
  end
end

num_trials(1) = sum([num_trials]);

clear data_in

%% Finalize Structure
   
for i = 1:length(ft_stats)
  [a,goodchans_idx] = intersect({stat_data.sensor_info.label}, ft_stats{i}.label);                      % Reintroduce bad channels
  [a,badchans_idx]  = setdiff  ({stat_data.sensor_info.label}, ft_stats{i}.label);
  goodchan_idx = sort(goodchans_idx);
  badchans_idx = sort(badchans_idx); 
  stat_data.stats(i).event_code           = event_codes(i);
  stat_data.stats(i).prob(goodchan_idx,:) = ft_stats{i}.prob;
  stat_data.stats(i).prob(badchans_idx,:) = inf;                % should be high because it is not sig
  stat_data.stats(i).mask(goodchan_idx,:) = ft_stats{i}.mask;
  stat_data.stats(i).mask(badchans_idx,:) = 0; 
  stat_data.stats(i).time                 = ft_stats{i}.time;
  stat_data.stats(i).parms.method         = ft_stats{i}.cfg.method;
  stat_data.stats(i).parms.statistic      = ft_stats{i}.cfg.statistic;
  stat_data.stats(i).parms.alpha          = ft_stats{i}.cfg.alpha;
  stat_data.stats(i).num_trials           = num_trials(i);
  stat_data.stats(i).norm_bl              = ft_stats{i}.norm_bl;
  stat_data.stats(i).type                 = ft_stats{i}.type;
  if isfield(ft_stats{i},'posclusters'),         stat_data.stats(i).posclusters         = ft_stats{i}.posclusters; end
  if isfield(ft_stats{i},'negclusters'),         stat_data.stats(i).negclusters         = ft_stats{i}.negclusters; end
  if isfield(ft_stats{i},'posclusterslabelmat'), stat_data.stats(i).posclusterslabelmat = ft_stats{i}.posclusterslabelmat; end
  if isfield(ft_stats{i},'negclusterslabelmat'), stat_data.stats(i).negclusterslabelmat = ft_stats{i}.negclusterslabelmat; end  
  if isfield(ft_stats{i},'posdistribution'),     stat_data.stats(i).posdistribution     = ft_stats{i}.posdistribution; end
  if isfield(ft_stats{i},'negdistribution'),     stat_data.stats(i).negdistribution     = ft_stats{i}.negdistribution; end
  if isfield(ft_stats{i},'stat'),                stat_data.stats(i).stat                = ft_stats{i}.stat; end
  if isfield(ft_stats{i},'dimord'),              stat_data.stats(i).dimord              = ft_stats{i}.dimord; end
  if isfield(ft_stats{i},'label'),               stat_data.stats(i).label               = ft_stats{i}.label; end
  if isfield(ft_stats{i},'cfg'),                 stat_data.stats(i).cfg                 = ft_stats{i}.cfg; end
  end

% % added by Thomas on July 4, 2008 so stats output includes all fields
%   stat_data.stats(1).everything           = ft_stats{1};
%   stat_data.stats(2).everything           = ft_stats{2};
%   stat_data.stats(3).everything           = ft_stats{3};
%   display('adding more fields to between-condition comparison results...');

