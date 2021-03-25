function [stat_data] = ts_statistics (data_in,varargin)
%rmpath /home/mmildev/matlab/localmods/images/images
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

if ~mmil_check_nargs(nargin, 2), return; end;

opt = mmil_args2parms(varargin,...
                      {'conditions' , [],[],...
                       'events'     , [],[],...
                       'type'       , 'both',{'between','within','both'},...
                       'method'     , 'montecarlo',{'montecarlo','analytic','stats','glm'},...
                       'statistic'  , 'indepsamplesT',{'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                                          'indepsamplesZcoh','depsamplesT','depsamplesF' ,...
                                          'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
                                          'anova1','kruskalwallis'},...
                       'actvsbl_method','stats',{'montecarlo','analytic','stats','glm'},...
                       'actvsbl_statistic','ttest',{'ttest','ttest2','actvsblT'},...                   
                       'correctm'   , 'cluster',{'no','max','cluster','bonferoni','fdr'},...
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
                       'channel'    , [],[],...
											 'chantype',      'all', {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                       'latency'    , 'all',[],...
 											 'blcorrection','no',{'no','absolute','relative','relchange','zscore','zscores','yes',''},...
											 'bl_window',[],[],...		
                       'bl_tile',0,{0,1},...
                       'avgoverchan', 'no',[],...
                       'avgovertime', 'no',[],...
                       'numrandomization',500,[],...
                       'minnbchan',[],[],...
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
                       'removemean','no',{'yes','no'},...
                       'detrend','no',{'yes','no'},...       
                       'precision','double',{'double','single'},...
                      },...
                      false);

if iscell(opt.events)
  data_in = ts_checkdata_header(data_in,'precision',opt.precision,'events',[opt.events{:}]);
else
  data_in = ts_checkdata_header(data_in,'precision',opt.precision,'events',opt.events);
end
data_in = ts_data_selection(data_in,varargin{:});

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
tla.removemean     = opt.removemean; 
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

ft_stats = {}; num_trials = [];
%% Run The Between Statistics
 %% Convert to Field Trip

fprintf('converting to FieldTrip:\n');
if auto_design, opt.design       = []; end
if strcmpi(opt.type,'between') || strcmpi(opt.type,'both')
  for i=1:length(conditions)
    fprintf('condition %s.\n',num2str(conditions(i)));
    switch datafield
        case 'epochs'
          %ft_data        = ts_data2fieldtrip(data_in,'condition',conditions(i),'channels',opt.channel,'chantype',opt.chantype,'dimord','chan_time');
  ft_data        = ts_data2fieldtrip(data_in,'condition',conditions(i),'channels',opt.channel,'chantype',opt.chantype,'dimord','chan_time');
          tla.keeptrials = 'yes';
          warning off all
          %% process the data because epoch_data is does not have any filters
          %% applied
          tl_data {i}    = timelockanalysis(tla,ft_data);
          if ~strcmp(class(tl_data{i}.trial),opt.precision)
            fprintf('%s: converting trials to %s precision.\n',mfilename,opt.precision);
            eval(sprintf('tl_data{i}.trial = %s(tl_data{i}.trial);',opt.precision));
          end
          warning on all
          if isfield(tl_data{i},'dof'), tl_data {i} = rmfield(tl_data{i},'dof'); end
          clear ft_data
        case 'timefreq'
          %%% CHECK ME!
          tl_data {i}    = ts_data2fieldtrip(data_in,'condition',conditions(i),'channels',opt.channel,'chantype',opt.chantype,'dimord','chan_freq_time');
    end
    if auto_design
      switch opt.statistic
          case {'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                'indepsamplesZcoh','ttest','ttest2','paired-ttest','anova1','kruskalwallis'}
           if strcmpi(datafield,'epochs')
              opt.design(1,end+1:end+data_in.(datafield)(conditions(i)).num_trials) = i*ones(1,data_in.(datafield)(conditions(i)).num_trials);
           elseif strcmpi (datafield,'timefreq')
              %% CHECK ME!
              opt.design(i) = i;
           end
          case {'depsamplesT','depsamplesF','depsamplesregrT'}
              error('you must supply your own design matrix for depenedent tests as a stat.design option in the setup file.');
      end
    end
  end

   %% Run FieldTrip Analysis
  fprintf('testing between conditions.\n');
  cfg = rmfield(opt,'conditions');
  cfg = rmfield(cfg,'events');
  if isstruct(opt.cfg)
     fields = fieldnames(opt.cfg);
     for f = fields
      cfg.(f) = opt.cfg.(f);
     end
  end
  if isempty(cfg.channel), cfg.channel='all'; end
  switch datafield
      case 'epochs'
        % baseline correction
        if ismember(opt.blcorrection,{'absolute','relative','relchange','zscore','zscores','yes'})
          for c=1:length(tl_data)
            [ft_chidx jnk] = match_str(tl_data{c}.label,{data_in.sensor_info.label});
            % select baseline interval
            if isnumeric(opt.bl_window) && ~isempty(opt.bl_window) && (opt.bl_window(1) >= tl_data{c}.time(1)) && (opt.bl_window(2) <= 0)
              tidx = tl_data{c}.time >= opt.bl_window(1) & tl_data{c}.time <= opt.bl_window(2);
            else
              tidx = tl_data{c}.time <= 0;
            end
            fprintf('performing baseline correction: %s\n',opt.blcorrection);	                    
            if isempty(find(mean(tl_data{c}.trial(:,ft_chidx,tidx==1))==0,1))
              tidx = find(tidx);
              blmean = repmat(nanmean(tl_data{c}.trial(:,:,tidx),3),[1 1 length(tl_data{c}.time)]);
              blstd = repmat(nanstd(tl_data{c}.trial(:,:,tidx),[],3),[1 1 length(tl_data{c}.time)]);				
              switch opt.blcorrection
                case 'absolute'
                  tl_data{c}.trial = tl_data{c}.trial - blmean;
                case 'relative'
                  tl_data{c}.trial = tl_data{c}.trial ./ blmean;
                case 'relchange'
                  tl_data{c}.trial = (tl_data{c}.trial - blmean) ./ blmean;
                case {'zscore','zscores','yes'}
                  tl_data{c}.trial = (tl_data{c}.trial - blmean) ./ blstd;
              end
              clear blmean blstd	
            else
              warning('baseline correction was aborted. Average baseline power is zero!\n');
            end
          end
        end		
         cfg.channel = match_str(tl_data{1}.label,{data_in.sensor_info.label});
        ft_stats{end+1}   = timelockstatistics (cfg,tl_data{:});
      case 'timefreq'
        %% CHECK ME!
        cfg.parameter = 'powspctrm';
        ft_stats{end+1}   = freqstatistics (cfg,tl_data{:});
  end
  ft_stats{end}.norm_bl(1:length(data_in.sensor_info)) = 0;
  num_trials(end+1) = 0;
end
event_codes   = events;

%% Run the Individual Statistics vs Baseline

if strcmpi(opt.type,'within') || strcmpi(opt.type,'both')
  if isempty(ft_stats)
    cfg = rmfield(opt,'conditions');
    cfg = rmfield(cfg,'events');
    if isstruct(opt.cfg)
       fields = fieldnames(opt.cfg);
       for f = fields
        cfg.(f) = opt.cfg.(f);
       end
    end
  end
  event_codes = [];
  event_codes(1) = -1;
	if strcmpi(opt.actvsbl_method,'montecarlo')
		indiv.method     = opt.actvsbl_method;
	  indiv.statistic  = 'actvsblT';
		indiv.numrandomization = opt.numrandomization;
		indiv.clusterstatistic = opt.clusterstatistic;
		indiv.clusteralpha = opt.clusteralpha;
		indiv.alpha = opt.alpha;
		indiv.tail = cfg.tail;
		indiv.clustertail = cfg.clustertail;
		indiv.clusterstatistic = opt.clusterstatistic;
		indiv.minnbchan = cfg.minnbchan;
		indiv.neighbours = cfg.neighbours;
	else
		indiv.method     = opt.actvsbl_method;%'stats';
		indiv.statistic  = opt.actvsbl_statistic;%'ttest';
	end
  indiv.feedback   = opt.feedback;
  indiv.channel    = cfg.channel;
  indiv.latency    = opt.latency;
  indiv.design     = [];
	
  for i = 1:length(conditions)
    fprintf('analyzing condition %s against baseline.\n',num2str(conditions(i)));
    switch datafield
        %%%%%%% AVG STATS
        case 'epochs'
            ft_data = ts_data2fieldtrip(data_in,'condition',conditions(i),'channels',opt.channel,'chantype',opt.chantype,'dimord','chan_time');
%            ft_data = ts_data2fieldtrip(data_in,'condition',conditions(i),'channels',opt.channel,'dimord','chan_time');
            tla.keeptrials = 'yes';
            warning off all
            if strcmpi(indiv.method,'montecarlo')
              fprintf('turning off baseline correction for activation vs baseline\n');
              tla.blc = 'no'; 
            end
            tl_data    = timelockanalysis(tla,ft_data);
            if ~strcmp(class(tl_data.trial),opt.precision)
              fprintf('%s: converting trials to %s precision.\n',mfilename,opt.precision);
              eval(sprintf('tl_data.trial = %s(tl_data.trial);',opt.precision));
            end
            warning on all
            if isfield(tl_data,'dof'), tl_data = rmfield(tl_data,'dof'); end
            clear ft_data;

			  		% create design matrix
						if strcmpi(indiv.method,'montecarlo')
  		  			ntrials=data_in.(datafield)(conditions(i)).num_trials;
			  	  	design = zeros(2,2*ntrials);
			  	  	design(1,1:ntrials)=1;
			  	  	design(1,ntrials+1:2*ntrials)=2;
			  	  	design(2,1:ntrials)=[1:ntrials];
			  	  	design(2,ntrials+1:2*ntrials)=[1:ntrials];
		  		  	indiv.design = design;
		  		  	indiv.ivar  = 1;                         
					  	indiv.uvar  = 2;   						

              % split baseline and activation for permutation test
							lat = indiv.latency; % activation time interval
							t   = tl_data.time;							
							if opt.bl_tile
								lat = indiv.latency;
								if ischar(lat), lat=[0 max(t)]; end
								if lat(1)<0, lat(1)=0; end
								id1 = nearest(t,lat(1)) : min(nearest(t,lat(end)),nearest(t,max(t)));
								n1 = numel(id1);
                dat1 = rmfield(tl_data,{'time','trial','avg','var'});
								dat1.time = tl_data.time(1,id1);
								dat0 = dat1;
                dat1.trial = tl_data.trial(:,:,id1);
                dat1.avg   = tl_data.avg(:,id1);
                dat1.var   = tl_data.var(:,id1);
								if isa(opt.bl_window,'numeric') && length(opt.bl_window)==2
									bas = opt.bl_window;
								else
									bas = [min(t) min(t)/2];
								end
								id0 = nearest(t,bas(1)):nearest(t,bas(2));
								n0 = length(id0);
								if diff(lat) > diff(bas)
									% tile baseline
									ntiles = floor(n1/n0);
									nelems = rem(n1,n0);
									dat0.trial = repmat(tl_data.trial(:,:,id0),[1 1 ntiles]);
                  dat0.avg = repmat(tl_data.avg(:,id0),[1 ntiles]);
                  dat0.var = repmat(tl_data.var(:,id0),[1 ntiles]);
									if nelems
                    dat0.trial = cat(3,dat0.trial,tl_data.trial(:,:,1:nelems)); 
                    dat0.avg   = cat(2,dat0.avg,tl_data.avg(:,1:nelems));
                    dat0.var   = cat(2,dat0.var,tl_data.var(:,1:nelems));
                  end
								else
									n0 = n1;
									id0 = 1:n0;
									dat0.trial = tl_data.trial(:,:,id0);
                  dat0.avg   = tl_data.avg(:,id0);
                  dat0.var   = tl_data.var(:,id0);                  
								end
								ntimes = n1;
              else
                n0 = length(find(t < 0));
                n1 = length(find(t > 0));
                if isa(lat,'char') % use the entire baseline
                  di = min(n0,n1);	
                  lat = [0 t(nearest(t,0)+di)];
                elseif isa(lat,'numeric') 
                  % slice latency if necessary
                  lat(1) 	 = max(lat(1),min(t));
                  lat(end) = min(lat(end),max(t));
                  if lat(end)<0		
                  % all in baseline					
                    di = min(n1,nearest(t,lat(1)):nearest(t,lat(end)));
                  elseif lat(1)<0
                  % in both baseline and activation
                    di = min(n0,length(intersect(find(t>0),find(t<lat(end)))));
                  else
                  % all in activation
                    di = min(n0,nearest(t,lat(1)):nearest(t,lat(end)));
                  end
                end
                indiv.latency = lat;

                id0 				= 1:di;
                dat0 			 	= rmfield(tl_data,{'time','trial','avg','var'});
                dat0.time  	= tl_data.time(1,id0);
                dat0.trial 	= tl_data.trial(:,:,id0);
                dat0.avg 	 	= tl_data.avg(:,id0);
                dat0.var 	 	= tl_data.var(:,id0);

                id1 				= max(nearest(t,lat(1)),nearest(t,0)) + [id0];
                dat1 			 	= rmfield(tl_data,{'time','trial','avg','var'});
                dat1.time  	= tl_data.time(1,id1);
                dat1.trial 	= tl_data.trial(:,:,id1);
                dat1.avg 	 	= tl_data.avg(:,id1);
                dat1.var 	 	= tl_data.var(:,id1);

                dat0.time 	= dat0.time + (dat1.time(1) - dat0.time(1));
                ntimes = length(id0);
              end

% 							% activation vs baseline correction
% 							fprintf('Performing baseline correction (subtracting the time-averaged powspctrm mean)\n');
% 							idx = union(id0,id1);
% 							avg = repmat(nanmean(tl_data.trial(:,:,idx),3),[1 1 ntimes]);
% %							sd  = repmat(nanstd(tl_data.trial(:,:,:,idx),[],3),[1 1 ntimes]);
% 							dat0.trial = (dat0.trial - avg);% ./ sd;
% 							dat1.trial = (dat1.trial - avg);%./ sd;
% 							clear avg sd

              if isempty(opt.blcwindow), opt.blcwindow = [tl_data.time(1) 0]; end
              clear tl_data

              tl_data{1} = dat0; clear dat0;
							tl_data{2} = dat1; clear dat1;
							
            else
            	for j=1:data_in.(datafield)(conditions(i)).num_trials
                	indiv.design = cat(2,indiv.design,1);
									indiv.ivar = 1;
            	end		
							tl_data = {tl_data};
            end	
%            [ft_chidx jnk] = match_str(tl_data{c}.label,{data_in.sensor_info.label});
%            indiv.channel = ft_chidx;
            ft_stats{end+1}     = timelockstatistics (indiv,tl_data{:});

						if ~isempty(opt.blcwindow)
            	bl_samples = find(data_in.(datafield)(conditions(i)).time >= opt.blcwindow(1) & data_in.(datafield)(conditions(i)).time <= opt.blcwindow(2));
						else
							bl_samples = data_in.(datafield)(conditions(i)).time < 0;
						end
            for j=1:length(data_in.sensor_info)
                t = 0;
                if ~isnan(data_in.(datafield)(conditions(i)).data(j,:,:))
                    t = kstest(data_in.(datafield)(conditions(i)).data(j,bl_samples,:));
                end
                ft_stats{end}.norm_bl(j) = t;
            end
       
        %%%%%%%% TIME FREQ STATS
        case 'timefreq'
            %%% CHECK ME!
            tl_data        = ts_data2fieldtrip(data_in,'condition',conditions(i),'chantype',opt.chantype,'dimord','chan_freq_time');
            indiv.design   = 1;
            ft_stats{end+1}  = freqstatistics(indiv,tl_data);            
    end
%     ft_stats{i+1}.type = 'baseline vs activation';    
    num_trials(end+1)   = data_in.(datafield)(conditions(i)).num_trials;
    indiv.design = [];
    event_codes(end+1) = events(i);
    clear tl_data;
  end
end
if strcmpi(opt.type,'between') || strcmpi(opt.type,'both'), num_trials(1) = sum([num_trials]); end

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
%   stat_data.stats(i).type                 = ft_stats{i}.type;
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

%addpath /home/mmildev/matlab/localmods/images/images
