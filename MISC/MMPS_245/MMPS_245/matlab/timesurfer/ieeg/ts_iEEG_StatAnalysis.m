function ts_stat_analysis (varargin)

% Use ts_stat_analysis ('option',value,...
% Note: this function replaced ts_iEEG_StatAnalysis
%
% Perform statistical analysis using FielTrip
%
% Required input:
%
%   f_path  - path to mat files 
%   stat_path - path to where stat_data will be saved                              
%   f_name  - list of epoch_data files for the study%
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
%  See also: ts_iEEG_statistics, timelockstatistics, freqstatistics
%
%  Created on ?? by Rajan Patel
%  Last Modified 13-Sep-2008 by Jason Sherfey

%% Check Options

if ~mmil_check_nargs(nargin, 4), return; end;

opt = mmil_args2parms(varargin,...
                      {'f_path'     , [],[],...
											 'stat_path'  , [],[],...
                       'f_name'    , [],[],...
                       'conditions' , [],[],...
                       'events'     , [],[],...
                       'method'     , [],{'montecarlo','analytic','stats','glm'},...
                       'statistic'  , [],{'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                                          'indepsamplesZcoh','depsamplesT','depsamplesF' ,...
                                          'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
                                          'anova1','kruskalwallis'},...
											 'correctm'   , [],{'no','max','cluster','bonferoni','fdr'},...
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
											 'frequency', 'all',[],...
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
											 'blcorrection','absolute',{'no','absolute','relative','relchange','zscore','zscores','yes',''},...
											 'bl_window',[],[],...
											 'test_type','both',{'within_trials','between_trials','both'},...
                       'detrend','no',{'yes','no'},...
                      },...
                      false);

%                      'channel'    , 'all',[],...
%											 'chantype',      [], {'mag','grad','grad1','grad2','eeg','other','meg','all'},...

%% Get input filenames

if ~iscell (opt.f_name)
  f_names = {opt.f_name};
else
  f_names = opt.f_name;
end

f_path = opt.f_path;

timefreq_flag = 0;
[path,stats_file] = fileparts(f_names{1});
if strfind(stats_file,'timefreq')
  stats_file = stats_file(1:strfind(stats_file,'timefreq')+8);
	timefreq_flag = 1;	
elseif strfind(stats_file,'epoch')
  stats_file = stats_file(1:strfind(stats_file,'epoch')-1);
end
clear path

if isfield(opt,'method') && ~isempty(opt.method)
  stats_file = [stats_file opt.method '.'];
else
  error('There was no ''method'' specified.');
end

if isfield(opt,'statistic') && ~isempty(opt.method)
  stats_file = [stats_file opt.statistic '.'];
else
  error('There was no ''statistic'' specified.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if timefreq_flag && strcmpi(opt.method,'montecarlo')
	% get timesurfer info from files with timefreq_data averages
	for i=1:length(f_names)
		tmp_str = fullfile(f_path,f_names{i});
		if ~exist(tmp_str,'file'),fprintf('ERROR: File not found - %s\n',tmp_str); return; end;
		load(tmp_str); % timefreq_data
		data(i) = timefreq_data;
		data(i).timefreq =rmfield(data(i).timefreq,{'power','data'});
		clear timefreq_data;
	end
	
	% Check channel selections
  if (~isempty(opt.chantype) && ~isempty(opt.channel))
    mmil_error(parms, 'Cannot specify BOTH chantype AND channels param values.');
  % Choose all channels by default
  elseif (isempty(opt.chantype) && isempty(opt.channel))
    opt.channel = 1:data(1).num_sensors;
  % Convert all inputs directly to channel indexes
  elseif (isempty(opt.channel))
    switch opt.chantype
     case {'mag' 'grad1' 'grad2' 'eeg', 'other'}
       opt.channel = find(strcmp(opt.chantype,{data(1).sensor_info.typestring}));
     case {'grad'}
       opt.channel = find(strncmp(opt.chantype,{data(1).sensor_info.typestring},...
         length(opt.chantype)));
     case 'meg'
       [a,opt.channel] = find(ismember({data(1).sensor_info.typestring}, ...
         {'mag', 'grad1', 'grad2'}));
     case 'all'
       opt.channel = setdiff(1:data(1).num_sensors, find(strcmp('other', {data(1).sensor_info.typestring})));
    end;
  end;
	channels = opt.channel;

	if isfield(opt,'events') && ~isempty(opt.events)
		tmp_str = 'event'; 
		events = opt.events;
	elseif isfield(opt,'condition') && ~isempty(opt.conditions)
		tmp_str = 'condition'; 
		events = opt.conditions; 
		conditions = opt.conditions;
	else error('Must supply either condition numbers or event codes in order for proper file names to be created.');
	end
	opt.conditions = 1:length(data);
	opt.events     = [];
	f_path = sprintf('%s/%s',f_path,tmp_str);
	
	% compute each channel's stats individually b/c of memory limitations
	for ch = 1:length(channels)
		% load TF trial data for one channel
	  for k=1:length(events),
			infile = sprintf('%s.channel_%03i.mat',f_names{k}(1:strfind(f_names{k},'.mat')-1),channels(ch));%data(1).sensor_info(channels(ch)).lognum);
	    infile = sprintf('%s%d/%s',f_path,events(k),infile);
	    load(infile); 
			data_in(k) = timefreq_data;
	    clear timefreq_data;
		end		
		if length(data_in) > 1 
			data_in = ts_combine_data(data_in); 
		end;
		opt.channel = 1; % since there is only 1 channel
		args = mmil_parms2args(opt);

		stat = ts_timefreq_statistics(data_in,args{:});
	
    % concatenate stat structures
    if ch==1
        stat_full=stat;
    else
			for i=1:length(stat.stats)
	      stat_full.stats(i).prob=cat(1,stat_full.stats(i).prob,stat.stats(i).prob);
	      stat_full.stats(i).mask=cat(1,stat_full.stats(i).mask,stat.stats(i).mask);
	      stat_full.stats(i).stat=cat(1,stat_full.stats(i).stat,stat.stats(i).stat);
	      try
					stat_full.stats(i).posclusters	=	cat(2,stat_full.stats(i).posclusters,stat.stats(i).posclusters);
		      for i=1:length(stat.stats(i).posclusters) % shift cluster numbers so they match up
		        stat.stats(i).posclusterslabelmat(stat.stats(i).posclusterslabelmat==i) = i+length(stat_full.stats(i).posclusters);
		      end
					stat_full.stats(i).posclusterslabelmat=cat(1,stat_full.stats(i).posclusterslabelmat,stat.stats(i).posclusterslabelmat);
					stat_full.stats(i).posdistribution=cat(1,stat_full.stats(i).posdistribution,stat.stats(i).posdistribution);
				end			
				try
					stat_full.stats(i).negclusters=cat(2,stat_full.stats(i).posclusters,stat.stats(i).negclusters);
		      for i=1:length(stat.stats(i).negclusters)
	          stat.stats(i).negclusterslabelmat(stat.stats(i).negclusterslabelmat==i)=i+length(stat_full.stats(i).negclusters);
	  	    end
					stat_full.stats(i).negclusterslabelmat=cat(1,stat_full.stats(i).negclusterslabelmat,stat.stats(i).negclusterslabelmat);
					stat_full.stats(i).negdistribution=cat(1,stat_full.stats(i).negdistribution,stat.stats(i).negdistribution);
				end
				try stat_full.stats(i).label=cat(1,stat_full.stats(i).label,stat.stats(i).label); end
  		end
	  end		
		stat_full.sensor_info = cat(1,stat_full.sensor_info,stat.sensor_info);	
	end
	
	[sens,idx] = unique({stat_full.sensor_info.label});
	stat_full.sensor_info = stat_full.sensor_info(idx);
	stat_full.num_sensors = length(stat_full.sensor_info);
	
	%% Sort Clusters
	for i=1:length(stat_full.stats)
		% sort positive clusters by p-value
	 try
		[sortp,sort_ind]=sort([stat_full.stats(i).posclusters.prob],'ascend');
		stat_full.stats(i).posclusters=stat_full.stats(i).posclusters(sort_ind);
		labelmattemp=stat_full.stats(i).posclusterslabelmat;
		for i=1:length(stat_full.stats(i).posclusters)
	    labelmattemp(stat_full.stats(i).posclusterslabelmat==sort_ind(i))=i;
		end
		stat_full.stats(i).posclusterslabelmat=labelmattemp;
	 end
	 try
		% sort negative clusters by p-value
		[sortp,sort_ind]=sort([stat_full.stats(i).negclusters.prob],'ascend');
		stat_full.stats(i).negclusters=stat_full.stats(i).negclusters(sort_ind);
		labelmattemp=stat_full.stats(i).negclusterslabelmat;
		for i=1:length(stat_full.stats(i).negclusters)
	    labelmattemp(stat_full.stats(i).negclusterslabelmat==sort_ind(i))=i;
		end
		stat_full.stats(i).negclusterslabelmat=labelmattemp;
	 end
	end

stat_data=stat_full;

else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(f_names)
  f_names{i} = fullfile (f_path,f_names{i});
  if ~exist(f_names{i},'file')
    fprintf('ERROR: File not found - %s\n',f_names{i});
    return;
  else
    if strcmp(who('-file',f_names{i}),'timefreq_data')
      fprintf('Loading timefreq data file for analysis: %s\n',f_names{i});
      load(f_names{i});
      data(i) = timefreq_data;
      clear timefreq_data 
    elseif strcmp(who('-file',f_names{i}),'epoch_data')
      fprintf('Loading epoch data file for analysis: %s\n',f_names{i});
      load(f_names{i});
      data(i) = epoch_data;
      clear epoch_data
    else
       fprintf('ERROR: Files must contain epoch_data - %s\n',f_names{i}); 
    end
  end
end

if isfield (opt,'events') && ~isempty(opt.events)
  events     = opt.events; 
elseif isfield(opt,'conditions') && ~isempty(opt.conditions)
  conditions = opt.conditions;
else
  error('Must supply either condition numbers or event codes in order for proper file names to be created.');
end

opt.conditions = 1:length(data);
opt.events     = [];

if length(data) > 1
  data_in = ts_combine_data(data);
else
  data_in = data;
end

%% Perform Analysis

opt = rmfield (opt,'f_path');
opt = rmfield (opt,'f_name');

args = mmil_parms2args(opt);
stat_data = ts_statistics(data_in,args{:});

a=1:10;
%pd=pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%save stat_data stat_data
%% Save Stat Structure

stats_path = opt.stat_path; %fullfile(f_path,'stats');

if ~exist(stats_path,'dir')
  fprintf('Creating: %s\n',stats_path);
  mkdir(stats_path);
end

stats_file = [stats_file 'stat.'];

if exist('conditions','var')
  stats_file = [stats_file 'cond'];
  for i = 1:length(conditions)
    stats_file = [stats_file num2str(conditions(i)) '.'];
  end
elseif exist ('events','var')
  stats_file = [stats_file 'event'];
  for i = 1:length(events)
    stats_file = [stats_file num2str(events(i)) '.'];
  end
end
stats_file = [stats_file 'mat'];

stats_file = fullfile (stats_path,stats_file);

if exist (stats_file,'file')
  fprintf('WARNING: Overwriting file - %s\n');
end

fprintf('Saving stats: %s\n',stats_file);
save (stats_file,'stat_data');
