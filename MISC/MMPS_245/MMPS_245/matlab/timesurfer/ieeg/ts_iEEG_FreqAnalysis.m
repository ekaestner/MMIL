function ts_iEEG_FreqAnalysis (varargin)

% Use ts_iEEG_FreqAnalysis ('epoch_dir',epoch_dir...);
%
% Uses epoch_data TimeSurfer data to perform frequency analysis using
% FieldTrip.
%
% Required input:
%
%    <None> - Will prompt user for necessary input if there is no default.
%
%    Note: If you sepcify a frequency analysis method there then maybe
%    required inputs, please see FieldTrip documentation for freqanalysis.
%
% Optional input:
%
%    f_path        - path to locate f_name
%    f_name        - name of epoch data file - epoch_data TimeSurfer format.
%    bl_name       - name of epoch data file for baseline events in f_path
%                    (if baseline events are diff than events to process)
%    freq_path     - name of directory to save time-frequency data in.
%    conditions    - conditions to analyze {default = all}
%    reject        - 'yes' or 'no' - time_frequency based artifact rejection {default = 'no'}
%    fa_method     - method option for freqanalysis: 'mtmfft',
%                    'mtmconvol', 'mtmwelch', 'wltconvol', 'tfr' {default = 'wltconvol'}
%
%      NOTE: You may specify any options available to the particular
%      frequanalysis method specified.  Please see
%      ts_freqanalysis_fieldtrip
%
% Created:  09/20/07 by Rajan Patel
% Rcnt Mod: 09/13/08 by Jason Sherfey
% Last Mod: 09/15/12 by Don Hagler
%
% See also ts_freqanalysis_fieldtrip, ts_freqanalysis_reject,
%          ts_freqanalysis_makebaseline ts_freqanalysis_bl_reject,
%          freqanalysis

%% Convert All Possible Args

opt          = mmil_args2parms( varargin, ...
               { 'f_path',[],[],...
                 'f_name',[],[],...
                 'bl_name',[],[],...
                 'freq_path','freq_analysis',[],...
                 'conditions',[],[],...
                 'event_codes',[],[],...
                 'combinations',[],[],...  %% FOR ADDING COMBINATIONS TO THIS SCRIPTS IN FUTURE
                 'neweventcodes',[],[],...
                 'combo_calc',[],[],...
                 'reject','no',{'yes','no'},...
								 'makebaseline','no',{'yes','no'},...
                 'fa_method','wltconvol',{'mtmfft','mtmconvol','mtmwelch','wltconvol','tfr'},...
                 'foi',           [1 2 3 4:1:60 62:2:200],    [],...
                 'toi',           [], [],...
                 'channels',      [],          [],...
                 'chantype',      [], {'mag','grad','grad1','grad2','eeg','other','meg','all'},...
                 'channelcmb',    {'all' 'all'},[],...
                 'keeptrials',    'no',        {'yes','no'},...
                 'savefieldtrip', false,       sort([false true]),...
                 'output' , 'pow', {'pow','powandcsd','fourier'},...
                 'width', [],[],...%[1 1 1 linspace(1,15,length(4:1:60)) (62:2:200)/4], [],...
                 'gwidth', [],[],...%[3], [],...
                 'taper', 'dpss', {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin',...
                                   'chebwin','flattopwin','gausswin','hamming','hann','kaiser','nuttallwin',...
                                   'parzenwin','rectwin','tukeywin','triang','dpss'},...
                 'tapsmofrq', [4], [],...
                 'keeptapers','no',{'yes','no'},...
                 'pad','maxperlen',[] ,...
                 'downsample',1, [],...
                 't_ftiwin',[],[],...
                 'reject',   'yes',  {'yes','no'},...
                 'reject_exclude','bychannel',{'inall','bychannel','groupchannels'},...
                 'baseline_start',-.08,[],...
                 'baseline_end',-.005,[],...                      
                 'bl_threshold',10,[],...
                 'trial_threshold',40,[],...
                 'chan_threshold',80,[],...
                 'baseline_events',[],[],...
                 'baseline_file',[],[],...
                 'verbose',      true,     sort([false true]), ...
                 'logfile',      [],       [],...
                 'logfid',       [1],        [], ...
								 'visualreject_flag', false, sort([false true])...
                }, ...
                false );

%% Check Options              
              
if isempty(opt.f_name) % Prompt user for epoch files
  [opt.f_path,opt.f_name] = uigetfile({'*.mat'},'Select epoch_data file(s) to process...','MultiSelect','on');
  if ~iscell(opt.fname), opt.fname = {opt.fname}; end
end

% Distribute path to all file names
if ~iscell(opt.f_name), opt.f_name = {opt.f_name}; end
for i = 1:length(opt.f_name)
  opt.f_name{i} = fullfile (opt.f_path,opt.f_name{i});
end

% Check input data
for i = 1:length(opt.f_name)
  if ~exist(opt.f_name{i},'file')
    mmil_error(opt,'Cannot find file: %s',opt.f_name{i})
  elseif ~strcmp('epoch_data',who(opt.f_name{i},'file'));
    mmil_error(opt,'File contains invalid data (all data must be epoch_data): %s',opt.f_name{i});    
  end
end

% Distribute path to all file names
if ~isempty(opt.bl_name) && ~iscell(opt.bl_name), opt.bl_name = {opt.bl_name}; end
for i = 1:length(opt.bl_name)
  opt.bl_name{i} = fullfile (opt.f_path,opt.bl_name{i});
end

% Check input data
for i = 1:length(opt.bl_name)
  if ~exist(opt.bl_name{i},'file')
    mmil_error(opt,'Cannot find file: %s',opt.bl_name{i})
  elseif ~strcmp('epoch_data',who(opt.bl_name{i},'file'));
    mmil_error(opt,'File contains invalid data (all data must be epoch_data): %s',opt.bl_name{i});    
  end
end

% Check for baseline_file
if ~isempty(opt.baseline_file)
    if ~exist(opt.baseline_file,'file')
        mmil_error(opt,'Cannot find baseline data file: %s',opt.baseline_file);
    elseif ~strcmp('baseline_data',who(opt.baseline_file,'file'))
        mmil_error(opt,'Baseline file does not contain the appropriate data): %s.',opt.baseline_file);
    else        
        load(opt.baseline_file);
        opt.baseline_data = baseline_data;
				opt = rmfield(opt,'baseline_file');
    end
end

% Check output directory exists
%opt.freq_path = fullfile(opt.f_path,opt.freq_path);
if ~exist(opt.freq_path,'dir')
  mmil_logstr(opt,'Creating directory: %s.',opt.freq_path);
  mkdir(opt.freq_path); 
end

opt.method = opt.fa_method;

% Turn on or off reject
if strcmpi(opt.reject,'yes')
  opt.reject_flag = true;
else
  opt.reject_flag = false;
end
opt = rmfield(opt,'reject');

% Turn on or off saving trial data
if strcmpi(opt.keeptrials,'yes')
	opt.trials_flag = true;
else
	opt.trials_flag = false;
end
opt = rmfield(opt,'keeptrials');

% Turn on or off making baseline_data
if strcmpi(opt.makebaseline,'yes')
	opt.makebaseline_flag = true;
else
	opt.makebaseline_flag = false;
end
opt = rmfield(opt,'makebaseline');
% conditionally override makebaseline option
if isfield(opt,'baseline_data')
	opt.makebaseline_flag = false;
elseif opt.reject_flag
	opt.makebaseline_flag = true;
end

% Convert seconds to ms for TimeSurfer function
opt.toi      = opt.toi * 1000;       
opt.t_ftiwin = opt.t_ftiwin * 1000;
opt.baseline_start = opt.baseline_start * 1000;
opt.baseline_end = opt.baseline_end * 1000;


%% Perform wavelet analysis

%args = mmil_parms2args (rmfield(opt,{'f_path','f_name','bl_name','freq_path','fa_method'}));
args = mmil_parms2args (rmfield(opt,{'f_path','bl_name','fa_method'}));

for i = 1:length(opt.f_name)
    mmil_logstr(opt,'Loading data file: %s',opt.f_name{i});
    load(opt.f_name{i});
    epoch_data_array(i) = epoch_data;
end

for i = 1:length(opt.bl_name)
    mmil_logstr(opt,'Loading data file: %s',opt.bl_name{i});
    load(opt.bl_name{i});
    epoch_data_array(end+1) = epoch_data;
end

if length(epoch_data_array) > 1
    epoch_data = ts_combine_data(epoch_data_array);
else
    epoch_data = epoch_data_array;
end

clear epoch_data_array

if ~isempty(opt.baseline_events)
   blconds = find(ismember([epoch_data.epochs.event_code],[opt.baseline_events]));
else
   blconds = 1:length(epoch_data.epochs);
end

if opt.makebaseline_flag
	[timefreq_data,badtrials,baseline_data,badbaselines] = ts_freqanalysis_fieldtrip(epoch_data,args{:});
else
	[timefreq_data,badtrials] = ts_freqanalysis_fieldtrip(epoch_data,args{:});
end
%[timefreq_data,badtrials,badbaselines] = ts_freqanalysis_fieldtrip(epoch_data,args{:});

 %% TO DO:
 %% ADD COMBINATIONS HERE: USE TS_COMBINE_CONDITIONS and then get rid or
 %% ts_iEEG_combos all together.  use combinations options already set up

%% Setup file names
% Determine file name to save as
[dum, file_tag, ext] = fileparts(opt.f_name{1});
clear dum ext;
file_tag = strrep(file_tag,'.epoch','');
if strfind(file_tag,'event')
    file_name = file_tag(1:strfind(file_tag,'event')-1);
elseif strfind(file_tag,'cond')
    file_name = file_tag(1:strfind(file_tag,'cond')-1);
elseif strfind(file_tag,'epoch_data')
		% filename generated by ts_process_fif_data
		file_name = [file_tag(1:strfind(file_tag,'epoch_data')+9) '.'];		
end
%file_name = strcat(file_name,...
%    'foi',strrep(num2str(opt.foi(1)),'.',''),'-',strrep(num2str(opt.foi(end)),'.',''),...
%    '.toi',strrep(num2str(opt.toi(1)/1000),'.',''),'-',strrep(num2str(opt.toi(end)/1000),'.',''));

str_foi = strrep(sprintf('foi%d-%d',opt.foi(1),opt.foi(end)),'.','');
str_toi = strrep(sprintf('toi%.4g-%.4g',opt.toi(1)/1000,opt.toi(end)/1000),'.','');
file_name = sprintf('%s%s.%s',file_name,str_foi,str_toi);

if opt.reject_flag
  switch opt.reject_exclude
       case 'inall'
           file_name = sprintf('%s.tfrejall%s',file_name,num2str(opt.trial_threshold));
       case 'bychannel'
           file_name = sprintf('%s.tfrejbychan%s',file_name,num2str(opt.trial_threshold));
       case 'groupchannels'
           file_name = sprintf('%s.tfrejgroups%s',file_name,num2str(opt.trial_threshold));
  end
  file_name = sprintf('%s.blthresh%s',file_name,num2str(opt.bl_threshold));
	file_name = sprintf('%s.chanthresh%s',file_name,num2str(opt.chan_threshold));		

	try
	% make tf_reject_data
		for cond=1:length(timefreq_data.timefreq)
			nchan = length(badtrials{cond});
			chans = 1:nchan;
			if ~isempty(opt.channels) && isnumeric(opt.channels) && max(opt.channels)<=length(badtrials{cond})
				nchan = length(opt.channels);
				chans = opt.channels;
			end	
			tf_reject_data(cond).event_code = timefreq_data.timefreq(cond).event_code;
			tf_reject_data(cond).original_num_trials = epoch_data.epochs(cond).num_trials;			
			for ch=1:nchan
				tf_reject_data(cond).sensor_info(ch) = timefreq_data.sensor_info(chans(ch));
				tf_reject_data(cond).badtrials{ch} = badtrials{cond}{chans(ch)};
			end
			badchan_idx = find([timefreq_data.sensor_info.badchan]);
			tf_reject_data(cond).badchans = {timefreq_data.sensor_info(badchan_idx).label};
		end
	%% save timefreq reject data
	reject_name = strcat(file_name,'.tf_reject_data.mat');
	reject_name = fullfile(opt.freq_path,reject_name);
	save(reject_name,'tf_reject_data');
	end
end


%% Save Baseline Stuff
if opt.makebaseline_flag
	bl_name   = strcat(file_name,'.baseline.mat');
	bl_name   = fullfile(opt.freq_path,bl_name);
	save(bl_name,'baseline_data');
end

if opt.reject_flag && exist('badbaselines','var')
    bl_reject = strcat(file_name,'.baseline.reject_info.txt');
    bl_reject = fullfile(opt.freq_path,bl_reject);
    mmil_logstr(opt,'Saving Bad Trial Info: %s',bl_reject);
    fid = fopen(bl_reject,'w');
    fprintf(fid,'Rejected Trials for Baseline Data.\n\n');
    fprintf(fid,'Original Number of Trials: %s\n',num2str(sum([epoch_data.epochs(blconds).num_trials])));
    fprintf(fid,'\nChannel\n\n');
    for b = 1:length(badbaselines)
        fprintf(fid,'%-7sNum Bad Trials: %s\n',baseline_data.sensor_info(b).label,num2str(length(badbaselines{b})));
    end
    fprintf(fid,'\nTotal Bad Trials across all channels: %s\n',num2str(length(unique(cell2mat(badbaselines)))));  
    fprintf(fid,'Percentage of Trials across all channels: %s%%\n\n',num2str((length(unique(cell2mat(badbaselines)))/sum([epoch_data.epochs.num_trials]))*100));  
    fprintf(fid,'Trials that were removed:\n');
    fprintf(fid,'Channel\n\n');
    for b = 1:length(badbaselines)
        fprintf(fid,'%-7sBad Trials: %s\n',baseline_data.sensor_info(b).label,num2str(badbaselines{b}));
    end
    fclose(fid);
end


%% Save timefreq data event by event

tf_name = strcat(file_name,'.timefreq.');
timefreq_conditions = timefreq_data.timefreq;
timefreq_data.timefreq = [];

for i = 1:length(timefreq_conditions)
    timefreq_data.timefreq = timefreq_conditions(i);

		timefreq_data.timefreq.power = single(timefreq_data.timefreq.power);
		timefreq_data.timefreq.data = single(timefreq_data.timefreq.data);		

    if strfind(file_tag,'event')
        tf_file_name = strcat(tf_name,'event',num2str(timefreq_data.timefreq.event_code));
    else
        tf_file_name = strcat(tf_name,'cond',num2str(i));
    end
    tf_reject_name = sprintf('%s.reject_info.txt',tf_file_name);
    tf_reject_name = fullfile(opt.freq_path,tf_reject_name);
    tf_file_name = sprintf('%s.mat',tf_file_name);
    tf_file_name = fullfile(opt.freq_path,tf_file_name);
    % Save File
    mmil_logstr(opt,'Saving timefreq_data: %s.',tf_file_name);
    save(tf_file_name,'timefreq_data');
    if opt.reject_flag
        mmil_logstr(opt,'Saving Bad Trial Info: %s',tf_reject_name);
        fid = fopen(tf_reject_name,'w');
        fprintf(fid,'Rejected Trials for event %s.\n\n',num2str(timefreq_data.timefreq.event_code));
        fprintf(fid,'Original Number of Trials: %s\n\n',num2str(epoch_data.epochs(i).num_trials));
        fprintf(fid,'Channel\n\n');
        for b = 1:length(badtrials{i})
            fprintf(fid,'%-7sNum Bad Trials: %s\n',timefreq_data.sensor_info(b).label,num2str(length(badtrials{i}{b})));
        end
        fprintf(fid,'\nTotal Bad Trials across all channels: %s\n',num2str(length(unique(cell2mat(badtrials{i})))));
        fprintf(fid,'Percentage of Trials across all channels: %s%%\n\n',num2str((length(unique(cell2mat(badtrials{i})))/epoch_data.epochs(i).num_trials)*100));
        fprintf(fid,'Trials that were removed:\n');
        fprintf(fid,'Channel\n\n');
        for b = 1:length(badtrials{i})
            fprintf(fid,'%-7sBad Trials: %s\n',timefreq_data.sensor_info(b).label,num2str(badtrials{i}{b}));
        end
        fclose(fid);
    end
end
