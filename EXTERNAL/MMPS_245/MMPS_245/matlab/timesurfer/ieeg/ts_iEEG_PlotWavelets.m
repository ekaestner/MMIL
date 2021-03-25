function ts_iEEG_PlotWavelets (varargin)

%  Use:  ts_iEEG_PlotWavelets ('f_name', f_name...)
%
%  Plots time-frequency representation of power or coherence.
%
%  Required inputs:
%
%    <None> - Will prompt for necessary file inputs.
%
%  Optional inputs: (also see FieldTrip documentation for multiplotTFR)
%
%    figpath      - path to save figures
%    f_path       - path to where wavelet files are located
%    f_name       - file names of wavelet files
%    condkey      - path/file name of key to condition numbers specifying name and
%                   event value that corresponds to each condition number
%    cls          - 'yes' or 'no' clear figure from screen - useful when
%                   processing a batch job
%    notes        - any additional strings to add to the plots
%
%  ts_iEEG_Plot Options that can also be specified:
%
%    xlim         - x axis limits {default = [min max]}
%    ylim         - y axis limits {default = [min max]}
%    zlim         - z axis limits or 'absmax' {default = [min max]}
%    baseline     - 'yes', 'no', 'absolute','relative','relchange','z-score' {default = 'no'}
%    baseline_start - default = first time point
%    baseline_end   - default = 0 
%    savefigs       - 'yes' or 'no' {default = 'no'}
%    channels       - cell array of channel names to be plotted
%    layoutfile     - a valid neuroscan montage or fieldtrip layout file
%    timedsfact     - downsampling factor - applies to time
%    freqdsfact     - downsampling factor - applied to freq
%
%  Created:        08/21/2007 by Rajan Patel
%  Last Modified:  10/04/2007 by Rajan Patel
%
%  See also ts_iEEG_FreqAnalysis, ts_iEEG_calc_wavelets, ts_iEEG_Plot.

%% Plot Settings

try
    options = varargin;
    for index = 1:length(options)
        if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; 
        end;
    end;
    if ~isempty( varargin ), opt=struct(options{:});
    else opt = []; end;
catch
    fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
    return;
end;

if isfield(opt,'badchanfile'), badchanfile = opt.badchanfile;
else                           badchanfile = [];             end
if ~exist(badchanfile,'file')
    fprintf('No bad channels selected...\n');
    badchans = [];
else
    fid=fopen(badchanfile,'rt');
    j = 1;
    while (~feof(fid))
     badchans{j} = fgetl(fid);
     j=j+1;
    end
    fclose(fid);
end

args = {};

if ~isempty(badchans),            args{end+1} = 'badchans';        args{end+1} = badchans; end
if isfield(opt,'xlim'),           args{end+1} = 'xlim';            args{end+1} = opt.xlim;  end
if isfield(opt,'ylim'),           args{end+1} = 'ylim';            args{end+1} = opt.ylim;  end
if isfield(opt,'zlim'),           args{end+1} = 'zlim';            args{end+1} = opt.zlim;  end
if isfield(opt,'showaxis'),       args{end+1} = 'showaxis';        args{end+1} = opt.showaxis; end
if isfield(opt,'cls'),            args{end+1} = 'clearfigs';       args{end+1} = opt.cls; end
if isfield(opt,'savefigs'),       args{end+1} = 'savefigs';        args{end+1} = opt.savefigs; 
    if isfield(opt,'figpath'),    args{end+1} = 'figpath';         args{end+1} = opt.figpath; end; end
if isfield(opt,'baseline'),       args{end+1} = 'baseline';        args{end+1} = opt.baseline; end
if isfield(opt,'baseline_start'), args{end+1} = 'baseline_start';  args{end+1} = opt.baseline_start; end
if isfield(opt,'baseline_end'),   args{end+1} = 'baseline_end';    args{end+1} = opt.baseline_end; end
if isfield(opt,'max_rows'),       args{end+1} = 'max_rows';        args{end+1} = opt.max_rows; end
if isfield(opt,'max_columns'),    args{end+1} = 'max_columns';     args{end+1} = opt.max_columns; end        
if isfield(opt,'paperorient'),    args{end+1} = 'paperorient';     args{end+1} = opt.paperorient; end  
if isfield(opt,'channels'),       args{end+1} = 'channels';        args{end+1} = opt.channels; end
if isfield(opt,'whichfigs'),      args{end+1} = 'whichfigs';       args{end+1} = opt.whichfigs; end
if isfield(opt,'showbadchans'),   args{end+1} = 'showbadchans';    args{end+1} = opt.showbadchans; end
if isfield(opt,'layoutfile'),     args{end+1} = 'layoutfile';      args{end+1} = opt.layoutfile; end
if isfield(opt,'laminar'),        args{end+1} = 'laminar';         args{end+1} = opt.laminar; end
if isfield(opt,'plotzscores'),    args{end+1} = 'plotzscores';     args{end+1} = opt.plotzscores; end
if isfield(opt,'freqcorr'),       args{end+1} = 'freqcorr';        args{end+1} = opt.freqcorr; end
if isfield(opt,'timedsfact'),     args{end+1} = 'timedsfact';      args{end+1} = opt.timedsfact; end
if isfield(opt,'freqdsfact'),     args{end+1} = 'freqdsfact';      args{end+1} = opt.freqdsfact; end

args{end+1} = 'showcolorbar';     args{end+1} = 'yes';

if isfield(opt,'condkey')
    condkey = opt.condkey;
    if exist(condkey,'file'), load(condkey);
    else fprintf('No key to condition numbers to load ...\n'); end
end;

if isfield(opt,'notes') && ~isempty(opt.notes)
  notes = opt.notes;
else
  notes = '';
end

if isfield(opt,'baseline_file') && ~isempty(opt.baseline_file)
  baseline_file = fullfile(opt.f_path,opt.baseline_file); 
  load(baseline_file);
  if ~exist('baseline_data','var')
    error('Baseline data file does not have baseline_data.'); 
  end
  args{end+1} = 'baseline_data';
  args{end+1} = baseline_data;
  if length(baseline_data.timefreq.num_trials) == 1
    notes = sprintf('Trials in baseline:\n%d\n%s',...
                   baseline_data.timefreq.num_trials,notes);
  else
    notes = sprintf('Trials in baseline:\nVaries\n%s',...
                    notes);    
  end
end

common_args = args;
clear args;

%% Set Condition Files

if isfield(opt,'f_path'), condpathnames = opt.f_path; end;
if isfield(opt,'f_name'), cond_files = opt.f_name;    
else
    [cond_files,condpathnames]=uigetfile('*.mat','Select Condition You Wish to Plot','MultiSelect','on');
end

%cd(condpathnames);

if ~iscell(cond_files), cond_files = {cond_files}; end

%% Plot

%[subject,date,task,condition,data_info,data_procs] = ts_iEEG_decon_fname (cond_files);

for j=1:length(cond_files)
    [subject,date,task,condition,data_info,data_procs] = ts_iEEG_decon_fname (cond_files{j});
    args = common_args;
    if ~exist(fullfile(condpathnames,cond_files{j}),'file')
      fprintf('ERROR: Could not find timefreq file: %s.\n',fullfile(condpathnames,cond_files{j}));
      return;
    else
      if strcmp(who('-file',fullfile(condpathnames,cond_files{j})),'timefreq_data'), load(fullfile(condpathnames,cond_files{j}));
      else fprintf('ERROR: Data in file must be timefreq_data: %s\n', fullfile(condpathnames,cond_files{j})); return; end
    end
    if ~isempty(condition), cond_num = str2num(strtok(condition{1},'cond')); end;
    if exist('cond_key','var') && ~isempty(find([cond_key.event_code] == timefreq_data.timefreq(1).event_code))
        cond_index = find([cond_key.event_code] == timefreq_data.timefreq(1).event_code);
        if    iscell(cond_key(cond_index).name), condition_name{j} = cond_key(cond_index).name{1};  
        else  condition_name{j} = cond_key(cond_index).name; end
        condition_name{j} = strrep(condition_name{j},'_','-');      % prevents _ making a character a subscript in the plot;
    else
        condition_name{j} = num2str(timefreq_data.timefreq(1).event_code);
    end    
    comp = {data_procs{1,:}};
    tag = '';
    if ~isempty(comp)
        for i=1:length(comp)
            if strfind(comp{i},'tfrejall')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('TF Rejected Across All\n%s',notes);
            elseif strfind(comp{i},'tfrejbychan')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('TF Rejected by Channel\n%s',notes);
            elseif strfind(comp{i},'tfrejgroups')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('TF Rejected by Chan Groups\n%s',notes);
            elseif strfind(comp{i},'blthresh')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('Baselines Rejected\n%s',notes);
            elseif strfind(comp{i},'ica')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('ICA Rejected\n%s',notes);
            elseif strfind(comp{i},'rej')
                tag = sprintf('%s_%s',comp{i},tag);
                notes = sprintf('Manual Rejections\n%s',notes);
            end
        end
    end
    args{end+1} = 'title';
    args{end+1} = sprintf('Subject: %s / Task: %s / Condition: %s',subject{1},task{1},condition_name{j});
    args{end+1} = 'figname';
    try
      args{end+1} = sprintf('%s_%s_%s_%swavelets.eps',subject{1},task{1},condition{1},tag);
    catch
      args{end+1} = sprintf('%s_%s_event%s_%swavelets.eps',subject{1},task{1},num2str(timefreq_data.timefreq(1).event_code),tag);
    end
    args{end+1} = 'footnotes';
    if length(timefreq_data.timefreq(1).num_trials) == 1
      args{end+1} = sprintf('Num Trials: %s\n%s',num2str(timefreq_data.timefreq(1).num_trials),notes);
    else
      args{end+1} = sprintf('Num Trials: Varies by Chan\n%s',notes);  
    end
    args{end+1} = 'conditionnames';
    args{end+1} = condition_name;
 
    args_list = args;
    args = {timefreq_data args_list{:}};
    ts_iEEG_Plot_TFR (args{:});
    clear args args_list;
end