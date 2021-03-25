function ts_iEEG_PlotAvg (varargin)

% Use: ts_iEEG_PlotAvg ('option',value ...);
%
% This will plot averages of iEEG data.  Data should be provided as a
% TimeSurfer avg_data structure saved in a matlab file.  Previous support
% for .eeg and .avg files has been removed.  Use appropriate data
% conversion utilities first.
%
% Use run_iEEG_PlotAvg to easily call this function.
%
% Visit neuroimaging.pbwiki.com to learn more about specified data
% structure as well as file naming convention.
%
% Require Input:
%
%    <none> - with no inputs specified in the function call the user will
%    be prompted to provide all necessary inputs during run time.
%
% Options:
%
%  figpath  - path to save files
%  f_path     - path to files
%  f_name     - file name of data file containing avg_data
%  condkey    - path/file name of key to condition numbers specifying name and
%               event value that corresponds to each condition number
%  conditions - the condition numbers to plot
%  eventvals  - event values to plot (specifying event values will override
%               condition numbers input)
%  cls        - 'yes' or 'no' - to close figures after plotting when doign a batch job
%               {default = 'no'}
%  badchanfile- list of bad channels - won't be placed on plot
%  stats      - 'yes' or 'no'
%  epochprefix - a string giving the prefix of the epoch files such that
%                the file name is epochprefix.cond#.mat for all the
%                conditions.
%  showbadchans- leave an empty space for bad channels {default = 'yes'}
%  alpha      - significance level (default = 0.05)
%  laminar    - 'yes' or 'no' to specify laminar data
%  layoutfile     - a valid neuroscan montage or fieldtrip layout file
%  notes      - a string containing any addition notes you wish to ass to
%                the legend
%
%  ts_iEEG_Plot Options that can also be specified:
%
%    xlim          =  x axis limits (default = [xmin xmax])
%    ylim          =  y axis limits (default = [ymin ymax])
%    showaxis      = 'yes' or 'no'  (default = 'no')
%    showzeroaxis  = 'yes' or 'no'  (default = 'yes')
%    savefigs      = 'yes' or 'no'  (default = 'no')
%    channels      = cell array of channels name to be plotted
%    statistics    = path and file_name of statistics file created by
%                    ts_iEEG_StatAnalysis that corresponds to the
%                    conditions/events being plotted.  The number of files
%                    listed should correspond to the number of plots being
%                    made.
%    alpha         = what is the significance value, (default = 0.05)
%
% Created      :  Rajan Patel 08/20/2007
% Last Modified:  Rajan Patel 09/20/2007
%
%
% See also ts_iEEG_Plot, run_iEEG_Process, ts_iEEG_StatAnalysis

%% Check Options

try
    options = varargin;
    for index = 1:length(options)
        if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
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

% Plotting Options

args = {};
if ~isempty(badchans),         args{end+1} = 'badchans';     args{end+1} = badchans; end
if isfield(opt,'xlim'),        args{end+1} = 'xlim';         args{end+1} = opt.xlim;  end
if isfield(opt,'ylim'),        args{end+1} = 'ylim';         args{end+1} = opt.ylim;  end
if isfield(opt,'showaxis'),    args{end+1} = 'showaxis';     args{end+1} = opt.showaxis; end
if isfield(opt,'showzeroaxis'),args{end+1} = 'showzeroaxis'; args{end+1} = opt.showzeroaxis; end
if isfield(opt,'cls'),         args{end+1} = 'clearfigs';    args{end+1} = opt.cls; end
if isfield(opt,'savefigs'),    args{end+1} = 'savefigs';     args{end+1} = opt.savefigs;
    if isfield(opt,'figpath'), args{end+1} = 'figpath';      args{end+1} = opt.figpath; end; end
if isfield(opt,'channels'),    args{end+1} = 'channels';     args{end+1} = opt.channels; 
else                                                         opt.channels = 'all';     end
if isfield(opt,'showbadchans'),args{end+1} = 'showbadchans'; args{end+1} = opt.showbadchans; end
if isfield(opt,'whichfigs'),   args{end+1} = 'whichfigs';    args{end+1} = opt.whichfigs; end;
if isfield(opt,'statsmarker'), args{end+1} = 'statsmarker'; args{end+1} = opt.statsmarker; end;
if isfield(opt,'markersize'),  args{end+1} = 'markersize'; args{end+1} = opt.markersize; end;
if isfield(opt,'layoutfile'),  args{end+1} = 'layoutfile'; args{end+1} = opt.layoutfile; end
if isfield(opt,'plot_colors'), args{end+1} = 'plot_colors'; args{end+1} = opt.plot_colors; end
if isfield(opt,'laminar'),     args{end+1} = 'laminar';     args{end+1} = opt.laminar; end
if isfield(opt,'conditions'),  plot_cond           = opt.conditions;
else                           plot_cond           = [];               end
if isfield(opt,'eventvals'),   plot_eventvals      = opt.eventvals;    end
if isfield(opt,'alpha') && ~isempty(opt.alpha),       alpha = opt.alpha;
else                                                  alpha = []; end
if isfield(opt,'baseline'),    args{end+1} = 'baseline';      args{end+1} = opt.baseline; end
if isfield(opt,'baseline_start'), args{end+1}='baseline_start'; args{end+1}=opt.baseline_start; end

if isfield(opt,'condkey')
    condkey = opt.condkey;
    if exist(condkey,'file'), load(condkey);
    else fprintf('No key to condition numbers to load ...\n'); end
end;

%% Select Average Files

avg_files = {};

if isfield(opt,'datafile')
  try
    for i = 1:length(opt.datafile)
      [opt.f_path fname fext] = fileparts(opt.datafile{i});
      avg_files{i} = [fname fext];
    end
  end
end

if isfield(opt,'f_name') && isempty(avg_files)
    avg_files    = opt.f_name;
end
if ~isempty(avg_files)
    if ~iscell(avg_files), avg_files={avg_files}; end;
    try
        avgpathnames = opt.f_path;
    catch
        fprintf('ERROR - Must specify path if file is specified');
        avg_files={};
    end
    if ~exist(fullfile(avgpathnames, avg_files{1}),'file')
        fprintf('ERROR - Data file not found: %s', fullfile(avgpathnames, avg_files{1}));
        avg_files={};
    end
end
if isempty(avg_files)
    [avg_files,avgpathnames,index]=uigetfile({'*.mat'},'Select data file','MultiSelect','off');
    avg_files={avg_files};
end


%% Load File

fprintf ('Loading average file...\n');
% if strcmp(who('-file',fullfile(avgpathnames,avg_files{1})),'avg_data')
if ismember('avg_data',who('-file',fullfile(avgpathnames,avg_files{1})))  
  load(fullfile(avgpathnames,avg_files{1}));
else 
  fprintf('Unexpected data in data file %s...\n', fullfile(avgpathnames,avg_files{1}));
  return;
end

try
    [subject, date, task, condition, data_info, data_procs] = ts_iEEG_decon_fname (avg_files);
catch
    clear legend_base subject date task condition data_info data_procs;
    subject={''}; date={''}; task={''}; condition={''};
 		data_info={''}; data_procs={''};
end
try	[legend_base] = ts_iEEG_const_legend (1,data_info, data_procs); catch legend_base={''}; end

if ~isstr(data_procs{1}), data_procs={''}; end;

% Determine what to plot.

if exist('plot_eventvals','var')
    plot_cond = [];
    for i=1:length(avg_data.averages)
        if ismember(avg_data.averages(i).event_code,plot_eventvals)
            plot_cond=cat(1,plot_cond,i);
        end
    end
    if length(plot_cond) ~= length(plot_eventvals)
        fprintf('ERROR: Invalid event values given.\n');
        return;
    end
end
if isempty(plot_cond)
    plot_cond = str2num(input('Enter condition numbers to plot:','s'));
end

try
  plotevents = [avg_data.averages(plot_cond).event_code];
catch
  error('ERROR: Invalid conditions.');
end

args{end+1} = 'conditions';
args{end+1} = plot_cond;
%plot_cond
%% Load Stats File

if isfield(opt,'statistics') && ~isempty(opt.statistics)
  if exist(opt.statistics,'file')
   if strcmp(who('-file',opt.statistics),'stat_data')
     fprintf('Loading statistics file: %s.\n',opt.statistics);
     load(opt.statistics);
   else
     error('Stats file contains invalid data: %s.',opt.statistics);
   end
   if ~isempty(alpha)
     for i = 1:length(stat_data.stats)          % modify the mask for the asked alpha
       stat_data.stats(i).mask = zeros(size(stat_data.stats(i).prob));
       stat_data.stats(i).mask(find(stat_data.stats(i).prob < alpha)) = 1;
       stat_data.stats(i).mask(find(stat_data.stats(i).prob >=  alpha)) = 0;
       stat_data.stats(i).parms.alpha = alpha;
     end
     % added by Thomas, July 4th, 2008; commented by Jason, July 21, 2008
     %stat_data.stats(1).mask=(stat_data.stats(1).posclusterslabelmat) | (stat_data.stats(1).negclusterslabelmat);
   end
%    if length(plot_cond) == length(intersect([stat_data.cond_num],plot_cond)) 
%      args{end+1} = 'stat_data';
%      args{end+1} = stat_data;
%    else
%      fprintf('The statistics file does not contain the correct conditions: %s.\n',opt.statistics);
%    end
   args{end+1} = 'stat_data';
   args{end+1} = stat_data;
  else
    fprintf('Could not find statistics file: %s.\n',opt.statistics);
  end
end


%% Set up titles, legends, ect.

for j=1:length(plotevents)
    if exist('cond_key','var')
        index = find([cond_key.event_code] == plotevents(j));     
        if    iscell(cond_key(index).name), condition_name{j} = cond_key(index).name{1};  
        else  condition_name{j} = cond_key(index).name; end
        condition_name{j} = strrep(condition_name{j},'_','-');      % prevents _ making a character a subscript in the plot
    else
        condition_name{j} = num2str(avg_data.averages(plot_cond(j)).event_code);
    end
    condition(j) = {num2str(avg_data.averages(plot_cond(j)).event_code)};  
end

foot_note = sprintf('Sampling Rate: %s\n',...
                   num2str(avg_data.sfreq));
                 

if exist('stat_data','var')
  foot_note = sprintf('%s%s\nSig: p < %s\n',foot_note,legend_base{1},num2str(stat_data.stats(1).parms.alpha));
else
  foot_note = sprintf('%s%s\n',foot_note,legend_base{1});
end

if isfield(opt,'notes') && ~isempty(opt.notes)
    foot_note = sprintf('%s%s',foot_note,opt.notes);
end


args{end+1} = 'conditionnames';
args{end+1} = condition_name;

args{end+1} = 'footnotes';
args{end+1} = foot_note;

if isfield(opt,'prefix')
  f_plot = [opt.prefix '_'];
else
  f_plot=[subject{1} '_' task{1} '_'];
end
for j=1:length(condition);
    f_plot = strcat(f_plot,'event',condition{j},'_');
end

tag = '';
if any(cell2mat(strfind(data_procs,'rej')))
 tag = 'rej_';
end

args{end+1} = 'figname';
args{end+1} = strcat(f_plot,tag,'avgs.eps');

args{end+1} = 'title';
if isfield(opt,'title')
  args{end+1} = opt.title;
else
  args{end+1} = sprintf('Average Data Plot - Subject: %s - Day %s Post Implant - Task: %s',subject{1},date{1},task{1});
end
if isfield(opt,'paperorient') && ~isempty(opt.paperorient)
  args{end+1} = 'paperorient';
  args{end+1} = opt.paperorient;
end

if isfield(opt,'max_rows')
  args{end+1} = 'max_rows';
  args{end+1} = opt.max_rows;
end
if isfield(opt,'max_columns')
  args{end+1} = 'max_columns';
  args{end+1} = opt.max_columns;
end
if isfield(opt,'resolution')
  args{end+1} = 'resolution';
  args{end+1} = opt.resolution;
end

%% Plot
args_list = args;
args = {avg_data args_list{:}};
ts_iEEG_Plot_Waveforms (args{:});
