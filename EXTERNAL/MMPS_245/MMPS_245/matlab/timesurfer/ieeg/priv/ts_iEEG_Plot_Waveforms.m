function ts_iEEG_PlotWaveforms (data,varargin)

% Use:  ts_iEEG_PlotWaveforms (data,'conditions',[1 4 5],...)
%
% Note: If you plot lots of conditions together, MATLAB can get bogged down 
%       considerably.  In addition more than 16 may cause anomalies with
%       the legends or cause the program to crash.
%
% Required Input:
%    
%   data - a TimeSurfer avg_data structure.
%
% Optional Parameters: 
%
%    conditions     - the conditions (not event numbers) to plot  
%                     {default is the smaller of: all conditions or first 4 conditions}
%    channels       - cell array of channel names to plot.  {default = all}
%    badchans       - list of channels that should not be plotted {default
%    = []} 
%    baseline       - 'yes','relative','relchange','absolute', 'z-score', 'no' 
%                      {default = 'no'} {'yes' = 'absolute'}
%                     'relative' and 'relchange' only valid for time-frequency data
%    baseline_start - start time for calculation of baseline {default = first time point}
%    baseline_end   - end time for calculation of baseline   {default = 0}
%    xlim           - [xmin xmax] - x-axis limits {default based on max and min of time}
%    ylim           - [ymin ymax] - y-axis limits {default based on max and min of data}
%    zlim           - [zmin zmax] - z-axis limits {default based on max and min of data}
%                     or 'absmax'
%    showaxis       - 'yes' or 'no' - show MATLAB axis {default = 'yes'}
%    showzeroaxis   - 'yes' or 'no' - show axis where x and y = zero {default = 'yes'} 
%    conditionnames - cell array of strings describing each condition that
%                     is plotted {default = condition numbers}
%    footnotes      - string describing data {default = ''}
%    figname        - template file name for the figure
%    title          - title placed in figure
%    layoutfile     - a valid neuroscan montage or fieldtrip layout file if
%                     no layout file is specified the plots are auto
%                     generated at run time
%    savefigs       - 'yes' or 'no' - save figure to file {default = 'no'}
%    clearfigs      - 'yes' or 'no' - clear figures off screen {default ='no'}
%    stat_data       - a valid stat_data data structure
%    laminar        - 'yes' or 'no'
%
% Advanced options:
%
%    max_rows, max_columns, max_group_size, plot_colors, figbackground,
%    zeroaxiswidth, zeroaxisstyle, axisfontsize, papertype, leftmargin,
%    rightmargin, topmargin, bottommargin, paperorient, shading, whichfigs,
%    rowpad, columnpad.
%
% Created:       09/25/2007 - Rajan Patel
% Last Modified: 01/11/2008 - Rajan Patel
%
% See also: ts_iEEG_PlotAvg, ts_iEEG_PlotWavelets

%% Check for errors in calling function

if nargin < 1
   help(mfilename);
   return;
end

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), opt=struct(options{:}); 
  else opt= []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

%% Data Options
                                                   
if isfield(data,'averages')                                                % what type of data structure
   datafield = 'averages';  
else
  fprintf('ERROR: Supplied data is not supported by this function.\n');
  return;
end

if isfield(opt,'conditions')
  conds = opt.conditions;
  for i = 1:length(conds)
   if conds(i) > length(data.(datafield))
    fprintf('ERROR: Condition %s is outside of the range of data supplied.\n',num2str(conds(i)));
    return
   end
  end
   if length(conds)>16
    fprintf('WARNING: More than 16 conditions may cause problems.\n');
   end
else
  num_conditions  = length(data.(datafield));
  num_condstoplot = min(4, num_conditions); 
  conds           = 1:num_condstoplot;                                     % default plot all or first 16 conditions
end

if isfield(opt, 'channels') && ~isempty(opt.channels)
    if ~iscell(opt.channels), opt.channels = {opt.channels}; end
    channels = opt.channels;  
    channels = intersect({data.sensor_info.label},channels);               % find the channels that match those in the data
    data = ts_data_selection(data,'chanlabel',channels);
    if isempty(channels)
        fprintf('ERROR: %s not found as channels in given data.\n',channels{:});
        return;
    elseif length(channels) ~= length(opt.channels)
        fprintf('WARNING: One or more channels not found in supplied data:\n%s\n',cell2mat(setdiff(opt.channels,channels)));
    end
else
    channels = {};
end

if isfield(opt, 'badchans')                                                % channels not to plot
  ignore_channels = opt.badchans;
else
  ignore_channels = [];
end

if isfield(opt,'showbadchans') && strcmpi(opt.showbadchans,'no')
 badchan_idx     = find([data.sensor_info.badchan]);                        % don't plot channels marked bad
 for i = 1:length(ignore_channels)
     for s = 1: length(data.sensor_info)
         if strcmp(data.sensor_info(s).label,ignore_channels{i})
             data.sensor_info(s).badchan = 1;
         end
     end
 end
 ignore_channels = [{data.sensor_info(badchan_idx).label}];
end

baseline = 0;
if isfield(opt,'baseline')
  if strcmpi(opt.baseline,'yes') || ...
     strcmpi(opt.baseline,'absolute'), baseline = 1; end
  if strcmpi(opt.baseline,'relative'), baseline = 2; end
  if strcmpi(opt.baseline,'relchange'),baseline = 3; end
  if strcmpi(opt.baseline,'z-score'),  baseline = 4; end
end

if isfield(opt,'baseline_start')
  baseline_start = opt.baseline_start;
else
  baseline_start = [];
end

if isfield(opt,'baseline_end')
 baseline_end = opt.baseline_end;
else
 baseline_end = [];
end

plot_stats = 0;
if isfield(opt,'stat_data') && ~isempty(opt.stat_data)
  plot_stats = 1;
end

%% Figure Options 

if isfield(opt,'xlim') && ~isempty(opt.xlim)                                                    
  xmin = opt.xlim(1);
  xmax = opt.xlim(2);
else
  xmin = [];
  xmax = [];
  for i=1:length(conds)                                                     
    xmin = min([xmin data.(datafield)(conds(i)).time]);                     % look for x min across all selected conditions
    xmax = max([xmax data.(datafield)(conds(i)).time]);                     % look for x max across all selected conditions
  end
end

if isfield(opt,'ylim') && ~isempty(opt.ylim)
  ymin = opt.ylim(1);
  ymax = opt.ylim(2);
  calc_ylim = 0;
else
  calc_ylim = 1;
  ymin = [];
  ymax = [];
end

abs_zlim=0;
calc_zlim = 1;
zmin = [];
zmax = [];
if isfield(opt,'zlim')
 if strcmpi(opt.zlim,'absmax')
  abs_zlim=1; 
 elseif ~isempty(opt.zlim)
  zmin = opt.zlim(1);
  zmax = opt.zlim(2);
  calc_zlim = 0;
 end
end

show_axis = 1 ;
if isfield(opt,'showaxis')
    if strcmpi(opt.showaxis,'no'), show_axis=0; end
end

center_axis=1;
if isfield(opt,'showzeroaxis')
    if strcmpi(opt.showzeroaxis,'no'), center_axis=0; end
end

center_axis_x_ticks = 1;
center_axis_y_ticks = 0;
if center_axis
  if isfield(opt,'showzeroaxisxticks')
    if strcmpi(opt.showzeroaxisxticks,'no'), center_axis_x_ticks = 0; end
  end
  if isfield(opt,'showzeroaxisyticks')
    if strcmpi(opt.showzeroaxisyticks,'yes'), center_axis_y_ticks = 1; end
  end
end

show_colorbar = 1;
if isfield(opt,'showcolorbar')
    if strcmpi(opt.showcolorbar,'no'), show_colorbar = 0; end
end

if isfield(opt,'conditionnames')
 condition_names = opt.conditionnames;
else
 condition_names = {};
end

if isfield(opt,'footnotes')
  foot_notes = opt.footnotes;
else
  foot_notes = '';
end

if isfield(opt,'figname')
 fig_name = opt.figname;
 [fig_path, fig_name, fig_ext, versn] = fileparts(fig_name);
 fig_name = [fig_path fig_name];
else
 fig_name = {};
end

if isfield(opt,'title')
 title_string = opt.title;
else
 title_string = '';
end

if isfield(opt,'resolution') && ~isempty(opt.resolution)
  resolution = num2str(opt.resolution);
else
  resolution = '300';
end

save_figs = 0;
if isfield(opt,'savefigs');
  opt.savefigs = lower(opt.savefigs);
  if ~strcmp(opt.savefigs,'no')
    save_figs = 1;
    switch opt.savefigs
      case 'bmp'
        print_command{1} = '-dbmp';
      case 'emf'
        print_command{1} = '-dmeta';
        print_command{2} = sprintf('-r%s',resolution);
      case 'eps'
        print_command{1} = '-depsc';
        print_command{2} = '-tiff';
        print_command{3} = sprintf('-r%s',resolution);
      case 'hdf'
        print_command{1} = '-dhdf';
        print_command{2} = sprintf('-r%s',resolution);
      case 'ill'
        print_command{1} = '-dill';
        print_command{2} = sprintf('-r%s',resolution);
        opt.paperorient  = 'portrait';
      case 'jpg'
        print_command{1} = '-djpeg';
        print_command{2} = sprintf('-r%s',resolution);
      case 'pbm'
        print_command{1} = '-dpbm';
      case 'pcx'
        print_command{1} = '-dpcx24b';
      case 'pdf'
        print_command{1} = '-dpdf';
      case 'png'
        print_command{1} = '-dpng';
        print_command{2} = sprintf('-r%s',resolution);
      case 'ppm'
        print_command{1} = '-dppm';
      case 'tif'
        print_command{1} = '-dtiff';
        print_command{2} = sprintf('-r%s',resolution);
      otherwise
        fprintf('Using default EPS format to save figures.\n');
        opt.savefigs = 'eps';
        print_command{1} = '-depsc';
        print_command{2} = sprintf('-tiff');
        print_command{3} = sprintf('-r%s',resolution);
    end
  end
end

clear_figs = 0;
if isfield(opt,'clearfigs')
 if strcmpi(opt.clearfigs,'yes'), clear_figs = 1; end
end

autolayout = 1;
if isfield(opt,'layoutfile') && ~isempty(opt.layoutfile)
  if exist(opt.layoutfile,'file')
    autolayout = 0;
  else
    fprintf('Could not find specified layout file: %s.\n',opt.layoutfile);
  end
end

laminar = 0;
channamepos   = [.5 1];
channamealign = 'center';
% if isfield(opt,'laminar') && strcmpi(opt.laminar,'yes')
%     laminar = 1;
%     channamepos = [0 .5];
%     channamealign = 'right';
% end

%% Hidden Options

if isfield(opt,'max_rows')                                                 % maximum rows of subplots
 max_rows = opt.max_rows;
else
 max_rows = 8;
end
 
if isfield(opt,'max_columns')                                              % maximum columns of subplots
 max_columns = opt.max_columns;
else
 max_columns = 8;
end

if isfield(opt, 'max_group_size')                                          % if a chan group is larger than this it gets its own page
 max_group_size = opt.max_group_size;
else
 max_group_size = 30;
end

if isfield(opt,'plot_colors')                                             % list of colors to toggle through for condition averages
 plot_colors = opt.plot_colors;
else
 plot_colors = 'brgcmwky';
end

if isfield(opt,'plotlinewidth')
  plot_line_width = opt.plotlinewidth;
else
  plot_line_width = 0.5;
end

if isfield(opt,'statslinewidth')
  stats_line_width = opt.statslinewidth;
else
 if plot_line_width < 1
  stats_line_width = 0.1;
 else
  stats_line_width = 0.1;
 end
end

if isfield(opt,'statsmarker')
  stats_marker = opt.statsmarker;
else
  stats_marker = '*';
end

if isfield(opt,'markersize');
  marker_size = opt.markersize;
else
  marker_size = 2;
%   marker_size = 3;
end

if isfield(opt,'figbackground')                                           % background color of figure
 fig_background = opt.figbackground;
else
 fig_background = 'w';
end

if isfield(opt,'zeroaxiswidth')                                           % width of zero axis line if used
 zero_line_width = opt.zeroaxiswidth;
else
 zero_line_width = 0.25;                                                          
end

if isfield(opt,'zeroaxisstyle')
 line_style = opt.zeroaxisstyle;                                          % style of zero axis line
else
 line_style = '-';                    
end

if isfield(opt,'axisfontsize')                                            % size of axis font in points
 axis_font = opt.axisfontsize;
else
 axis_font  = 4;      
end

if isfield(opt,'papertype')                                               % paper size when saving figure
 papertype = opt.papertype;
else
 papertype = 'usletter';
end

if isfield(opt,'leftmargin')                                              % margins for paper
 leftmargin = opt.leftmargin;
else
 leftmargin = 0.5;
end

if isfield(opt,'rightmargin')
 rightmargin = opt.rightmargin;
else
 rightmargin = 0.5;
end

if isfield(opt,'topmargin')
 topmargin = opt.topmargin;
else
 topmargin   = 0.5;
end

if isfield(opt,'bottommargin')
 bottommargin = opt.bottommargin;
else
 bottommargin = 0.5;
end

if isfield(opt,'paperorient')                                              % printing landscape
 paperorient = opt.paperorient;
else
 paperorient = 'portrait';
end

if isfield(opt,'shading')                                                  % type of shading for time frequency plots
 shading_style = opt.shading;                                              % is a pcolor option
else
 shading_style = 'interp';
end

if isfield(opt,'whichfigs')                                                % if there is more than one page of figures which one to plot
    whichfigs = opt.whichfigs;
else
    whichfigs = [];
end

if isfield(opt,'rowpad')                                                   % these control how close plots are to each other           
  rowpad = opt.rowpad;
else
  rowpad = 0.01;
end

if isfield(opt,'columnpad')
  columnpad = opt.columnpad;
else
  columnpad = 0.01;
end

if isfield(opt,'showonlyfirstaxis');
  showonlyfirstaxis = opt.showonlyfirstaxis;
  if strcmpi(showonlyfirstaxis,'no')
    rowpad = rowpad + 0.02;
    columnpad = columnpad  + 0.02;
  end
else
  showonlyfirstaxis = 'yes';
end

invalid_stats = 0;

%% Set up for plotting.

% Load a Layout or Montage File

if ~autolayout
  [dum1,dum2,ext,dum3] = fileparts(opt.layoutfile);
  clear dum1 dum2 dum3
  ext = lower(ext);
  switch (ext)
    case '.lay'
      figs = ts_iEEG_readFTlayout(opt.layoutfile,data.sensor_info);
    case '.asc'
      figs = ts_iEEG_readNSmontage(opt.layoutfile,data.sensor_info);
    otherwise
      fprintf('The layout file must be a valid NeuroScan .asc montage file or a FieldTrip .lay file.\n');
      autolayout = 1;
  end
end

% Use sensor info to determine number and layout of figures.
if autolayout
  if length(channels) == 1                                                   % set up for single plot
    figs.plots.row   = 1;
    figs.plots.column= 1;
    [figs.plots.name,figs.plots.index] = intersect({data.sensor_info.label},channels);
    figs.name        = channels{1};                                        % labels figure with channel name as opposed to page #
  else                                                                       % set up plotting for multiple channels
    if ~isempty(channels)
      ignore_channels = [setdiff({data.sensor_info.label},channels)...     % if given list of channels set up ignore_channels to those not given
        ignore_channels];                                                  % and preserve any bad channels
    end
    figs = ts_iEEG_makefigs (data.sensor_info,max_rows,max_columns,max_group_size,ignore_channels,laminar);
    if ~isempty(channels), try figs.name = cell2mat(channels); end; end
    if ~isempty(whichfigs), figs = figs(whichfigs); end
  end

  % Set up some constants for subplot command


  for i = 1:length(figs)
    figs(i).max_rows    = 0;
    figs(i).max_columns = 0;
    for j=1:length(figs(i).plots)
      figs(i).max_rows    = max(figs(i).max_rows,figs(i).plots(j).row);       % adjust size based on figure
      figs(i).max_columns = max(figs(i).max_columns,figs(i).plots(j).column);
    end

    figs(i).width = ((.85-(figs(i).max_columns*columnpad))/figs(i).max_columns);  % allowable width for each plot
    % this accounts for padding between plots
    figs(i).height = ((.9-(figs(i).max_rows*rowpad))/figs(i).max_rows);           % allowable height for each plot with a
    % pad between rows
    
    % Prevent wierd aspect ratios.
    if figs(i).width > (1.75*figs(i).height)
      figs(i).width  = figs(i).width * .75;
    end
    if figs(i).height > (1.75*figs(i).width)
      figs(i).height = figs(i).height *.75;
    end

    figs(i).ypos = [];
      for j = 1:figs(i).max_rows
        figs(i).ypos = cat(1,figs(i).ypos,(0.9-(j*figs(i).height))-(rowpad*(j-1)));
      end                                                                           
    figs(i).ypos = repmat(figs(i).ypos,[1 figs(i).max_columns]);

    figs(i).xpos = [];
    for j = 1:figs(i).max_columns
      figs(i).xpos=cat(2,figs(i).xpos,((j-1)*figs(i).width)+(columnpad*(j-1))+.05); % where to place each plot on the row
    end
    figs(i).xpos = repmat(figs(i).xpos,[figs(i).max_rows 1]);

    figs(i).width  = repmat(figs(i).width,[figs(i).max_rows figs(i).max_columns]);
    figs(i).height = repmat(figs(i).height,[figs(i).max_rows figs(i).max_columns]);
  end
end

%% Plotting



scrsz = get(0,'ScreenSize');                                               % initial figure size based on screen size
plot_colors = strrep(plot_colors,fig_background,'');                       % remove background color from plot colors
bctag ='';
        if baseline
            fprintf('subtracting baseline...\n');
            if 1
              if isempty(baseline_end),   baseline_end = 0;      end
              if isempty(baseline_start), baseline_start = -inf; end
              tot   = 0;
              ntrls = [];
              for j = 1:length(conds)
                t     = data.(datafield)(cond(j)).time;
                tix   = find(t >= baseline_start & t <= baseline_end);
                ntrls = [ntrls data.(datafield)(cond(j)).num_trials];
                tot   = tot + data.(datafield)(cond(j)).data(:,tix)*ntrls(end);
              end
              tot = tot / sum(ntrls);
              mean_baseline = mean(tot,2);
              clear tot
              for j = 1:length(conds)
                data.(datafield)(conds(j)).data = data.(datafield)(conds(j)).data - mean_baseline*ones(1,length(t));
              end
            else
              for j=1:length(conds)
                  bctag='_bc';
                  num_samples = length(data.(datafield)(1).time);
                  if isempty(baseline_start), baseline_start = data.(datafield)(conds(j)).time(1); end
                  if isempty(baseline_end),   baseline_end   = 0; end
                  samples = find(data.(datafield)(conds(j)).time >= baseline_start & data.(datafield)(conds(j)).time <= baseline_end);
                  mean_baseline = ...
                      mean(data.(datafield)(conds(j)).data(:,samples(1):samples(end)),2);
                  data.(datafield)(conds(j)).data = ...
                      data.(datafield)(conds(j)).data - mean_baseline*ones(1,num_samples);
              end
            end
        end
        for i = 1:length(figs)
            if calc_ylim
                for j = 1:length(figs(i).plots)
                  if data.sensor_info(figs(i).plots(j).index).badchan, continue; end
                    for k=1:length(conds)
                        xmin_idx = nearest(data.(datafield)(conds(k)).time,xmin);
                        xmax_idx = nearest(data.(datafield)(conds(k)).time,xmax);
                        ymin = min([ymin data.(datafield)(conds(k)).data(figs(i).plots(j).index,xmin_idx:xmax_idx)]);% look for y min/max across channels from this page
                        ymax = max([ymax data.(datafield)(conds(k)).data(figs(i).plots(j).index,xmin_idx:xmax_idx)]);% and all conditions
                    end
                end
            end
            limits_tag = sprintf('_sec%s-%s',num2str(xmin),num2str(xmax));
            if ~isempty(fig_name)
                if plot_stats, this_figures_name = sprintf('%s_stats_%s_%s',fig_name,...
                              opt.stat_data.stats(1).parms.method,opt.stat_data.stats(1).parms.statistic); 
                else this_figures_name = sprintf('%s',fig_name); end
                this_figures_name = sprintf('%s%s%s_%s',this_figures_name,limits_tag,bctag,figs(i).name);
            else
                this_figures_name = 'avg_plot';
                for j=i:length(conds)
                    this_figures_name = [this_figures_name '_cond' num2str(conds(j))];
                end
                if plot_stats, this_figures_name = sprintf('%s_stats_%s_%s',this_figures_name,...
                              opt.stat_data.stats(1).parms.method,opt.stat_data.stats(1).parms.statistic); end
                this_figures_name = sprintf('%s%s%s_%s',this_figures_name,limits_tag,bctag,figs(i).name);
            end         
            figure('Name',this_figures_name,...                                    % set up figure window
                'NumberTitle','off',...
                'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],...
                'Color',fig_background);                                           % background color            
            leg_ypos  = .85;
            leg_xpos  = .90;
            for j = 1:length(figs(i).plots)
                subplot('Position',...
                  [figs(i).xpos(figs(i).plots(j).row,figs(i).plots(j).column) ...
                   figs(i).ypos(figs(i).plots(j).row,figs(i).plots(j).column)...
                   figs(i).width(figs(i).plots(j).row,figs(i).plots(j).column) ...
                   figs(i).height(figs(i).plots(j).row,figs(i).plots(j).column)]);
                hold on                                                            % only add to plot, make sure nothing is replaced
                curr_color = 1;
                for k=1:length(conds)
                  if data.sensor_info(figs(i).plots(j).index).badchan ~= 1
                    plot_data = data.(datafield)(conds(k)).data(figs(i).plots(j).index,data.(datafield)(conds(k)).time >= xmin & data.(datafield)(conds(k)).time <= xmax);
                    plot_time = data.(datafield)(conds(k)).time(data.(datafield)(conds(k)).time >= xmin & data.(datafield)(conds(k)).time <= xmax);
                    plot(plot_time,...                            % time
                         plot_data,...   % data for current channel
                         plot_colors(curr_color),...                               % set color of plot
                         'linewidth',plot_line_width,...
                         'Clipping','off');   
                    if plot_stats                                                    % statistics plotting
                      stat_index = find([opt.stat_data.stats.event_code] == data.(datafield)(conds(k)).event_code);               % make sure there are stats for that condition
                      if ~isempty(stat_index) 
                          if ismember(data.sensor_info(figs(i).plots(j).index).label,{opt.stat_data.sensor_info.label})
                            idx  = find(ismember({opt.stat_data.sensor_info.label},data.sensor_info(figs(i).plots(j).index).label));
                            sigs = find(opt.stat_data.stats(stat_index).mask(idx,:) == 1);
                            ttd  = data.(datafield)(conds(k)).time;
                            tts  = opt.stat_data.stats(stat_index).time(sigs);                            
%                             [dum,s2dt] = intersect(ttd,tts);
%                             if isempty(s2dt)
                              s2dt = [];
                              for tix = 1:length(tts)
                                s2dt = [s2dt nearest(ttd,tts(tix))];
                              end
%                             end
                            sigts = opt.stat_data.stats(stat_index).time(sigs);
                            sigtd = data.(datafield)(conds(k)).data(figs(i).plots(j).index,s2dt);
%                             if length(sigts) ~= length(sigtd)
%                               sigtd = [];
%                               for kk = 1:length(sigts)
%                                 sigtd = [sigtd nearest(data.(datafield)(conds(k)).time,sigts(kk))];
%                               end
%                               sigtd = unique(sigtd);
%                             end
                            if length(sigts) ~= length(s2dt)
                              s2dt = [];
                              for kk = 1:length(sigts)
                                s2dt = [s2dt nearest(data.(datafield)(conds(k)).time,sigts(kk))];
                              end
                              s2dt = unique(s2dt);
                              sigtd = data.(datafield)(conds(k)).data(figs(i).plots(j).index,s2dt);
                            end
                            scatter(sigts,sigtd,...
                                    'SizeData',marker_size,...
                                    'MarkerEdgeColor',plot_colors(curr_color),...
                                    'MarkerFaceColor',plot_colors(curr_color),...
                                    'Marker',stats_marker,...
                                    'Clipping','off');
                          end
                      else
                        if j == 1
                          fprintf('WARNING: Could not locate event %f in the statistics data provided.\n',data.(datafield)(conds(k)).event_code);
                          invalid_stats = 1;
                        end
                      end  % isempty(stat_index)                      
                    end    % plot_stats
                    clear plot_time plot_data
                  end % skip plotting for bad channels
                    if  leg_ypos < 0.07, leg_xpos =  0;                            % when more then max_rows conditions to plot
                        leg_ypos = .9;  end                       % place other legends on right side                      
                    if j==length(figs(i).plots)
                        if ~isempty(condition_names)
                            annotation('textbox','Position',[leg_xpos leg_ypos .1 .1],...
                                'VerticalAlignment','top',...
                                'HorizontalAlignment','left',...
                                'FontSize',8,...
                                'Color',plot_colors(curr_color),...
                                'FitHeightToText','on',...
                                'LineStyle','none',...
                                'String',...
                                sprintf('%s\nNum Trials: %s', condition_names{k},num2str(data.(datafield)(conds(k)).num_trials))); 
                        else                                                        % make up a legend
                            annotation('textbox','Position',[leg_xpos leg_ypos .1 .1],...
                                'VerticalAlignment','top',...
                                'HorizontalAlignment','left',...
                                'FontSize',8,...
                                'Color',plot_colors(curr_color),...
                                'FitHeightToText','on',...
                                'LineStyle','none',...
                                'String',sprintf('Condition %s\nNum Trials: %s',num2str(conds(k)),num2str(data.(datafield)(conds(k)).num_trials)));
                        end                       
                     leg_ypos  = leg_ypos - 0.05;                                   % move y position for each legend down 
                    end
                    curr_color = curr_color + 1;
                    if curr_color > length(plot_colors), curr_color = 1; end       % cycle through the 8 colors over and over
                end % for loop of conditions                
                xrange = xmin:xmax;
                title(figs(i).plots(j).name,'FontSize',8,...     
                    'VerticalAlignment','Middle','Clipping','off');                                             % show channel name
                box off;                                                           % remove MATLAB box from plot
                try xlim([xmin xmax]); end                                                 % set up axis could also use:
                try ylim([ymin ymax]); end                                                % axis([xmin xmax ymin ymax]);
                if show_axis && data.sensor_info(figs(i).plots(j).index).badchan ~= 1 % display MATLAB axis
                  if strcmpi(showonlyfirstaxis,'no')
                    axis on;
                    set(gca,'FontSize',axis_font);
                    x_tick_markers = get(gca,'xtick');
                    y_tick_markers = get(gca,'ytick');
                  else
                   if j==1
                    axis on;
                    set(gca,'FontSize',axis_font);
                    x_tick_markers = get(gca,'xtick');
                    y_tick_markers = get(gca,'ytick');
                   else
                     axis off;
                   end
                  end
                else
                    axis off;
                    x_tick_markers = get(gca,'xtick');
                    y_tick_markers = get(gca,'ytick');
                end
                if center_axis && data.sensor_info(figs(i).plots(j).index).badchan ~= 1    % display manually created zero axis
                    plot(xlim,zeros(2),line_style,'color','black','linewidth',zero_line_width,'Clipping','on');
                    if center_axis_x_ticks, scatter(x_tick_markers,zeros(1,length(x_tick_markers)),'+k','linewidth',zero_line_width,'Clipping','on'); end   % tick marks on zero axis
                    plot(zeros(2),ylim,line_style,'color','black','linewidth',zero_line_width,'Clipping','on');
                    if center_axis_y_ticks, scatter(zeros(1,length(y_tick_markers)),y_tick_markers,'+k','linewidth',zero_line_width,'Clipping','on'); end
                end
                if plot_stats && data.sensor_info(figs(i).plots(j).index).badchan ~= 1
                  stat_index = find([opt.stat_data.stats.event_code] == -1);               % make sure there are stats for that condition
                  if ~isempty(stat_index) %&& ~(invalid_stats)
                    if ismember(data.sensor_info(figs(i).plots(j).index).label,{opt.stat_data.sensor_info.label})
                      idx = find(ismember({opt.stat_data.sensor_info.label},data.sensor_info(figs(i).plots(j).index).label));
                    
                      sigs = find(opt.stat_data.stats(stat_index).mask(idx,:)==1);
                      % shade regions of significant differences
                      mask = squeeze(opt.stat_data.stats(stat_index).mask(idx,:)) > 0;
                      t    = opt.stat_data.stats(stat_index).time;
                      n = length(mask);
                      mask  = [0 mask(2:end-1).*(((mask(1:n-2)==0) + (mask(3:n)==0)) > 0) 0];
                      edges = find(mask);
                      for k = 1:length(edges)/2
                        ix1 = t(edges(2*k-1));
                        ix2 = t(edges(2*k));
                        iy1 = ymin;
                        iy2 = ymax;
                        fill([ix1 ix1 ix2 ix2],[iy1 iy2 iy2 iy1],'b','EdgeColor','k',...
                          'FaceAlpha',.1,'EdgeAlpha',0);%.1);
                      end                    
  %                     scatter(opt.stat_data.stats(stat_index).time(sigs),zeros(1,length(opt.stat_data.stats(stat_index).time(sigs))),...
  %                                   'SizeData',marker_size,...
  %                                   'MarkerEdgeColor','k',...
  %                                   'MarkerFaceColor','k',...
  %                                   'Marker',stats_marker);
                    end
                  else
                    if j == 1
                      fprintf('WARNING: Could not find the statistics between the plotted conditions in the given statistics data.\n');
                    end
                  end  % isempty(stat_index)
                end
                axis square
            end %% each channel
            if calc_ylim, ymin = []; ymax = []; end                                % reset y lim for new page
            annotation('textbox','Position',[0 .9 1 .1],...                        % set up figure's title box
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'Color','k',...
                'FontSize',10,...
                'FontWeight','bold',...
                'FitHeightToText','on',...
                'LineStyle','none',...
                'String',title_string);
            annotation('textbox','Position',[leg_xpos leg_ypos .1 .1],...          % set up figure's footnotes
                'VerticalAlignment','top',...
                'HorizontalAlignment','center',...
                'Color','k',...
                'FontSize',6,...
                'FitHeightToText','on',...
                'LineStyle','none',...
                'String',foot_notes);
            hold off
            if save_figs
                set(gcf,'PaperOrientation',paperorient,...                         % set paper properties first in order to get
                    'PaperType',papertype,...                                      % proper paper size for figure
                    'PaperUnits','inches','PaperPositionMode','auto');
                papersize = get(gcf,'PaperSize');
                fig_width = papersize(1) - leftmargin - rightmargin;               % fill page
                fig_height= papersize(2) - topmargin - bottommargin;
                fig_left  = (papersize(1)-fig_width)/2;                            % center figure
                fig_bottom= (papersize(2)-fig_height)/2;
                set(gcf,'PaperPosition', [fig_left, fig_bottom, fig_width, fig_height]);
                this_figures_name = [this_figures_name '.' opt.savefigs];
                if isfield(opt,'figpath') && ~isempty(opt.figpath), this_figures_dir = opt.figpath;
                else this_figures_dir  = fullfile(pwd,'images');
                end
                if ~exist(this_figures_dir,'dir'), mkdir(this_figures_dir); end
                this_figures_name = fullfile(this_figures_dir,this_figures_name);
                fprintf('Saving figure: %s\n',this_figures_name);
                print(print_command{:}, this_figures_name);                %%% print to an eps file
            end
            if clear_figs, close; end                                              % don't keep figures on screen option
        end % each page